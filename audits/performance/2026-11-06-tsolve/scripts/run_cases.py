"""Replay the transcendental-fallback corpus with one backend.

usage: run_cases.py {symox,sympy} OUT.json [--timeout S] [--jobs N]
                    [--corpus tsolve|unrad] [--trace]

``--corpus unrad`` replays the radical corpus of the unrad lot
(``2026-11-06-unrad/scripts/cases.py``). ``--trace`` (symox only) records,
per case, which Python fallbacks of ``solve`` ran: ``_masked_solution``,
``_inverted_solution``, ``_unrad_solution`` (``_solve_relations``),
``_solve_piecewise`` and ``tsolve_step`` (``symox.solvers._tsolve_steps``),
and ``_unrad_resultant``, and the Python computations the native host calls
back (``nsimplify`` for ``mpmath.identify``, ``equals``).

The backend is imported once and every case runs in a forked child (with a
wall-clock timeout); the result is recorded as ``str`` of the solution list
(or of the exception) together with its ``srepr`` for the comparison.
"""
import json
import os
import select
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
if '--corpus' in sys.argv and sys.argv[sys.argv.index('--corpus') + 1] == 'unrad':
    sys.path.insert(0, os.path.join(HERE, '..', '..', '2026-11-06-unrad', 'scripts'))
else:
    sys.path.insert(0, HERE)
from cases import CASES, POSITIVE, SYMBOLS  # noqa: E402

TRACE = []


def install_trace():
    import symox
    import symox._solve_relations as rel
    import symox.solvers._tsolve_steps as steps
    import symox.solvers._algebraic as alg

    def spy(module, name):
        original = getattr(module, name)

        def wrapper(*args, **kwargs):
            TRACE.append(name)
            return original(*args, **kwargs)
        setattr(module, name, wrapper)
    for name in ('_masked_solution', '_inverted_solution', '_unrad_solution'):
        spy(rel, name)
    spy(steps, 'tsolve_step')
    spy(alg, '_unrad_resultant')
    if hasattr(symox, '_solve_piecewise'):
        spy(symox, '_solve_piecewise')
    # the computations that the native `_tsolve` host takes from Python
    # (`mpmath.identify` through `nsimplify`, the tail of `Expr.equals`)
    spy(symox, 'nsimplify')
    import symox.core.expr as cexpr
    for cls in cexpr.Expr.__mro__:
        if 'equals' in cls.__dict__:
            original = cls.__dict__['equals']

            def equals(self, *args, _original=original, **kwargs):
                TRACE.append('equals')
                return _original(self, *args, **kwargs)
            setattr(cls, 'equals', equals)
            break


def namespace(backend):
    if backend == 'symox':
        import symox as mod
    else:
        import sympy as mod
    ns = {k: getattr(mod, k) for k in dir(mod) if not k.startswith('_')}
    for name in SYMBOLS.split():
        ns[name] = mod.Symbol(name)
    for name in POSITIVE.split():
        ns[name] = mod.Symbol(name, positive=True)
    return mod, ns


def run_one(mod, ns, expression, symbols, flags):
    eq = eval(expression, ns)
    sym = eval(symbols, ns)
    start = time.time()
    try:
        r = mod.solve(eq, sym, **flags)
        out = {'ok': True, 'str': str(r), 'srepr': mod.srepr(r) if not isinstance(r, (set, tuple)) else str(r)}
    except Exception as e:  # the error class is part of the observation
        out = {'ok': False, 'error': type(e).__name__, 'str': str(e)[:300]}
    out['seconds'] = round(time.time() - start, 3)
    if TRACE:
        out['python_fallbacks'] = sorted(set(TRACE))
    return out


def main():
    backend, outpath = sys.argv[1], sys.argv[2]
    timeout = 60.0
    jobs = 8
    args = sys.argv[3:]
    if '--timeout' in args:
        timeout = float(args[args.index('--timeout') + 1])
    if '--jobs' in args:
        jobs = int(args[args.index('--jobs') + 1])
    mod, ns = namespace(backend)
    if '--trace' in args and backend == 'symox':
        install_trace()
    results = {}
    pending = list(CASES)
    running = {}
    while pending or running:
        while pending and len(running) < jobs:
            cid, expression, symbols, flags = pending.pop(0)
            r, w = os.pipe()
            pid = os.fork()
            if pid == 0:
                os.close(r)
                try:
                    out = run_one(mod, ns, expression, symbols, flags)
                except Exception as e:
                    out = {'ok': False, 'error': 'Setup' + type(e).__name__, 'str': str(e)[:300]}
                with os.fdopen(w, 'w') as f:
                    f.write(json.dumps(out))
                os._exit(0)
            os.close(w)
            running[r] = (pid, cid, time.time(), [])
        ready, _, _ = select.select(list(running), [], [], 0.2)
        for r in ready:
            data = os.read(r, 1 << 20)
            if data:
                running[r][3].append(data)
                continue
            pid, cid, _, chunks = running.pop(r)
            os.close(r)
            os.waitpid(pid, 0)
            try:
                results[cid] = json.loads(b''.join(chunks).decode())
            except ValueError:
                results[cid] = {'ok': False, 'error': 'Crash', 'str': ''}
        now = time.time()
        for r, (pid, cid, started, _) in list(running.items()):
            if now - started > timeout:
                os.kill(pid, 9)
                os.waitpid(pid, 0)
                os.close(r)
                running.pop(r)
                results[cid] = {'ok': False, 'error': 'Timeout', 'str': '', 'seconds': timeout}
    ordered = {cid: results[cid] for cid, *_ in CASES}
    with open(outpath, 'w') as f:
        json.dump(ordered, f, indent=1, sort_keys=False)


if __name__ == '__main__':
    main()
