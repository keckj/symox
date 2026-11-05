"""Run the odesolve cases with sympy or symox and write one JSON line per case.

usage: runner.py sympy|symox|symox-stock OUT.jsonl [name-prefix]

`symox-stock` removes the native particular-solution kernels
(`ode_uc_trialset`, `ode_uc_solve`, `ode_vop_solve`) from `_core`, so the
reference bodies of `nonhomogeneous` run.
For each case: dsolve(eq, f(x), hint=h), dsolve(..., hint=h, simplify=False)
and dsolve(eq, f(x)) (default hint), as str, or 'EXC <type>: <msg>'.
"""
import json
import os
import signal
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen_eqs  # noqa: E402

lib = sys.argv[1]
if lib == 'sympy':
    import sympy as L
else:
    import symox as L
    if lib == 'symox-stock':
        # the reference bodies of the particular-solution helpers
        from symox import _core
        for _k in ('ode_uc_trialset', 'ode_uc_solve', 'ode_vop_solve'):
            if hasattr(_core, _k):
                delattr(_core, _k)

ns = {k: getattr(L, k) for k in dir(L) if not k.startswith('_')}
ns['f'] = L.Function('f')
ns['x'] = L.Symbol('x')
ns['a'], ns['b'] = L.symbols('a b')
ns['k'] = L.Symbol('k', real=True)
ns['p'] = L.Symbol('p', positive=True)


class TO(Exception):
    pass


def alarm(*a):
    raise TO()


signal.signal(signal.SIGALRM, alarm)


def call(fn, limit=30):
    signal.alarm(limit)
    t = time.perf_counter()
    try:
        r = str(fn())
    except TO:
        r = 'TIMEOUT'
    except Exception as e:  # noqa: BLE001
        r = 'EXC %s: %s' % (type(e).__name__, str(e)[:200])
    finally:
        signal.alarm(0)
    return r, time.perf_counter() - t


prefix = sys.argv[3] if len(sys.argv) > 3 else ''
# ODE_SKIP: comma separated case names (symox exceeds the memory cap inside
# integration on vop_016, before and after this change, and on vop_013 with
# simplify=False after the dotprodsimp fix; SymPy times out on both).
# The Rust allocation failure aborts the process, so such cases are skipped.
skip = set(filter(None, os.environ.get('ODE_SKIP', '').split(',')))
with open(sys.argv[2], 'w') as out:
    for name, src, hint in gen_eqs.cases():
        if not name.startswith(prefix) or name in skip:
            continue
        eq = eval(src, ns)
        fx = ns['f'](ns['x'])
        rec = {'name': name, 'eq': src, 'hint': hint}
        rec['hint_sol'], rec['t_hint'] = call(lambda: L.dsolve(eq, fx, hint=hint))
        rec['nosimp'], _ = call(lambda: L.dsolve(eq, fx, hint=hint, simplify=False))
        rec['default'], rec['t_def'] = call(lambda: L.dsolve(eq, fx))
        out.write(json.dumps(rec) + '\n')
        out.flush()
