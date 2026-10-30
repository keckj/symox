"""Replay recorded base-case calls and store results and timings.

usage: replay_calls.py <sympy|symox|symox-stock> CALLS.json OUT.json
symox-stock runs the Python bodies (native entry points disabled).
"""
import json
import signal
import sys
import time
import importlib

backend, src, out = sys.argv[1:4]
S = importlib.import_module('sympy' if backend == 'sympy' else 'symox')
prde = importlib.import_module(S.__name__ + '.integrals.prde')
risch = importlib.import_module(S.__name__ + '.integrals.risch')
if backend == 'symox-stock':
    from symox.integrals import _prde_native as nat
    for name in ('is_log_deriv_k_t_radical_in_field', 'parametric_log_deriv_heu', 'limited_integrate'):
        setattr(nat, name, lambda *a, **k: None)

x = S.Symbol('x')
DE = risch.DifferentialExtension(extension={'D': [S.Poly(1, x)]})


def dec(p):
    return S.Poly(S.sympify(p[0], locals={'x': x}), x, domain=p[1])


def enc(v):
    if v is None:
        return None
    if isinstance(v, S.Poly):
        return ['Poly', S.srepr(v.as_expr()), str(v.get_domain())]
    if isinstance(v, (tuple, list)):
        return [type(v).__name__] + [enc(i) for i in v]
    return [type(v).__name__, S.srepr(v)]


def call(kind, args):
    if kind == 'is_log_deriv':
        return lambda: prde.is_log_deriv_k_t_radical_in_field(dec(args[0]), dec(args[1]), DE)
    if kind == 'heu':
        return lambda: prde.parametric_log_deriv_heu(*[dec(a) for a in args], DE)
    G = [(dec(a), dec(d)) for a, d in args[2]]
    return lambda: prde.limited_integrate(dec(args[0]), dec(args[1]), list(G), DE)


def _alarm(*a):
    raise TimeoutError


signal.signal(signal.SIGALRM, _alarm)
data = json.load(open(src))
results = []
for kind, args in data['records']:
    f = call(kind, args)
    best = None
    res = None
    for rep in range(3):
        signal.alarm(20)
        t0 = time.perf_counter()
        try:
            r = enc(f())
        except TimeoutError:
            r = ['EXC', 'Timeout']
        except Exception as exc:
            r = ['EXC', type(exc).__name__]
        dt = time.perf_counter() - t0
        signal.alarm(0)
        best = dt if best is None else min(best, dt)
        res = r
        if r == ['EXC', 'Timeout']:
            break
    results.append([kind, res, best])
stats = None
if backend == 'symox':
    from symox.integrals import _prde_native as nat
    stats = dict(nat.STATS)
json.dump({'backend': backend, 'results': results, 'stats': stats}, open(out, 'w'))
print(backend, 'replayed', len(results), 'calls; total %.3f s' % sum(r[2] for r in results))
