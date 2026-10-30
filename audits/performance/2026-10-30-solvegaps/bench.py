"""Ten solve/solveset/unrad calls of the solvegaps task.

usage: python bench.py sympy|symox out.json
For each call, the first (cold) time and the
median of 5 further (warm) runs are recorded, with str(result).
"""
import json
import statistics
import sys
import time

lib = sys.argv[1]
if lib == 'sympy':
    import sympy as S
    from sympy.solvers.solvers import unrad
else:
    import symox as S
    from symox.solvers.solvers import unrad

x, y = S.symbols('x y')
xr, yr = S.Symbol('x', real=True), S.Symbol('y', real=True)
R = S.Rational
CASES = [
    ('solveset 11534', lambda: S.solveset(-yr + xr/S.sqrt(1 - xr**2), xr, S.S.Reals)),
    ('solve x**(2*y)-3', lambda: S.solve(x**(2*y) - 3, x)),
    ('solve x**(2*y)-3*x**y+2', lambda: S.solve(x**(2*y) - 3*x**y + 2, x)),
    ('solve x**(3*y)-3', lambda: S.solve(x**(3*y) - 3, x)),
    ('solve sqrt(x)+sqrt(x+1)-3', lambda: S.solve(S.sqrt(x) + S.sqrt(x + 1) - 3, x)),
    ('solveset sqrt(x)+sqrt(x+1)-3', lambda: S.solveset(S.sqrt(x) + S.sqrt(x + 1) - 3, x, S.S.Reals)),
    ('solveset sqrt(x-1)-x+3', lambda: S.solveset(S.sqrt(x - 1) - x + 3, x, S.S.Reals)),
    ('solve x**(1/3)+sqrt(x)-1', lambda: S.solve(x**R(1, 3) + S.sqrt(x) - 1, x)),
    ('unrad cov', lambda: unrad(x**R(1, 3) + S.sqrt(x) - 1, x, cov=[])),
    ('solve x*exp(x)-1', lambda: S.solve(x*S.exp(x) - 1, x)),
]
out = []
for name, fn in CASES:
    t0 = time.perf_counter()
    try:
        r = str(fn())
    except Exception as e:  # recorded, not hidden
        r = 'error: %s: %s' % (type(e).__name__, e)
    cold = time.perf_counter() - t0
    warm = []
    for _ in range(5):
        t0 = time.perf_counter()
        try:
            fn()
        except Exception:
            pass
        warm.append(time.perf_counter() - t0)
    out.append(dict(name=name, result=r, cold=cold, warm=statistics.median(warm)))
    print('%-32s cold %8.4f s  warm %8.4f s  %s' % (name, cold, statistics.median(warm), r[:90]))
json.dump(out, open(sys.argv[2], 'w'), indent=1)
