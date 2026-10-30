"""Time a few public checksol calls: one cold call, then the median of 5.
usage: bench.py sympy|symox out.json"""
import json
import statistics
import sys
import time
import warnings

lib = sys.argv[1]
M = __import__(lib)
checksol = __import__(lib + '.solvers.solvers', fromlist=['checksol']).checksol
x, y, z = M.symbols('x y z')
f = M.Function('f')
CALLS = [
    ('x**2-4 at 2', lambda: checksol(x**2 - 4, x, 2)),
    ('Eq(x**2,4) at -2', lambda: checksol(M.Eq(x**2, 4), x, -2)),
    ('system 2x2', lambda: checksol([x + y - 3, 2*x - y], {x: 1, y: 2})),
    ('exp(x)-y at log(y)', lambda: checksol(M.exp(x) - y, x, M.log(y))),
    ('x*exp(x)-1 at LambertW(1)', lambda: checksol(x*M.exp(x) - 1, x, M.LambertW(1))),
    ('trig identity sol', lambda: checksol(x - 2, {x: 2*M.sin(y)**2 + 2*M.cos(y)**2})),
    ('sqrt(x)+1 at 1', lambda: checksol(M.sqrt(x) + 1, x, 1)),
    ('z**x-y minimal', lambda: checksol(z**x - y, {x: M.log(y)/M.log(z)}, minimal=True)),
    ('Derivative f(x)', lambda: checksol(M.Derivative(f(x), x) - 1, {f(x): x})),
    ('unverified warn', lambda: checksol(x*M.exp(x) - y, x, M.LambertW(y), warn=True)),
]
out = {}
for name, call in CALLS:
    M.core.cache.clear_cache() if hasattr(M.core, 'cache') else None
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        t0 = time.perf_counter()
        r = call()
        cold = time.perf_counter() - t0
        warm = []
        for _ in range(5):
            t0 = time.perf_counter()
            call()
            warm.append(time.perf_counter() - t0)
    out[name] = {'r': r, 'cold_ms': cold * 1e3, 'warm_ms': statistics.median(warm) * 1e3}
json.dump(out, open(sys.argv[2], 'w'), indent=1)
for k, v in out.items():
    print('%-28s %-6s %8.2f %8.2f' % (k, v['r'], v['cold_ms'], v['warm_ms']))
