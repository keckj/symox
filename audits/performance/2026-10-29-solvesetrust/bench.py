"""solveset timings: first call in a fresh process and the median of
repeated calls. Usage: bench.py {symox|sympy} OUT.json"""
import json
import statistics
import sys
import time

if sys.argv[1] == 'sympy':
    from sympy import Eq, Mod, S, Symbol, log, solveset, sqrt
else:
    from symox import Eq, Mod, S, Symbol, log, solveset, sqrt

x = Symbol('x')
CASES = [
    ('mod_linear', lambda: solveset(3 - Mod(5*x - 8, 7), x, S.Integers)),
    ('mod_square_160', lambda: solveset(Mod(x**2, 160) - 9, x, S.Integers)),
    ('mod_dlog_41', lambda: solveset(Mod(7**x, 41) - 15, x, S.Integers)),
    ('radical_sqrt', lambda: solveset(sqrt(x) - x + 2, x, S.Reals)),
    ('radical_two_roots', lambda: solveset(sqrt(x + 1) + sqrt(x) - 3, x, S.Reals)),
    ('eq_cancelled_pole', lambda: solveset(Eq(1/x + x, 1/x), x, S.Reals)),
    ('x_log_x', lambda: solveset(x*log(x), x, S.Reals)),
    ('identically_zero', lambda: solveset(-x**2 - 2*x + (x + 1)**2 - 1, x, S.Reals)),
    ('ineq_pole_reals', lambda: solveset(x + 1/x > -2 + 1/x, x, S.Reals)),
    ('ineq_pole_integers', lambda: solveset(x + 1/x > -2 + 1/x, x, S.Integers)),
]
REPEAT = 21
out = {}
for name, call in CASES:
    t = time.perf_counter()
    result = call()
    first = time.perf_counter() - t
    times = []
    for _ in range(REPEAT):
        t = time.perf_counter()
        call()
        times.append(time.perf_counter() - t)
    out[name] = {'first_s': first, 'median_s': statistics.median(times), 'result': str(result)}
json.dump(out, open(sys.argv[2], 'w'), indent=1)
