"""Timing harness: python bench.py <label>; run with PYTHONPATH set to the
symox python tree under test. Prints JSON lines: case, first (fresh process
is handled by the driver), median of repeated calls, result string."""
import json, sys, time, statistics
t0 = time.perf_counter()
import symox
from symox import Symbol, symbols, cos, pi, integrate, Integral
try:
    from symox.utilities.cache import clear_cache
except ImportError:
    clear_cache = lambda: None
x, y, n = symbols('x y n'); th = Symbol('theta')
CASES = {
    'cos3_def': lambda: integrate(cos(3*th)/(5 - 4*cos(th)), (th, 0, 2*pi)),
    'cos2_def': lambda: integrate(cos(2*th)/(5 - 4*cos(th)), (th, 0, 2*pi)),
    'cos2_indef': lambda: integrate(cos(2*th)/(5 - 4*cos(th)), th),
    'cos1_def': lambda: integrate(cos(th)/(5 - 4*cos(th)), (th, 0, 2*pi)),
    'inv2cos_final_false': lambda: Integral(1/(2 + cos(x)), x).doit(final=False),
    'x_pow_ny': lambda: integrate(x**(n*y), y),
    'two_pow_ny': lambda: integrate(2**(n*y), y),
    'three_x_pow_ny': lambda: integrate(3*x**(n*y), y),
    'x_pow_ny_plus3': lambda: integrate(x**(n*y + 3), y),
    'x_pow_2y': lambda: integrate(x**(2*y + 1), y),
}
case = sys.argv[2]
f = CASES[case]
s = time.perf_counter(); r = f(); first = time.perf_counter() - s
reps = []
for _ in range(5):
    clear_cache()
    s = time.perf_counter(); f(); reps.append(time.perf_counter() - s)
print(json.dumps(dict(label=sys.argv[1], case=case, first=first, median=statistics.median(reps), result=str(r)[:200])))
