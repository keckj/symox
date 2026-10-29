"""Before/after timing of public holonomic calls.

Usage: PYTHONPATH=<tree>/python python bench.py <label> <mode>
mode 'first': one call per case in this fresh process (run once per process);
mode 'repeat': median of 7 repeated calls after one warm-up.
"""
import sys, time, json, statistics
from symox import *
from symox.holonomic import expr_to_holonomic
x, a = symbols('x a')
CASES = {
    'sin_x0_0': lambda: expr_to_holonomic(sin(x), x),
    'sin_pow3': lambda: expr_to_holonomic(sin(x), x)**3,
    'exp_pow_neg': lambda: expr_to_holonomic(exp(x), x)**-2,
    'sqrt_x': lambda: expr_to_holonomic(sqrt(x), x),
    'x2_plus_x': lambda: expr_to_holonomic(x**2 + x, x),
    'sinc': lambda: expr_to_holonomic(sin(x)/x, x),
    'erfc': lambda: expr_to_holonomic(erfc(x), x),
    'erf_rewrite_erfc': lambda: erf(x).rewrite(erfc),
    'sin_plus_cos': lambda: expr_to_holonomic(sin(x), x) + expr_to_holonomic(cos(x), x),
    'float_cos': lambda: expr_to_holonomic(1.2*cos(3.1*x), x),
    'float_param_sin': lambda: expr_to_holonomic(1.5*a*sin(x), x),
}
label, mode, case = sys.argv[1], sys.argv[2], sys.argv[3]
f = CASES[case]
def run():
    t = time.perf_counter()
    try:
        f(); ok = 'ok'
    except Exception as e:
        ok = type(e).__name__
    return time.perf_counter() - t, ok
if mode == 'first':
    t, ok = run()
else:
    run(); ts = []
    for _ in range(7):
        t, ok = run(); ts.append(t)
    t = statistics.median(ts)
print(json.dumps({'label': label, 'mode': mode, 'case': case, 'seconds': t, 'status': ok}))
