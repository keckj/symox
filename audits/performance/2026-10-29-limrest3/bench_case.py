"""bench_case.py <case>: one fresh process, one public call.

Prints one JSON line: first-call time (ms), median of 7 batches of repeated
calls (us per call) and the printed result.
"""
import json
import statistics
import sys
import time

from symox import *  # noqa: F401,F403
from symox.series.gruntz import gruntz  # noqa: F401

x, n, a = symbols('x n a')
m = Symbol('m', positive=True)
y = Symbol('y', positive=True)
z = Symbol('z', negative=True)
ar = Symbol('a', real=True)

CASES = {
    # held (Limit/Integral) on the base, computed now
    'ei': lambda: limit((Ei(x - exp(-exp(x))) - Ei(x)) * exp(-x) * exp(exp(x)) * x, x, oo),
    'erf': lambda: limit((erf(x - exp(-exp(x))) - erf(x)) * exp(exp(x)) * exp(x**2), x, oo),
    'ugam': lambda: limit(uppergamma(x, 1) / gamma(x), x, oo),
    'subfac': lambda: limit(subfactorial(x), x, oo),
    'stirl': lambda: limit(((x + S(1)/2) * log(x) - x + log(2*pi)/2 - log(factorial(x)) + S(1)/(12*x)) * x**3, x, oo),
    'harm_m1': lambda: limit(harmonic(n, m + 1), n, oo),
    'harm_z1': lambda: limit(harmonic(x, z + 1), x, oo),
    'atan_lim': lambda: limit(ar*(log(x) - log(ar**2*x**2 + 1)/2) - atan(ar*x)/x, x, oo),
    'atan_int': lambda: integrate(atan(ar*x)/x**2, (x, 1, oo)),
    # already computed on the base (controls)
    'c_sinx': lambda: limit(sin(x)/x, x, 0),
    'c_e': lambda: limit((1 + 1/x)**x, x, oo),
    'c_gruntz': lambda: limit(exp(x)*(sin(1/x + exp(-x)) - sin(1/x)), x, oo),
    'c_ei': lambda: limit(Ei(x)*exp(-x)*x, x, oo),
    'c_erf': lambda: limit((1 - erf(x))*exp(x**2)*x, x, oo),
    'c_harm3': lambda: limit(harmonic(x, 3), x, oo),
    'c_lgam': lambda: limit(loggamma(x)/(x*log(x)), x, oo),
    'c_atan': lambda: limit(atan(x)*x - pi*x/2, x, oo),
}

if __name__ == '__main__':
    k = sys.argv[1]
    f = CASES[k]
    t = time.perf_counter()
    r = f()
    first = (time.perf_counter() - t) * 1e3
    reps = max(1, min(50, int(0.2 / max(first / 1e3, 1e-6))))
    batches = []
    for _ in range(7):
        t = time.perf_counter()
        for _ in range(reps):
            f()
        batches.append((time.perf_counter() - t) / reps * 1e6)
    print(json.dumps({'case': k, 'first_ms': first, 'repeat_us': statistics.median(batches),
                      'reps': reps, 'result': str(r)}))
