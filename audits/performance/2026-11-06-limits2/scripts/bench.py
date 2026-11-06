"""bench.py OUT.json: time the public oscillating limits and the bounds
arithmetic of the symox found on PYTHONPATH (best of REPEAT runs, a fresh
symbol each run so no cached result is reused)."""
import json
import sys
import time

from symox import (E, Rational, Symbol, cos, exp, frac, limit, log, oo, pi,
                   sin, tan)
from symox.calculus.accumulationbounds import AccumBounds as AB

REPEAT = 5

LIMITS = {
    'sin(x), oo': lambda x: limit(sin(x), x, oo),
    'sin(1/x), 0': lambda x: limit(sin(1/x), x, 0),
    'cos(1/x), 0, +-': lambda x: limit(cos(1/x), x, 0, '+-'),
    'exp(sin(x)), oo': lambda x: limit(exp(sin(x)), x, oo),
    'sin(x)**2, oo': lambda x: limit(sin(x)**2, x, oo),
    '1/(2 + sin(x)), oo': lambda x: limit(1/(2 + sin(x)), x, oo),
    'sin(x) + cos(x), oo': lambda x: limit(sin(x) + cos(x), x, oo),
    'frac(x), oo': lambda x: limit(frac(x), x, oo),
    'tan(x), oo': lambda x: limit(tan(x), x, oo),
    'log(2 + sin(x)), oo': lambda x: limit(log(2 + sin(x)), x, oo),
    'x*sin(x), oo': lambda x: limit(x*sin(x), x, oo),
    'sin(x)**x, oo (held)': lambda x: limit(sin(x)**x, x, oo),
}


def arithmetic(x):
    a, b = AB(1, 2), AB(-1, 1)
    for k in range(1, 30):
        (a + k, a*-k, k/a, b**2, a**AB(1, 2), 2**a, a/b, a - x, AB(0, 1)**oo,
         a < Rational(3, 2), a.intersection(AB(0, k)), sin(AB(0, pi)),
         exp(AB(0, k)), log(AB(1, k + 1)))


CASES = dict(LIMITS, **{'bounds arithmetic (29 x 14 ops)': arithmetic})


def best(f):
    t = []
    for i in range(REPEAT):
        y = Symbol('x%d' % i)
        start = time.perf_counter()
        f(y)
        t.append(time.perf_counter() - start)
    return min(t)


out = {name: round(best(f)*1e3, 3) for name, f in CASES.items()}
json.dump(out, open(sys.argv[1], 'w'), indent=1)
for name, ms in out.items():
    print('%-36s %9.3f ms' % (name, ms))
