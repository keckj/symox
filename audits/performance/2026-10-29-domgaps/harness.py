"""Time public calls touched by the domain element comparison change.

Usage: PYTHONPATH=<tree>/python python harness.py <label>
Prints one JSON line per case: first call (fresh process) and the median of
repeated calls, in microseconds. Errors are recorded by type name, since
several cases raise before the change.
"""
import json
import statistics
import sys
import time

t0 = time.perf_counter()
from symox import QQ, ZZ, symbols, sympify, sqrt, Matrix
from symox.polys.matrices import DomainMatrix
from symox.polys.matrices.normalforms import smith_normal_form, invariant_factors
x, y = symbols('x y')
FX, PX, FXY = QQ.frac_field(x), QQ[x], QQ.frac_field(x, y)


def dm(rows, K):
    return DomainMatrix([[K.from_sympy(sympify(v)) for v in r] for r in rows], (len(rows), len(rows[0])), K)


Z3 = [['x', 0, 1], [0, 0, 0], [1, 0, 'x']]
FULL3 = [['x', 1, 'x + 1'], [1, 'x', 2], ['x**2', 0, 1]]
e, f, g = FX.from_sympy(x / (x + 1)), FX.from_sympy(sympify(2)), PX.from_sympy(x**2 + 1)
CASES = {
    'charpoly QQ(x) 2x2 zero column': lambda: dm([['x', 0], [1, 0]], FX).charpoly(),
    'charpoly QQ(x) 3x3 zero column': lambda: dm(Z3, FX).charpoly(),
    'charpoly QQ(x) full 3x3': lambda: dm(FULL3, FX).charpoly(),
    'charpoly QQ[x] zero column': lambda: dm(Z3, PX).charpoly(),
    'smith QQ(x) rank deficient': lambda: smith_normal_form(dm(Z3, FX)),
    'invariant_factors QQ(x) full': lambda: invariant_factors(dm(FULL3, FX)),
    'rref QQ(x) 3x3': lambda: dm(FULL3, FX).rref(),
    'det QQ(x,y) 2x2': lambda: dm([['x', 'y'], ['y', 'x']], FXY).det(),
    'Matrix.charpoly QQ(x) entries': lambda: Matrix([[x, 0, 1], [0, 0, 0], [1, 0, x]]).charpoly(),
    'QQ(x) elem == 0 x1000': lambda: [e == 0 for _ in range(1000)],
    'QQ(x) elem == 2 x1000': lambda: [f == 2 for _ in range(1000)],
    'QQ[x] elem != 1 x1000': lambda: [g != 1 for _ in range(1000)],
    'QQ(x) elem < elem x1000': lambda: [e < f for _ in range(1000)],
}


def run(fn):
    try:
        fn()
        return None
    except Exception as exc:
        return type(exc).__name__


label = sys.argv[1]
print(json.dumps({'label': label, 'case': 'import', 'first_us': (time.perf_counter() - t0) * 1e6}))
for name, fn in CASES.items():
    t = time.perf_counter(); err = run(fn); first = (time.perf_counter() - t) * 1e6
    times = []
    for _ in range(15):
        t = time.perf_counter(); run(fn); times.append((time.perf_counter() - t) * 1e6)
    print(json.dumps({'label': label, 'case': name, 'first_us': round(first, 1),
                      'median_us': round(statistics.median(times), 1), 'error': err}))
