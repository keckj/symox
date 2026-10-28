"""Run one polycap1 case: first call in this fresh process, then repeats.

usage: python bench.py <case> [repeats]   (prints one JSON line)
The library is symox from PYTHONPATH, or the reference sympy when
BENCH_LIB=sympy (run from /tmp). Run under timeout and a memory cap.
"""
import json, sys, time

import os
if os.environ.get('BENCH_LIB') == 'sympy':
    import sympy as lib
else:
    import symox as lib
from_lib = lib.__name__
x = lib.Symbol('x')
R, I = lib.Rational, lib.I

CASES = {
    'count_5000': lambda: lib.count_roots(x**5000 - x - 1),
    'count_5000_ab': lambda: lib.count_roots(x**5000 - x - 1, 0, 2),
    'count_5000_rat': lambda: lib.count_roots(x**5000 - x - 1, R(-3, 2), R(7, 5)),
    'count_5000_float': lambda: lib.count_roots(x**5000 - x - 1, -1.5, lib.Float('1.25')),
    'count_3000_sparse': lambda: lib.count_roots(x**3000 - 3*x**7 + 2*x - 5),
    'count_cheb120': lambda: lib.count_roots(lib.chebyshevt(120, x)),
    'count_cheb120_ab': lambda: lib.count_roots(lib.chebyshevt(120, x), R(1, 3), 1),
    'count_poly_sq': lambda: lib.count_roots(lib.expand((x**40 - 2)**2*(x**3 - x + R(1, 7)))),
    'count_complex_60': lambda: lib.count_roots(x**60 - x - 1, -2 - 2*I, 2 + 2*I),
    'factor_300': lambda: lib.factor_list(x**300 + x + 1),
    'factor_200': lambda: lib.factor_list(x**200 + x + 1),
    'factor_250': lambda: lib.factor_list(x**250 + x + 1),
    'factor_301': lambda: lib.factor_list(x**301 + x + 1),
    'factor_expr_300': lambda: lib.factor(x**300 + x + 1),
    'factor_sparse_256': lambda: lib.factor_list(x**256 - x**5 + 3),
    'factor_sparse_210': lambda: lib.factor_list(x**210 + 2*x**3 - 7),
}


def run(case, repeats):
    f = CASES[case]
    times, out = [], None
    for _ in range(1 + repeats):
        t = time.perf_counter()
        try:
            out = f()
        except Exception as e:  # recorded, never counted as a speed-up
            out = 'ERROR %s: %s' % (type(e).__name__, e)
        times.append(time.perf_counter() - t)
    first, rest = times[0], sorted(times[1:])
    med = rest[len(rest) // 2] if rest else None
    print(json.dumps({'case': case, 'lib': from_lib, 'first': first,
                      'repeat_median': med, 'result': str(out)[:400]}))


if __name__ == '__main__':
    run(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 2)
