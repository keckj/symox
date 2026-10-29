"""Before/after timings of resultant(..., includePRS=True) over ZZ.

usage: PYTHONPATH=<tree>/python python bench.py <case> <repeats>
Prints: case first_call_seconds median_repeat_seconds sha256(str(result))[:16]
"""
import hashlib
import statistics
import sys
import time

sys.set_int_max_str_digits(0)
from symox import Poly, expand, resultant, subresultants, symbols, factor_list

x = symbols('x')


def big(n, m):
    return expand((x + 1)**n), expand((x - 2)**m + 1)


CASES = {
    'prs_200_150': lambda: resultant(*big(200, 150), x, includePRS=True),
    'prs_150_200': lambda: resultant(*reversed(big(200, 150)), x, includePRS=True),
    'prs_120_90': lambda: resultant(*big(120, 90), x, includePRS=True),
    'prs_60_45': lambda: resultant(*big(60, 45), x, includePRS=True),
    'prs_40_30': lambda: resultant(*big(40, 30), x, includePRS=True),
    'subres_60_45': lambda: subresultants(*big(60, 45), x),
    'poly_prs_60_45': lambda: Poly(big(60, 45)[0], x).resultant(Poly(big(60, 45)[1], x), includePRS=True),
    'prs_small': lambda: resultant(x**10 - 2*x**3 + 1, x**7 + 5*x - 11, x, includePRS=True),
    'prs_abnormal': lambda: resultant(expand((x - 1)**3*(x + 2)*(x**2 + 3)), expand((x - 1)**2*(x**2 + 5)), x, includePRS=True),
    'prs_dense_30': lambda: resultant(sum((k + 1)*x**k for k in range(31)), sum((-1)**k*(2*k + 3)*x**k for k in range(25)), x, includePRS=True),
    'factor_2000': lambda: factor_list(x**2000 + x + 1),
}


def run(case, repeats):
    f = CASES[case]
    t = time.perf_counter()
    try:
        r = f()
        digest = hashlib.sha256(str(r).encode()).hexdigest()[:16]
    except Exception as e:
        digest = '%s: %s' % (type(e).__name__, str(e)[:80])
    first = time.perf_counter() - t
    times = []
    for _ in range(repeats):
        t = time.perf_counter()
        try:
            f()
        except Exception:
            pass
        times.append(time.perf_counter() - t)
    med = statistics.median(times) if times else float('nan')
    print('%s %.4f %.4f %s' % (case, first, med, digest), flush=True)


if __name__ == '__main__':
    run(sys.argv[1], int(sys.argv[2]))
