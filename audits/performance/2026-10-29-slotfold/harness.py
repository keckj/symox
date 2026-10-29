"""piecewise_fold / Sum / Product / Integral constructor timings.

usage: python harness.py sympy|symox out.jsonl   (sympy: run from /tmp)
Each case: first call (fresh process, after warm-up of unrelated objects)
and the median of 7 repeated batches. Results printed with str().
"""
import json, sys, time
mod = sys.argv[1]
t0 = time.perf_counter()
if mod == 'sympy':
    import sympy as S
else:
    import symox as S
import_s = time.perf_counter() - t0
from_names = ('Piecewise Sum Product Integral piecewise_fold sin symbols').split()
Piecewise, Sum, Product, Integral, piecewise_fold, sin, symbols = (getattr(S, n) for n in from_names)
x, y, n, k = symbols('x y n k')
CASES = [
    ('fold_mul', lambda: piecewise_fold(Piecewise((x, x > 0), (0, True))*n)),
    ('sum_free', lambda: Sum(Piecewise((x, x > 0), (0, True))*n, (n, 1, 3))),
    ('sum_bound', lambda: Sum(Piecewise((n, n > 2), (0, True))*x, (n, 1, k))),
    ('integral_add', lambda: Integral(Piecewise((x, x < 1), (x**2, True))*y + 1, (x, 0, 2))),
    ('product_two', lambda: Product(Piecewise((n, n > 1), (1, True))*Piecewise((2, n < 3), (n, True)), (n, 1, 4))),
    ('fold_common', lambda: piecewise_fold(Piecewise((x, x < 0), (1, True)) + Piecewise((x**2, x < 0), (2, True)))),
    ('fold_nested', lambda: piecewise_fold(Piecewise((Piecewise((x, y > 0), (1, True)), x < 1), (0, True)))),
    ('fold_func', lambda: piecewise_fold(sin(Piecewise((x, x > 0), (0, True))))),
    ('fold_noeval', lambda: piecewise_fold(Piecewise((x, x > 0), (0, True))*Piecewise((1, x > 0), (2, True)), evaluate=False)),
    ('integral_two', lambda: Integral(Piecewise((1, x > y), (0, True))*Piecewise((x, y > 0), (0, True)), (x, 0, 1))),
]
out = []
for name, fn in CASES:
    t = time.perf_counter(); r = fn(); first = time.perf_counter() - t
    N = 200
    batches = []
    for _ in range(7):
        t = time.perf_counter()
        for _ in range(N): fn()
        batches.append((time.perf_counter() - t) / N)
    batches.sort()
    out.append({'case': name, 'first_us': first*1e6, 'median_us': batches[3]*1e6, 'result': str(r)})
with open(sys.argv[2], 'w') as f:
    for row in out: f.write(json.dumps(row) + '\n')
for row in out: print('%-14s first %9.1f us  median %8.1f us  %s' % (row['case'], row['first_us'], row['median_us'], row['result']))
