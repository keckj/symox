"""Time one call: first call in this fresh process, then the median of repeats.

usage: python bench.py sympy|symox REPEATS EXPR
Prints: module, expr, first-call seconds, median repeat seconds, digest.
Repeats use a new argument each time is not possible for cached results, so
the repeats call the same expression; SymPy caches some of these values, so
its repeat times may reflect its cache (reported as measured).
"""
import hashlib
import statistics
import sys
import time

sys.set_int_max_str_digits(0)
mod, repeats, case = sys.argv[1], int(sys.argv[2]), sys.argv[3]
t0 = time.perf_counter()
if mod == 'sympy':
    import sympy as S
else:
    import symox as S
ns = {name: getattr(S, name) for name in dir(S) if not name.startswith('_')}
t = time.perf_counter()
value = eval(case, ns)
first = time.perf_counter() - t
times = []
for _ in range(repeats):
    t = time.perf_counter()
    eval(case, ns)
    times.append(time.perf_counter() - t)
rep = statistics.median(times) if times else float('nan')
if getattr(value, 'is_Rational', False) or isinstance(value, int):
    text = str(value)
    tag = 'len=%d sha=%s' % (len(text), hashlib.sha256(text.encode()).hexdigest()[:16])
else:
    tag = 'UNEVAL'
print('%s\t%s\t%.4f\t%.4f\t%s' % (mod, case, first, rep, tag), flush=True)
