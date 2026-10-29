"""Time one call: first call in this fresh process, then the median of repeats.

usage: python bench.py sympy|symox REPEATS EXPR
Prints: module, expr, first-call seconds, median repeat seconds, digest of
str(result) (sha256, first 16 hex digits) and its length.
The big operands are module-level Python ints (A, B: about 10**5 digits with
a common factor G) so that their construction is not timed.
"""
import functools
import hashlib
import operator
import statistics
import sys
import time

sys.set_int_max_str_digits(0)
mod, repeats, case = sys.argv[1], int(sys.argv[2]), sys.argv[3]
if mod == 'sympy':
    import sympy as S
else:
    import symox as S
ns = {name: getattr(S, name) for name in dir(S) if not name.startswith('_')}
G = 13**30000
ns.update(functools=functools, operator=operator, G=G,
          A=G * 17**40000 + G, B=G * 19**38000 * 3 + 5 * G, x=S.Symbol('x'))
t = time.perf_counter()
value = eval(case, ns)
first = time.perf_counter() - t
times = []
for _ in range(repeats):
    t = time.perf_counter()
    eval(case, ns)
    times.append(time.perf_counter() - t)
rep = statistics.median(times) if times else float('nan')
text = str(value)
tag = 'len=%d sha=%s' % (len(text), hashlib.sha256(text.encode()).hexdigest()[:16])
print('%s\t%s\t%.4f\t%.4f\t%s' % (mod, case, first, rep, tag), flush=True)
