"""Time risch_integrate in a fresh process: first call, then the median of 3 repeated calls.

usage: python bench.py <symox|sympy> <index>
prints: index first-seconds repeated-seconds sha1-prefix-of-srepr(result) integrand sha1-prefix-of-str(result)
"""
import sys, time, importlib, hashlib, statistics
mod = sys.argv[1]
i = int(sys.argv[2])
S = importlib.import_module(mod)
risch = importlib.import_module(mod + '.integrals.risch')
x = S.Symbol('x')
INTEGRANDS = ['x*log(x)', 'log(x)**2', 'log(x)**3/x', 'log(log(x))', 'exp(x)*x',
              'x*exp(x)*log(x)', '1/(x*log(x))', 'log(x)**2/x**2',
              'exp(x)/(1+exp(x))**2', 'log(x)*log(x+1)']
f = S.sympify(INTEGRANDS[i], locals={'x': x})
t0 = time.perf_counter()
r = risch.risch_integrate(f, x)
first = time.perf_counter() - t0
reps = []
for _ in range(3):
    t0 = time.perf_counter()
    risch.risch_integrate(f, x)
    reps.append(time.perf_counter() - t0)
print(i, '%.6f' % first, '%.6f' % statistics.median(reps),
      hashlib.sha1(S.srepr(r).encode()).hexdigest()[:10], INTEGRANDS[i],
      hashlib.sha1(str(r).encode()).hexdigest()[:10])
