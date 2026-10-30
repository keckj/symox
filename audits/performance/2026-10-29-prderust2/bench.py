"""Time one risch_integrate call in a fresh process.

usage: python bench.py <symox|sympy> <index>   (prints: index seconds result-srepr-hash)
Run by run.sh, which starts a fresh interpreter for every (config, index, repetition).
"""
import sys, time, importlib, hashlib
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
dt = time.perf_counter() - t0
print(i, '%.6f' % dt, hashlib.sha1(S.srepr(r).encode()).hexdigest()[:10], INTEGRANDS[i])
