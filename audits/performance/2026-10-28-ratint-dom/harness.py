"""ratint timings: first call in a fresh process and median of repeated calls.

usage: python harness.py <label> <case-index>   (one process per case and label)
"""
import json, statistics, sys, time
import symox as sx
from symox.integrals.rationaltools import ratint

x = sx.Symbol('x'); a = sx.Symbol('a')
CASES = ['1/(x**3+x+1)', '1/(x**3-x+1)', '1/((x**2+1)*(x**3+2))', '1/(x**5+1)',
         '1/(x**4+2*x**2+2)', '1/(x**2+1)', '(x+1)/(x**2+x+1)**2', '1/(x**2+a)',
         '1/(x**2+I)**2', 'I/(x**2+1)**2', '(x+I)/(x**3-I)']
label, i = sys.argv[1], int(sys.argv[2])
f = sx.sympify(CASES[i], locals={'x': x, 'a': a})
t0 = time.perf_counter()
try:
    r = ratint(f, x); out = len(str(r))
except Exception as e:
    out = type(e).__name__
first = time.perf_counter() - t0
reps = []
if isinstance(out, int):
    for _ in range(5):
        t0 = time.perf_counter(); ratint(f, x); reps.append(time.perf_counter() - t0)
print(json.dumps({'label': label, 'case': CASES[i], 'first_s': first,
                  'median_s': statistics.median(reps) if reps else None, 'str_len': out}))
