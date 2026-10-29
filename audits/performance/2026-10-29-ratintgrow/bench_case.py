"""One fresh-process measurement: python bench_case.py <case> -> JSON line."""
import json, statistics, sys, time

case = sys.argv[1]
from symox import symbols, sympify, srepr
from symox.integrals.rationaltools import ratint
import re

x = symbols('x')
CASES = {
    'x6p8': '1/(x**6+8)', 'x4m2x2p4': '1/(x**4-2*x**2+4)', 'x4p2': '1/(x**4+2)',
    'x8p1': '1/(x**8+1)', 'x8p16': '1/(x**8+16)', 'x4pxp1': '1/(x**4+x+1)',
    'x2_x6p8': 'x**2/(x**6+8)', 'x4p1': '1/(x**4+1)', 'x4mx2p1': '1/(x**4-x**2+1)',
    'xm1_4p1': '1/((x-1)**4+1)', 'x6p27': '1/(x**6+27)', 'x6p1': '1/(x**6+1)',
    'x3pxp1': '1/(x**3+x+1)', 'x3m3xp1': '1/(x**3-3*x+1)', 'x5p1': '1/(x**5+1)',
}
expr = sympify(CASES[case], locals={'x': x})
f = lambda: ratint(expr, x)
t0 = time.perf_counter(); result = f(); first = time.perf_counter() - t0
n = 1
while True:
    t0 = time.perf_counter()
    for _ in range(n): f()
    if time.perf_counter() - t0 > 0.05 or n > 100000: break
    n *= 2
reps = []
for _ in range(7):
    t0 = time.perf_counter()
    for _ in range(n): f()
    reps.append((time.perf_counter() - t0) / n)
text = re.sub(r"dummy_index=\d+", "", srepr(result))
print(json.dumps({'case': case, 'first_s': first, 'repeat_s': statistics.median(reps), 'result': text[:4000]}))
