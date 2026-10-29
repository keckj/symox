"""One fresh-process measurement: python bench_case.py <case> -> JSON line."""
import json, statistics, sys, time

case = sys.argv[1]
from symox import symbols, sympify, srepr
from symox.integrals.rationaltools import ratint
import re

x = symbols('x')
CASES = {
    'lrt_x4p1': '1/(x**4+1)', 'lrt_x4p4': '1/(x**4+4)', 'lrt_x4pq': '1/(x**4+1/4)',
    'lrt_x4mx2p1': '1/(x**4-x**2+1)', 'lrt_mixed': '1/((x**2+1)*(x**2+4)*(x-1))',
    'lrt_x4_x6p1': 'x**4/(x**6+1)', 'lrt_x6p1': '1/(x**6+1)',
    'lrt_x3p1_x_1': '1/(x**3+x+1)', 'lrt_x3m2': '1/(x**3-2)', 'lrt_x5p1': '1/(x**5+1)',
    'lrt_x6m2': '1/(x**6-2)', 'lrt_x8p1': '1/(x**8+1)', 'lrt_x2_x4px1': 'x**2/(x**4+x+1)',
    'lrt_x5p2': '1/(x**5+2)', 'lrt_x7p1': '1/(x**7+1)',
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
