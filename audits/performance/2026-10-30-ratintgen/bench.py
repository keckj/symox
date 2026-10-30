"""Public ratint/integrate timings for the ratintgen families.

Usage: bench.py sympy|symox OUT.json
Run SymPy from /tmp; run symox from the worktree root with PYTHONPATH=python.
Each case: first call in a fresh process (subprocess), then the median of 7
repeated calls in one process (after one warm-up call).
"""
import json, subprocess, sys, statistics, time

CASES = [
    ('integrate', '1/(x**2+I)**2'), ('integrate', 'I/(x**2+1)**2'),
    ('integrate', '(x+I)/(x**3-I)'), ('integrate', '1/(x**2+2*I*x-5)'),
    ('integrate', '1/(x**2+a)**2'), ('integrate', '1/((x-a)*(x-b))'),
    ('integrate', 'x/(x**2+a*x+b)**2'), ('integrate', 'a/(x-a)**3'),
    ('integrate', '1/(x-sqrt(2))**2/(x+1)'), ('integrate', '1/(x**2+sqrt(2))**2'),
    ('ratint', '1/(x**2+a)**3'), ('ratint', '(x**2+1)/(x**2+I)**3'),
]

def one(engine, kind, expr, repeat):
    lib = 'sympy' if engine == 'sympy' else 'symox'
    m = __import__(lib)
    from importlib import import_module
    ratint = import_module(lib + '.integrals.rationaltools').ratint
    x, a, b = m.symbols('x a b')
    e = m.sympify(expr, locals={'x': x, 'a': a, 'b': b})
    fn = (lambda: m.integrate(e, x)) if kind == 'integrate' else (lambda: ratint(e, x))
    t0 = time.perf_counter(); r = fn(); first = time.perf_counter() - t0
    times = []
    for _ in range(repeat):
        t0 = time.perf_counter(); fn(); times.append(time.perf_counter() - t0)
    return first, (statistics.median(times) if times else None), str(r)

if __name__ == '__main__':
    if len(sys.argv) == 5:
        first, rep, r = one(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]))
        print(json.dumps([first, rep, r])); sys.exit()
    engine, out = sys.argv[1], sys.argv[2]
    rows = []
    for kind, expr in CASES:
        fresh = json.loads(subprocess.check_output([sys.executable, __file__, engine, kind, expr, '0']))
        warm = json.loads(subprocess.check_output([sys.executable, __file__, engine, kind, expr, '7']))
        rows.append({'call': '%s(%s, x)' % (kind, expr), 'fresh_s': fresh[0], 'repeat_median_s': warm[1], 'result': warm[2]})
        print('%-40s fresh %8.1f ms  repeat %8.2f ms' % (rows[-1]['call'], fresh[0]*1e3, warm[1]*1e3), flush=True)
    json.dump(rows, open(out, 'w'), indent=1)
