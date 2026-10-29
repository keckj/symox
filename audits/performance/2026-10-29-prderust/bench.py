"""prde benchmark: risch_integrate on 10 integrands that reach prde.

usage: python bench.py {symox|sympy} OUT.json [REPEATS]
Each integrand runs in a fresh subprocess: the first call is timed once,
then REPEATS further calls give the repeated-call median.
"""
import json, os, statistics, subprocess, sys, time

CASES = ['log(x)**2', 'x*log(x)', 'log(x)**3/x', 'x**2*log(x)**2', 'x*log(x+1)',
         'log(x**2+1)', 'x*log(x)**3', 'log(x)**4', '(x+1)*log(x)**2', 'log(x)**2/x**2']

CHILD = r'''
import sys, time, json
lib = sys.argv[1]; src = sys.argv[2]; n = int(sys.argv[3])
mod = __import__(lib)
from importlib import import_module
risch = import_module(lib + '.integrals.risch')
x = mod.Symbol('x'); e = mod.sympify(src, locals={'x': x})
t = time.perf_counter(); r = risch.risch_integrate(e, x); first = time.perf_counter() - t
rep = []
for _ in range(n):
    t = time.perf_counter(); risch.risch_integrate(e, x); rep.append(time.perf_counter() - t)
print(json.dumps({'first': first, 'rep': rep, 'result': str(r)}))
'''

def main():
    lib, out = sys.argv[1], sys.argv[2]
    n = int(sys.argv[3]) if len(sys.argv) > 3 else 7
    runs = int(os.environ.get('BENCH_RUNS', '3'))
    data = {}
    for c in CASES:
        firsts, reps, res = [], [], None
        for _ in range(runs):
            p = subprocess.run([sys.executable, '-c', CHILD, lib, c, str(n)],
                               capture_output=True, text=True, timeout=600, cwd='/tmp')
            d = json.loads(p.stdout.strip().splitlines()[-1])
            firsts.append(d['first']); reps.extend(d['rep']); res = d['result']
        data[c] = {'first_median': statistics.median(firsts), 'repeat_median': statistics.median(reps),
                   'firsts': firsts, 'repeats': reps, 'result': res}
        print('%-18s first %8.2f ms  repeat %8.2f ms' % (c, 1e3 * data[c]['first_median'], 1e3 * data[c]['repeat_median']), flush=True)
    json.dump(data, open(out, 'w'), indent=1)

if __name__ == '__main__':
    main()
