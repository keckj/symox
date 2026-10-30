"""usage: bench.py sympy|symox OUT.json [runs]

classify_ode(eq, func) and dsolve(eq, func) on 12 ODE corpus cases.
For each (case, operation) and run: a fresh process measures the first call
(the module is imported before timing), then the median of 5 repeated calls
in the same process. The reported value is the median over the runs.
BENCH_CASES (comma-separated case names) and BENCH_OPS (classify,dsolve)
override the defaults.
"""
import json, os, subprocess, sys, statistics

CASES = ['separable_01', 'separable_03', '1st_exact_06', 'linear_01', 'bernoulli_01',
         'riccati_01', '1st_homogeneous_coeff_best_01', 'almost_lin_02', 'linear_coeff_01',
         'separable_reduced_02', 'rational_riccati_02', 'fact_02']
if os.environ.get('BENCH_CASES'):
    CASES = os.environ['BENCH_CASES'].split(',')
OPS = tuple(os.environ.get('BENCH_OPS', 'classify,dsolve').split(','))
CORPUS = '/home/keckj/Documents/symox-worktrees/.tooling/sspeed/cases.json'

CHILD = r'''
import json, sys, time, statistics
mod, name, op = sys.argv[1:4]
M = __import__(mod)
ns = {n: getattr(M, n) for n in dir(M) if not n.startswith('_')}
c = {c['name']: c for c in json.load(open(%r))}[name]
eq = eval(c['eq'], ns); func = eval(c['func'], ns)
call = (lambda: M.classify_ode(eq, func)) if op == 'classify' else (lambda: M.dsolve(eq, func))
t = time.perf_counter(); call(); first = time.perf_counter() - t
rep = []
for _ in range(5):
    t = time.perf_counter(); call(); rep.append(time.perf_counter() - t)
print(json.dumps({'first': first, 'repeat': statistics.median(rep)}))
''' % CORPUS


def main():
    mod, out = sys.argv[1], sys.argv[2]
    runs = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    py = sys.executable
    res = {}
    for name in CASES:
        for op in OPS:
            samples = []
            for _ in range(runs):
                try:
                    p = subprocess.run([py, '-c', CHILD, mod, name, op], capture_output=True, text=True, timeout=300)
                    samples.append(json.loads(p.stdout.strip().splitlines()[-1]))
                except Exception as exc:
                    samples.append({'error': type(exc).__name__})
            ok = [s for s in samples if 'first' in s]
            res['%s/%s' % (name, op)] = {
                'samples': samples,
                'first_median': statistics.median([s['first'] for s in ok]) if ok else None,
                'repeat_median': statistics.median([s['repeat'] for s in ok]) if ok else None,
            }
            print(name, op, res['%s/%s' % (name, op)]['first_median'], res['%s/%s' % (name, op)]['repeat_median'], flush=True)
    json.dump(res, open(out, 'w'), indent=1)


main()
