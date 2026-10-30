import json, sys, time, signal, os
import sympy as M
ns = {n: getattr(M, n) for n in dir(M) if not n.startswith('_')}
cases = {c['name']: c for c in json.load(open('/home/keckj/Documents/symox-worktrees/.tooling/sspeed/cases.json'))}
names = sys.argv[1].split(',')
for n in names:
    c = cases[n]
    pid = os.fork()
    if pid == 0:
        signal.alarm(120)
        eq = eval(c['eq'], ns); func = eval(c['func'], ns)
        t = time.time()
        try:
            r = [str(h) for h in M.classify_ode(eq, func)]
        except Exception as e:
            r = 'ERR %s' % type(e).__name__
        print(json.dumps({'name': n, 'classify': r, 't': round(time.time()-t, 2)}), flush=True)
        os._exit(0)
    _, st = os.waitpid(pid, 0)
    if st: print(json.dumps({'name': n, 'classify': 'timeout/killed', 'status': st}), flush=True)
