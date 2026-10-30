"""Check every dsolve mismatch of a symox ODE corpus run against ref.jsonl with SymPy's checkodesol (both solutions)."""
import json, os, sys, signal
SS = '/home/keckj/Documents/symox-worktrees/.tooling/sspeed'
import sympy as M
from sympy import checkodesol
ns = {n: getattr(M, n) for n in dir(M) if not n.startswith('_')}
ref = {json.loads(l)['name']: json.loads(l) for l in open(SS + '/ref.jsonl')}
port = [json.loads(l) for l in open(sys.argv[1])]
cases = {c['name']: c for c in json.load(open(SS + '/cases.json'))}
out = open(sys.argv[2], 'w')
def check(eqs, sol, func):
    if sol is None or sol.startswith('ERR') : return 'n/a'
    try: s = M.sympify(sol, locals=ns)
    except Exception as e: return 'parse:%s' % type(e).__name__
    sols = list(s) if isinstance(s, (list, tuple)) else [s]
    res = []
    for q in sols:
        try: r = checkodesol(eqs, q, func, solve_for_func=False)
        except Exception as e: r = ('exc', type(e).__name__)
        res.append(r[0] if isinstance(r, tuple) and r[0] in (True, False) else str(r)[:40])
    return res
def run(c, key):
    eq = eval(c['eq'], ns); func = eval(c['func'], ns)
    rs = ref[c['name']]
    return {'ref': check(eq, rs.get(key), func), 'port': check(eq, P[c['name']].get(key), func)}
P = {r['name']: r for r in port}
for r in port:
    s = ref[r['name']]
    for key in ('hinted', 'default'):
        if r.get(key) == s.get(key) or 'fatal' in s: continue
        rd, wr = os.pipe(); pid = os.fork()
        if pid == 0:
            os.close(rd); signal.alarm(int(os.environ.get('TO', '90')))
            try: res = run(cases[r['name']], key)
            except BaseException as e: res = {'fatal': repr(e)[:100]}
            os.write(wr, json.dumps(res, default=str).encode()); os._exit(0)
        os.close(wr); data = b''
        while True:
            ch = os.read(rd, 65536)
            if not ch: break
            data += ch
        os.close(rd); os.waitpid(pid, 0)
        try: res = json.loads(data)
        except Exception: res = {'fatal': 'timeout'}
        res.update(name=r['name'], key=key, refsol=str(s.get(key))[:200], portsol=str(r.get(key))[:200])
        out.write(json.dumps(res) + '\n'); out.flush()
