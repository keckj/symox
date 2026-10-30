"""Differential: native checksol vs logged SymPy 1.14 checksol results.
usage: diff_checksol.py out.json [legacy] < pairs.jsonl"""
import json
import sys
import time
import collections

import symox
from symox import *  # noqa
from symox import _core

ns = dict(vars(symox))
legacy = 'legacy' in sys.argv[2:]
recs = {}
for line in sys.stdin:
    r = json.loads(line)
    recs[(r['f'], r['x'], r['v'], tuple(sorted(r['flags'].items())))] = r
out = []
prog = open(sys.argv[1] + '.progress', 'w')
for i, r in enumerate(recs.values()):
    prog.write('%d %s\n' % (i, r['test']))
    prog.flush()
    try:
        f, x, v = eval(r['f'], ns), eval(r['x'], ns), eval(r['v'], ns)
    except Exception as e:
        out.append(dict(r, n='parse:' + type(e).__name__))
        continue
    fl = r['flags']
    t0 = time.perf_counter()
    try:
        n = _core._checksol_native(f, x, v, fl.get('numerical', True), fl.get('force', True),
                                   fl.get('simplify', True), fl.get('minimal', False), *(() if not legacy else ()))
    except Exception as e:
        n = 'err:' + type(e).__name__
    out.append(dict(r, n=n, t=time.perf_counter() - t0))
json.dump(out, open(sys.argv[1], 'w'))
c = collections.Counter((str(o['r']), str(o['n'])) for o in out)
print(sorted(c.items()))
bad = [o for o in out if o['n'] is False and o['r'] is not False]
print('native False, SymPy not False:', len(bad))
