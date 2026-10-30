"""Replay the logged SymPy 1.14 checksol pairs of the solvegaps wave through
the internal native check (_core._checksol_native) and the public checksol.
usage: replay_pairs.py out.json < pairs.jsonl"""
import collections
import json
import signal
import sys
import time

import symox
from symox import _core
from symox.solvers.solvers import checksol

ns = dict(vars(symox))
from symox.functions.elementary.piecewise import ExprCondPair  # noqa: E402
ns['ExprCondPair'] = ExprCondPair


class Timeout(Exception):
    pass


def _alarm(*a):
    raise Timeout()


signal.signal(signal.SIGALRM, _alarm)
recs = {}
for line in sys.stdin:
    r = json.loads(line)
    recs[(r['f'], r['x'], r['v'], tuple(sorted(r['flags'].items())))] = r
out = []
for r in recs.values():
    try:
        f, x, v = eval(r['f'], ns), eval(r['x'], ns), eval(r['v'], ns)
    except Exception as e:
        out.append(dict(r, n='parse:' + type(e).__name__, p='parse'))
        continue
    fl = r['flags']
    signal.alarm(60)
    try:
        try:
            n = _core._checksol_native(f, x, v, fl.get('numerical', True), fl.get('force', True),
                                       fl.get('simplify', True), fl.get('minimal', False))
        except Timeout:
            n = 'timeout'
        except Exception as e:
            n = 'err:' + type(e).__name__
        t0 = time.perf_counter()
        try:
            p = checksol(f, x, v, **fl)
        except Timeout:
            p = 'timeout'
        except Exception as e:
            p = 'err:' + type(e).__name__
    finally:
        signal.alarm(0)
    out.append(dict(r, n=n, p=p, t=time.perf_counter() - t0))
json.dump(out, open(sys.argv[1], 'w'))
for key in ('n', 'p'):
    c = collections.Counter((str(o['r']), str(o[key])) for o in out)
    print(key, sorted(c.items()))
