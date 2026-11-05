"""usage: score.py RUNDIR [OTHER_RUNDIR]

Prints cls/hint/def scores (equal_output + equal_error_observation) of a
replay-ode.py run directory; with a second run, lists rows whose verdict
changed between them."""
import json, sys
from collections import Counter
from pathlib import Path


def rows(d):
    out = {}
    for p in sorted(Path(d).glob('[0-9][0-9].json')):
        for r in json.loads(p.read_text()):
            out[(r['id'], r['method'])] = r
    return out


GOOD = ('equal_output', 'equal_error_observation')


def score(rs):
    c = Counter(m for (i, m), r in rs.items() if r['verdict'] in GOOD)
    return {'cls': c['classify'], 'hint': c['hint'], 'def': c['default'], 'n': len({i for i, m in rs})}


a = rows(sys.argv[1])
print(json.dumps(score(a)))
if len(sys.argv) > 2:
    b = rows(sys.argv[2])
    print(json.dumps(score(b)))
    for k in sorted(a):
        va, vb = a[k]['verdict'], b[k]['verdict']
        if va != vb:
            print(k[0], k[1], va, '->', vb)
