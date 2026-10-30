"""Per-case comparison of two replay-ode.py output directories (the corpus-r11
script): python ode_compare.py <base dir> <after dir>.

Prints, per method, the verdict counts of each run and the cases whose symox
observation (output text or error class) differs between the two runs.
Resource-limited observations are reported separately: they are inconclusive.
"""
import glob
import json
import sys
from collections import Counter


def load(d):
    out = {}
    for p in sorted(glob.glob(d + '/[0-9][0-9].json')):
        for r in json.load(open(p)):
            out[(r['id'], r['method'])] = r
    return out


def obs(r):
    s = r.get('symox') or {}
    if s.get('ok'):
        return ('ok', s.get('text'))
    return ('err', s.get('error'))


a, b = load(sys.argv[1]), load(sys.argv[2])
assert set(a) == set(b), 'different case sets'
for run, d in (('base', a), ('after', b)):
    print(run, {m: dict(Counter(r['verdict'] for k, r in d.items() if k[1] == m))
                for m in ('classify', 'hint', 'default')})
diff = Counter(); limited = Counter()
for k in sorted(a):
    ra, rb = a[k], b[k]
    if 'resource' in ra['verdict'] or 'resource' in rb['verdict']:
        if ra['verdict'] != rb['verdict']:
            limited[k[1]] += 1
            print('inconclusive', k, ra['verdict'], rb['verdict'])
        continue
    if obs(ra) != obs(rb):
        diff[k[1]] += 1
        print('DIFF', k, str(obs(ra))[:150], '|', str(obs(rb))[:150])
print('per-case symox differences:', dict(diff), 'resource-limited changes:', dict(limited))
