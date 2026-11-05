"""Compare encoded prde results of two backends.

usage: compare.py REFERENCE.json OTHER.json [SHOW]

Rows are matched on their index. Prints the number of compared calls per
kind, the result and Dummy-count differences (the first SHOW of them) and
the total times. Works on the output of replay.py and random_towers.py
(``.json`` or ``.json.gz``).
"""
import collections
import gzip
import json
import sys


def load(path):
    rows = json.load(gzip.open(path, 'rt') if path.endswith('.gz') else open(path))
    # random_towers: [i, kind, spec, result, s, dummies, native];
    # replay: [i, kind, result, s, dummies].
    tower = bool(rows) and isinstance(rows[0][2], dict)
    out = {}
    for r in rows:
        if tower:
            out[r[0]] = (r[1], r[3], r[4], r[5] if len(r) > 5 else None, r[2])
        else:
            out[r[0]] = (r[1], r[2], r[3], r[4] if len(r) > 4 else None, None)
    return out


A, B = load(sys.argv[1]), load(sys.argv[2])
show = int(sys.argv[3]) if len(sys.argv) > 3 else 10
kinds = collections.Counter()
diffs = []
for k in sorted(A):
    if k not in B:
        continue
    kind, ra, ta, da, spec = A[k]
    _, rb, tb, db, _ = B[k]
    kinds[kind] += 1
    if ra != rb:
        diffs.append(('result', k, kind, ra, rb, spec))
    elif da is not None and db is not None and da != db:
        diffs.append(('dummies', k, kind, da, db, spec))
print('compared', sum(kinds.values()), dict(kinds))
print('differences', dict(collections.Counter(d[0] for d in diffs)))
for d in diffs[:show]:
    print('  ', str(d)[:400])
common = [k for k in A if k in B]
print('time %.3f s / %.3f s' % (sum(A[k][2] for k in common), sum(B[k][2] for k in common)))
