"""usage: compare.py REF.jsonl A.jsonl [B.jsonl] [-v]

Counts, per case family and field (hint_sol, nosimp, default), the cases
where A equals REF; with B, also where A equals B and the cases where A and
B disagree with each other.
"""
import collections
import json
import sys

args = [a for a in sys.argv[1:] if not a.startswith('-')]
verbose = '-v' in sys.argv


def load(p):
    return {j['name']: j for j in map(json.loads, open(p))}


ref = load(args[0])
a = load(args[1])
b = load(args[2]) if len(args) > 2 else None
FIELDS = ('hint_sol', 'nosimp', 'default')
same_ref = collections.Counter()
total = collections.Counter()
same_b = collections.Counter()
for name, r in ref.items():
    if name not in a or (b is not None and name not in b):
        continue
    fam = name.split('_')[0]
    for fld in [g for g in FIELDS if g in r]:
        total[fam, fld] += 1
        if a[name][fld] == r[fld]:
            same_ref[fam, fld] += 1
        elif verbose:
            print('REF', name, fld, '\n  R', r[fld][:300], '\n  A', a[name][fld][:300])
        if b is not None:
            if a[name][fld] == b[name][fld]:
                same_b[fam, fld] += 1
            elif verbose:
                print('B', name, fld, '\n  A', a[name][fld][:300], '\n  B', b[name][fld][:300])
for key in sorted(total):
    line = '%-6s %-9s A==REF %3d/%3d' % (key[0], key[1], same_ref[key], total[key])
    if b is not None:
        line += '   A==B %3d/%3d' % (same_b[key], total[key])
    print(line)
print('all A==REF %d/%d' % (sum(same_ref.values()), sum(total.values())))
if b is not None:
    print('all A==B %d/%d' % (sum(same_b.values()), sum(total.values())))
