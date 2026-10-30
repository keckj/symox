"""Compare two det_cases.py outputs: python cmp_det.py ref.json other.json [methods]"""
import json, lzma, sys


def load(path):
    opener = lzma.open if path.endswith('.xz') else open
    with opener(path, 'rt') as f:
        return json.load(f)


from collections import Counter
a = load(sys.argv[1]); b = load(sys.argv[2])
meths = sys.argv[3].split(',') if len(sys.argv) > 3 else None
c = Counter(); diffs = []
for k in sorted(a):
    for m in a[k]:
        if meths and m not in meths:
            continue
        if b[k][m][0] == 'SKIPPED':
            c[(m, 'skipped')] += 1; diffs.append((k, m, a[k][m], b[k][m])); continue
        same = a[k][m] == b[k][m]
        c[(m, 'same' if same else 'diff')] += 1
        if not same:
            diffs.append((k, m, a[k][m], b[k][m]))
for kk in sorted(c):
    print(kk, c[kk])
for k, m, x, y in diffs:
    print('--', k, m)
    print('   ref  ', str(x)[:400])
    print('   other', str(y)[:400])
