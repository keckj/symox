"""Median first-call / repeat times per case from raw_{sympy,before,after}.tsv."""
import statistics
import sys
from pathlib import Path

here = Path(__file__).resolve().parent


def load(name):
    rows = {}
    for line in (here / name).read_text().splitlines():
        parts = line.split('\t')
        if len(parts) < 5:
            rows.setdefault(parts[1], []).append(None)
            continue
        rows.setdefault(parts[1], []).append((float(parts[2]), float(parts[3]), parts[4]))
    return rows


def med(rows, i):
    vals = [r[i] for r in rows if r]
    return statistics.median(vals) if vals else None


data = {k: load('raw_%s.tsv' % k) for k in ('sympy', 'before', 'after')}
print('| call | SymPy 1.14 first | before first | after first | after repeat | after/SymPy | before state |')
print('|---|---|---|---|---|---|---|')
for case in data['after']:
    s, b, a = data['sympy'].get(case, []), data['before'][case], data['after'][case]
    sf, bf, af, ar = med(s, 0), med(b, 0), med(a, 0), med(a, 1)
    bstate = 'unevaluated' if b and b[0] and b[0][2] == 'UNEVAL' else 'exact'
    if a and a[0] and s and s[0] and a[0][2] != s[0][2]:
        bstate += ' MISMATCH'
    fmt = lambda v: '-' if v is None else '%.4f s' % v
    ratio = '-' if not sf or not af else '%.3g' % (af / sf)
    print('| `%s` | %s | %s | %s | %s | %s | %s |' % (case, fmt(sf), fmt(bf), fmt(af), fmt(ar), ratio, bstate))
