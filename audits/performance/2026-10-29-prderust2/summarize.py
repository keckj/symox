"""Median per (config, integrand) from raw.txt; prints a markdown table."""
import collections, statistics, sys
rows = collections.defaultdict(list); names = {}; hashes = collections.defaultdict(set)
for line in open(sys.argv[1] if len(sys.argv) > 1 else 'raw.txt'):
    p = line.split()
    if len(p) < 5 or not p[1].isdigit():
        continue
    label, i, sec, h, name = p[0], int(p[1]), float(p[2]), p[3], p[4]
    rows[(label, i)].append(sec); names[i] = name; hashes[(label, i)].add(h)
print('| # | integrand | before ms | after ms | after/before | SymPy ms | after/SymPy | n |')
print('|---|---|---|---|---|---|---|---|')
for i in sorted(names):
    def m(l):
        v = rows.get((l, i)); return statistics.median(v) * 1e3 if v else None
    b, a, s = m('before'), m('after'), m('sympy')
    fmt = lambda v: 'error' if v is None else '%.1f' % v
    rb = '%.2f' % (a / b) if a and b else '-'
    rs = '%.2f' % (a / s) if a and s else '-'
    n = len(rows.get(('after', i), []))
    print('| %d | `%s` | %s | %s | %s | %s | %s | %d |' % (i, names[i], fmt(b), fmt(a), rb, fmt(s), rs, n))
