"""Medians per (config, integrand) from raw.txt; prints markdown tables (first call, repeated call)."""
import collections, statistics, sys
first = collections.defaultdict(list); rep = collections.defaultdict(list); names = {}; hashes = collections.defaultdict(set); texts = collections.defaultdict(set)
for line in open(sys.argv[1] if len(sys.argv) > 1 else 'raw.txt'):
    p = line.split()
    if len(p) < 6 or not p[1].isdigit():
        continue
    label, i = p[0], int(p[1])
    first[(label, i)].append(float(p[2])); rep[(label, i)].append(float(p[3])); names[i] = p[5]; hashes[(label, i)].add(p[4]); texts[(label, i)].add(p[6] if len(p) > 6 else None)
for title, rows in (('first call in a fresh process', first), ('repeated call (median of 3 per process)', rep)):
    print('\n%s\n' % title)
    print('| # | integrand | before ms | after ms | after/before | SymPy ms | after/SymPy | same srepr as before | same srepr as SymPy | same str as SymPy | n |')
    print('|---|---|---|---|---|---|---|---|---|---|---|')
    for i in sorted(names):
        def m(l):
            v = rows.get((l, i)); return statistics.median(v) * 1e3 if v else None
        b, a, s = m('before'), m('after'), m('sympy')
        fmt = lambda v: 'error' if v is None else '%.1f' % v
        rb = '%.2f' % (a / b) if a and b else '-'
        rs = '%.2f' % (a / s) if a and s else '-'
        same = 'yes' if hashes[('before', i)] == hashes[('after', i)] and len(hashes[('after', i)]) == 1 else 'NO'
        ref = 'yes' if hashes[('sympy', i)] == hashes[('after', i)] and len(hashes[('after', i)]) == 1 else 'NO'
        txt = 'yes' if texts[('sympy', i)] == texts[('after', i)] and len(texts[('after', i)]) == 1 else 'NO'
        print('| %d | `%s` | %s | %s | %s | %s | %s | %s | %s | %s | %d |' % (i, names[i], fmt(b), fmt(a), rb, fmt(s), rs, same, ref, txt, len(rows.get(('after', i), []))))
