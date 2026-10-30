"""Compare replay outputs: symox (native), symox-stock (Python bodies) and
SymPy 1.14, call by call, with timings per kind.

usage: compare_calls.py CALLS.json SYMPY.json SYMOX.json STOCK.json
"""
import json
import sys

calls = json.load(open(sys.argv[1]))
sy, sx, st = (json.load(open(p))['results'] for p in sys.argv[2:5])
ncorpus = calls['corpus']
summary = {}
diffs = []
for i, (rec, a, b, c) in enumerate(zip(calls['records'], sy, sx, st)):
    kind = rec[0]
    part = 'corpus' if i < ncorpus else 'random'
    key = (kind, part)
    s = summary.setdefault(key, {'n': 0, 'same_sympy': 0, 'same_stock': 0, 'none': 0, 'exc': 0,
                                 't_sympy': 0.0, 't_symox': 0.0, 't_stock': 0.0})
    s['n'] += 1
    s['same_sympy'] += a[1] == b[1]
    s['same_stock'] += c[1] == b[1]
    s['none'] += a[1] is None or (a[1] is not None and a[1][:2] == ['tuple', None])
    s['exc'] += a[1] is not None and a[1][0] == 'EXC'
    s['t_sympy'] += a[2]
    s['t_symox'] += b[2]
    s['t_stock'] += c[2]
    if a[1] != b[1] or c[1] != b[1]:
        diffs.append({'i': i, 'kind': kind, 'args': rec[1], 'sympy': a[1], 'symox': b[1], 'stock': c[1]})
print('%-14s %-7s %5s %10s %10s %6s %5s %10s %10s %10s' % (
    'kind', 'part', 'n', 'eq SymPy', 'eq stock', 'None', 'exc', 'SymPy s', 'native s', 'stock s'))
for (kind, part), s in sorted(summary.items()):
    print('%-14s %-7s %5d %10d %10d %6d %5d %10.3f %10.3f %10.3f' % (
        kind, part, s['n'], s['same_sympy'], s['same_stock'], s['none'], s['exc'],
        s['t_sympy'], s['t_symox'], s['t_stock']))
print('differences:', len(diffs))
for d in diffs[:40]:
    print(json.dumps(d)[:600])
