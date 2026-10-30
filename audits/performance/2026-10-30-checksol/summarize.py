"""Medians over the bench_<label>_<i>.json runs: summarize.py raw"""
import glob
import json
import statistics
import sys

d = sys.argv[1]
labels = ('sympy', 'before', 'after')
runs = {lab: [json.load(open(p)) for p in sorted(glob.glob('%s/bench_%s_*.json' % (d, lab)))] for lab in labels}
print('| call | SymPy cold / warm (ms) | before cold / warm (ms) | after cold / warm (ms) | before r | after r | SymPy r |')
print('|---|---|---|---|---|---|---|')
for name in runs['sympy'][0]:
    cells = []
    for lab in labels:
        c = statistics.median(r[name]['cold_ms'] for r in runs[lab])
        w = statistics.median(r[name]['warm_ms'] for r in runs[lab])
        cells.append('%.2f / %.2f' % (c, w))
    rs = [runs[lab][0][name]['r'] for lab in ('before', 'after', 'sympy')]
    print('| %s | %s | %s | %s | %s | %s | %s |' % (name, *cells, *rs))
