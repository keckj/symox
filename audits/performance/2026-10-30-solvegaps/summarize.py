"""summarize.py <outdir>: median over the runs of each label, as Markdown."""
import glob
import json
import os
import statistics
import sys

out = sys.argv[1]
labels = ('sympy', 'before', 'after')
data = {}
for lab in labels:
    runs = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(out, lab + '-*.json')))]
    data[lab] = runs
print('| call | SymPy cold / warm (ms) | before cold / warm (ms) | after cold / warm (ms) | after == SymPy |')
print('|---|---|---|---|---|')
for i, case in enumerate(data['sympy'][0]):
    cells = []
    for lab in labels:
        cold = statistics.median(r[i]['cold'] for r in data[lab]) * 1e3
        warm = statistics.median(r[i]['warm'] for r in data[lab]) * 1e3
        cells.append('%.2f / %.2f' % (cold, warm))
    same = data['after'][0][i]['result'] == case['result']
    print('| %s | %s | %s |' % (case['name'], ' | '.join(cells), 'yes' if same else 'no'))
print()
for lab in ('before', 'after'):
    for i, case in enumerate(data['sympy'][0]):
        r = data[lab][0][i]['result']
        if r != case['result']:
            print('- %s, %s: `%s` (SymPy `%s`)' % (lab, case['name'], r[:160], case['result'][:160]))
