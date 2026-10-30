"""compare.py TAG [TAG2]: statement harness results of symox TAG (and TAG2)
against SymPy 1.14, per file: number of evaluated assertion sides whose str
equals the SymPy one; with TAG2, the keys that changed status."""
import glob, json, os, sys

raw = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'raw', 'tl')
tags = sys.argv[1:]
tot = {t: [0, 0] for t in tags}
changes = []
for f in sorted(glob.glob(os.path.join(raw, 'sympy_*.json'))):
    name = os.path.basename(f)[6:-5]
    ref = json.load(open(f))
    row = [name, str(len(ref))]
    res = {}
    for t in tags:
        p = os.path.join(raw, '%s_%s.json' % (t, name))
        d = json.load(open(p)) if os.path.exists(p) else {}
        res[t] = d
        same = sum(1 for k, v in ref.items() if k in d and d[k][0] == v[0])
        row.append(str(same))
        tot[t][0] += same
        tot[t][1] += len(ref)
    print('  '.join(row))
    if len(tags) == 2:
        a, b = res[tags[0]], res[tags[1]]
        for k, v in sorted(ref.items()):
            sa = k in a and a[k][0] == v[0]
            sb = k in b and b[k][0] == v[0]
            if sa != sb:
                changes.append(('FIXED' if sb else 'BROKEN', k, v[0][:60], a.get(k, ['-'])[0][:60], b.get(k, ['-'])[0][:60]))
print('total', {t: '%d/%d' % tuple(v) for t, v in tot.items()})
for c in changes:
    print(' | '.join(c))
