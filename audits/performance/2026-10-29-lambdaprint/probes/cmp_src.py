import json, sys
a = json.load(open(sys.argv[1])); b = json.load(open(sys.argv[2]))
n = 0; bad = {'src': 0, 'doc': 0, 'sig': 0, 'name': 0, 'err': 0}
verbose = len(sys.argv) > 3
for k in sorted(a):
    n += 1
    x, y = a[k], b.get(k)
    if not isinstance(x, dict) or not isinstance(y, dict):
        if x != y:
            bad['err'] += 1
            if verbose: print('ERR', k, repr(x)[:100], '|', repr(y)[:100])
        continue
    for f in ('src', 'doc', 'sig', 'name'):
        if x[f] != y[f]:
            bad[f] += 1
            if verbose and f in sys.argv[3]: print(f.upper(), k, '\n  S:', repr(x[f])[:300], '\n  X:', repr(y[f])[:300])
print('total', n, bad)
