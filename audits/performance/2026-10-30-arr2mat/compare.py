"""compare.py <sympy.jsonl> <symox.jsonl> [ndiffs]: per-pass agreement with
SymPy 1.14, the number of results computed on the native path and by the
stock fallback, and how many of the stock results agree."""
import json
import sys
from collections import Counter

ref = [json.loads(l) for l in open(sys.argv[1])]
got = [json.loads(l) for l in open(sys.argv[2])]
assert len(ref) == len(got)
agree, total, stock, native, diffs = Counter(), Counter(), Counter(), Counter(), []
for r, g in zip(ref, got):
    assert r['case'] == g['case']
    if r['input'] != g['input']:
        diffs.append((r['case'], 'input', r['input'], g['input']))
        continue
    for k, v in r.items():
        if k in ('case', 'input'):
            continue
        w = g.get(k)
        on_stock = w is not None and w.endswith(' [stock]')
        if on_stock:
            stock[k] += 1
            w = w[:-8]
        total[k] += 1
        if v == w:
            agree[k] += 1
            if not on_stock:
                native[k] += 1
        else:
            diffs.append((r['case'], k, v, w))
print('pass        agree/total  native-path  stock-fallback  stock-agree')
for k in total:
    print('%-11s %4d/%-4d    %4d         %4d            %4d' % (k, agree[k], total[k], total[k] - stock[k], stock[k], agree[k] - native[k]))
print('inputs differing:', sum(1 for d in diffs if d[1] == 'input'))
for d in diffs[:int(sys.argv[3]) if len(sys.argv) > 3 else 40]:
    print('---', d[0], d[1]); print('  sympy:', d[2]); print('  symox:', d[3])
