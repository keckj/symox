"""Median over processes of first-call and repeated-call times per label/case."""
import json, statistics, sys
from collections import defaultdict
d = defaultdict(lambda: defaultdict(list)); res = {}
for line in open(sys.argv[1] if len(sys.argv) > 1 else 'raw.jsonl'):
    label, js = line.split(' ', 1); r = json.loads(js)
    if r.get('timeout'):
        d[label][r['case']].append(None); continue
    d[label][r['case']].append((r['first'], r['median'])); res[(label, r['case'])] = r['result'][:50]
def med(v, i):
    v = [x[i] for x in v if x]
    return statistics.median(v) if v else float('nan')
print('| case | base first | after first | ratio | base repeat | after repeat | ratio | sympy repeat | after/sympy |')
print('|---|---|---|---|---|---|---|---|---|')
for c in d['after']:
    b, a, s = d['base'][c], d['after'][c], d['sympy'][c]
    bf, af, br, ar, sr = med(b, 0), med(a, 0), med(b, 1), med(a, 1), med(s, 1)
    print('| %s | %.4f | %.4f | %.1fx | %.4f | %.4f | %.1fx | %.4f | %.2f |' % (c, bf, af, bf / af, br, ar, br / ar, sr, ar / sr))
print()
for c in d['after']:
    if res.get(('base', c)) != res.get(('after', c)):
        print('result differs: %s base=%r after=%r sympy=%r' % (c, res.get(('base', c)), res.get(('after', c)), res.get(('sympy', c))))
