"""Median per (case, mode, binary) from raw.jsonl; checks base/new results agree."""
import json, statistics, sys
from collections import defaultdict
rows = [json.loads(l) for l in open(sys.argv[1])]
t, res = defaultdict(list), defaultdict(set)
for r in rows:
    t[r['case'], r['mode'], r['binary']].append(r['seconds'])
    res[r['case'], r['mode']].add(r['result'])
cases = list(dict.fromkeys(r['case'] for r in rows))
print('| case | mode | base ms | new ms | base/new | same result |')
print('|---|---|---:|---:|---:|---|')
for c in cases:
    for m in ('first', 'repeat'):
        b = statistics.median(t[c, m, 'base']) * 1e3
        n = statistics.median(t[c, m, 'new']) * 1e3
        print(f'| {c} | {m} | {b:.3f} | {n:.3f} | {b / n:.2f}x | {len(res[c, m]) == 1} |')
