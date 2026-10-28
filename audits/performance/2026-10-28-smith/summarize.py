"""Medians per case and side from raw.jsonl; checks base and new results are equal."""
import json, statistics, sys
from collections import defaultdict
rows = [json.loads(l) for l in open(sys.argv[1])]
d = defaultdict(lambda: defaultdict(list)); res = defaultdict(dict)
for r in rows:
    d[r['case']][r['side']].append(r)
    res[r['case']].setdefault(r['side'], set()).add(r['result'])
print('| case | first base (ms) | first new (ms) | ratio | repeat base (us) | repeat new (us) | ratio | same result |')
print('|---|---|---|---|---|---|---|---|')
for case, sides in d.items():
    fb = statistics.median(x['first_s'] for x in sides['base']) * 1e3
    fn = statistics.median(x['first_s'] for x in sides['new']) * 1e3
    rb = statistics.median(x['repeat_s'] for x in sides['base']) * 1e6
    rn = statistics.median(x['repeat_s'] for x in sides['new']) * 1e6
    same = res[case]['base'] == res[case]['new'] and len(res[case]['base']) == 1
    print(f'| {case} | {fb:.3f} | {fn:.3f} | {fb/fn:.2f}x | {rb:.1f} | {rn:.1f} | {rb/rn:.2f}x | {same} |')
