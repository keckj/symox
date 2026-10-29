"""summarize.py raw.jsonl: medians per case and side, markdown table."""
import json
import statistics
import sys
from collections import defaultdict

d = defaultdict(lambda: defaultdict(list))
res = defaultdict(dict)
order = []
for line in open(sys.argv[1]):
    r = json.loads(line)
    if r['data'] is None:
        continue
    c = r['data']['case']
    if c not in order:
        order.append(c)
    d[c][r['side'] + '_f'].append(r['data']['first_ms'])
    d[c][r['side'] + '_r'].append(r['data']['repeat_us'])
    res[c][r['side']] = r['data']['result']
print('| case | first base (ms) | first new (ms) | repeat base (us) | repeat new (us) | base result | new result |')
print('|---|---|---|---|---|---|---|')
for c in order:
    m = {k: statistics.median(v) for k, v in d[c].items()}
    print('| %s | %.2f | %.2f | %.1f | %.1f | `%s` | `%s` |' % (
        c, m['base_f'], m['new_f'], m['base_r'], m['new_r'], res[c]['base'][:60], res[c]['new'][:60]))
