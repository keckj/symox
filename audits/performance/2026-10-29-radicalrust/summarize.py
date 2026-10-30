"""Medians per case from raw/<label>-first.jsonl and raw/<label>-repeat.jsonl."""
import json
import os
import statistics
import sys

here = os.path.dirname(os.path.abspath(__file__))
labels = sys.argv[1:] or ['before', 'after', 'sympy']


def load(label, mode):
    out = {}
    with open(os.path.join(here, 'raw', '%s-%s.jsonl' % (label, mode))) as f:
        for line in f:
            if line.startswith('{'):
                r = json.loads(line)
                key = 'repeat_median_s' if mode == 'repeat' else 'first_s'
                out.setdefault(r['case'], [r['expr']]).append(r[key])
    return {k: (v[0], statistics.median(v[1:])) for k, v in out.items()}


first = {lb: load(lb, 'first') for lb in labels}
rep = {lb: load(lb, 'repeat') for lb in labels}
print('| # | call | ' + ' | '.join('%s first (ms)' % lb for lb in labels)
      + ' | ' + ' | '.join('%s repeat (ms)' % lb for lb in labels) + ' |')
print('|---' * (2 + 2 * len(labels)) + '|')
for i in sorted(first[labels[0]]):
    expr = first[labels[0]][i][0]
    cells = ['%.2f' % (1e3 * first[lb][i][1]) for lb in labels]
    cells += ['%.3f' % (1e3 * rep[lb][i][1]) for lb in labels]
    print('| %d | `%s` | %s |' % (i, expr, ' | '.join(cells)))
