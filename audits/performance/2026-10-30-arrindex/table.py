"""Markdown table from raw/{sympy,before,after}.jsonl: median of the
per-process medians, in ms."""
import json
import statistics
import sys

d = sys.argv[1] if len(sys.argv) > 1 else 'raw'
rows = {}
for mode in ('sympy', 'before', 'after'):
    for line in open('%s/%s.jsonl' % (d, mode)):
        r = json.loads(line)
        rows.setdefault(r['case'], {}).setdefault(mode, []).append(r)


def med(rs):
    v = [r['median_ms'] for r in rs if r['median_ms'] is not None]
    return statistics.median(v) if v else None


def fmt(v, rs):
    return '%.4f' % v if v is not None else rs[0]['result'].split(':')[0]


print('| # | call | SymPy 1.14 | before | after | SymPy/after | before/after | result vs SymPy |')
print('|---|---|---|---|---|---|---|---|')
for n, (case, modes) in enumerate(rows.items()):
    s, b, a = (med(modes[m]) for m in ('sympy', 'before', 'after'))
    same = 'same' if modes['sympy'][0]['result'] == modes['after'][0]['result'] else 'differs'
    ratio = lambda x, y: '%.1fx' % (x/y) if x and y else '-'
    print('| %d | `%s` | %s | %s | %s | %s | %s | %s |' % (n, case, fmt(s, modes['sympy']), fmt(b, modes['before']),
                                                     fmt(a, modes['after']), ratio(s, a), ratio(b, a), same))
