"""Print a markdown table from raw-{base,new,sympy}.jsonl."""
import json, os, sys
here = os.path.dirname(os.path.abspath(__file__))
rows = {}
for label in ('base', 'new', 'sympy'):
    path = os.path.join(here, 'raw-%s.jsonl' % label)
    if not os.path.exists(path):
        continue
    for line in open(path):
        d = json.loads(line)['data']
        rows.setdefault(d['case'], {})[label] = d
def fmt(d, key):
    if d is None:
        return '-'
    if 'timeout' in d:
        return '>%ss' % d['timeout']
    if str(d.get('result', '')).startswith('ERROR'):
        return 'error'
    v = d.get(key)
    return '-' if v is None else '%.3f' % v
print('| case | base first | new first | sympy first | base rep | new rep | sympy rep | base/new | sympy/new | same result as sympy |')
print('|---|---|---|---|---|---|---|---|---|---|')
for case, r in rows.items():
    b, n, s = r.get('base'), r.get('new'), r.get('sympy')
    def ratio(a, c):
        try:
            if any(str(d.get('result', '')).startswith('ERROR') or 'timeout' in d for d in (a, c)):
                return '-'
            return '%.1fx' % (a['first'] / c['first'])
        except Exception:
            return '-'
    same = '-' if not (n and s) else str(n.get('result') == s.get('result'))
    print('| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |' % (
        case, fmt(b, 'first'), fmt(n, 'first'), fmt(s, 'first'), fmt(b, 'repeat_median'),
        fmt(n, 'repeat_median'), fmt(s, 'repeat_median'), ratio(b, n), ratio(s, n), same))
