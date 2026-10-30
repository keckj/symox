"""Print the README table from sympy.json, before.json and after.json."""
import json
s, b, a = (json.load(open(f)) for f in ('sympy.json', 'before.json', 'after.json'))
fmt = lambda v: 'error' if v is None else '%.3g' % v
print('| case | SymPy 1.14 | before | after | after/SymPy | first call SymPy / after |')
print('|---|---|---|---|---|---|')
for k in a:
    sa, sb, sf = (d[k]['repeat_median'] for d in (s, b, a))
    print('| %s | %s | %s | %s | %.2f | %.3f / %.3f |' % (k, fmt(sa), fmt(sb), fmt(sf), sf / sa,
          s[k]['first_median'], a[k]['first_median']))
