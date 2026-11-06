"""Three-way comparison of random_param.py runs.

usage: compare_rp.py SYMPY.json NATIVE.json STOCK.json [SHOW]

Each call is classed by comparing the native and stock symox results with
SymPy's: ``same`` (native == SymPy), ``pre-existing`` (native == stock !=
SymPy: the stock symox body differs too), ``native-only`` (stock == SymPy
!= native), ``both-differ`` (all three distinct), and the timeout cases
(``['timeout']`` on any side, not a result). Dummy-count differences are
reported the same way, for calls whose results agree. Times are summed on
the calls where no side timed out.
"""
import collections
import json
import sys


def load(p):
    return {r[0]: r for r in json.load(open(p))}


A, N, T = (load(p) for p in sys.argv[1:4])
show = int(sys.argv[4]) if len(sys.argv) > 4 else 20
cls = collections.Counter()
dcls = collections.Counter()
rows = collections.defaultdict(list)
tm = [0.0, 0.0, 0.0]
native = collections.Counter()
TO = ['timeout']
for k in sorted(A):
    a, n, t = A[k], N[k], T[k]
    ra, rn, rt = a[3], n[3], t[3]
    native[a[1]] += 1
    if n[6]:
        native[a[1] + ' native'] += 1
    if TO in (ra, rn, rt):
        c = 'timeout sympy=%s native=%s stock=%s' % (ra == TO, rn == TO, rt == TO)
    elif rn == ra:
        c = 'same'
    elif rn == rt:
        c = 'pre-existing'
    elif rt == ra:
        c = 'native-only'
    else:
        c = 'both-differ'
    cls[c] += 1
    rows[c].append(k)
    if TO not in (ra, rn, rt):
        tm[0] += a[4]; tm[1] += n[4]; tm[2] += t[4]
        if c == 'same' and a[5] != n[5]:
            dcls['native %d sympy %d stock %d' % (n[5], a[5], t[5])] += 1
            rows['dummies'].append(k)
print('calls', len(A), dict(native))
for c, v in sorted(cls.items()):
    print('%-45s %4d  %s' % (c, v, rows[c][:show]))
print('dummy differences (results equal):', dict(dcls), rows['dummies'][:show])
print('time without timeouts: sympy %.2f s, native %.2f s, stock %.2f s' % tuple(tm))
