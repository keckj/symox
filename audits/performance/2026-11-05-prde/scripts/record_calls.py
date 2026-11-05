"""Record the prde calls made by risch_integrate on the integrand corpus and
by the SymPy test_prde.py / test_risch.py copies (symox, Python bodies: the
wrappers defeat the native protocol check, so every call runs the stock
body and is recorded with its DifferentialExtension state).

usage: record_calls.py OUT.json corpus|tests [TESTFILE...]
"""
import json
import os
import signal
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import codec

import symox as S
from symox.integrals import prde, risch, rde

KINDS = ('is_log_deriv_k_t_radical_in_field', 'parametric_log_deriv_heu', 'parametric_log_deriv',
         'limited_integrate', 'limited_integrate_reduce', 'param_rischDE', 'param_poly_rischDE',
         'prde_cancel_liouvillian')
records, seen, stats = [], set(), {}
depth = [0]


def wrap(name):
    f = prde.__dict__[name]

    def w(*a, **k):
        try:
            DE = a[{'is_log_deriv_k_t_radical_in_field': 2, 'parametric_log_deriv_heu': 4,
                    'parametric_log_deriv': 4, 'limited_integrate': 3, 'limited_integrate_reduce': 3,
                    'param_rischDE': 3, 'param_poly_rischDE': 4, 'prde_cancel_liouvillian': 3}[name]]
            args = [codec.value(S, v) for i, v in enumerate(a) if v is not DE]
            kw = {kk: codec.value(S, vv) for kk, vv in k.items()}
            rec = [name, args, kw, codec.de(S, DE), DE.case]
        except Exception as e:  # pragma: no cover
            rec = [name, None, None, None, 'encode-error %r' % e]
        depth[0] += 1
        t0 = time.perf_counter()
        try:
            return f(*a, **k)
        finally:
            dt = time.perf_counter() - t0
            depth[0] -= 1
            key = json.dumps(rec, sort_keys=True)
            st = stats.setdefault((name, rec[4]), [0, 0.0])
            st[0] += 1
            st[1] += dt
            if key not in seen:
                seen.add(key)
                records.append(rec + [dt, SOURCE[0]])
    prde.__dict__[name] = w
    for m in (risch, rde):
        if name in m.__dict__:
            m.__dict__[name] = w


SOURCE = ['?']
for n in KINDS:
    wrap(n)


def _alarm(*a):
    raise TimeoutError


if sys.argv[2] == 'corpus':
    from integrands import INTEGRANDS
    x = S.Symbol('x')
    signal.signal(signal.SIGALRM, _alarm)
    for e in INTEGRANDS:
        SOURCE[0] = e
        signal.alarm(60)
        try:
            risch.risch_integrate(S.sympify(e, locals={'x': x}), x)
        except BaseException:
            pass
        signal.alarm(0)
else:
    import pytest

    class P:
        def pytest_runtest_setup(self, item):
            SOURCE[0] = item.nodeid
    pytest.main(['-q', '-p', 'no:cacheprovider', '-x' if False else '-q'] + sys.argv[3:], plugins=[P()])

with open(sys.argv[1], 'w') as fh:
    json.dump({'records': records,
               'stats': [[k[0], k[1], v[0], v[1]] for k, v in sorted(stats.items())]}, fh)
for k, v in sorted(stats.items()):
    print('%-36s %-10s %6d calls %8.3f s' % (k[0], k[1], v[0], v[1]))
print('distinct records', len(records))
