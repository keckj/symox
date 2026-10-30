"""Inclusive time and call counts of the prde targets inside risch_integrate.
usage: prof_targets.py <symox|symox-stock|symox-noprecheck|sympy> [nowrap]
symox-noprecheck keeps the new prde entry points and removes only the
splitfactor precheck.
symox-stock is the state before this lot: the three new native entry points
return None and the splitfactor precheck is removed.
With nowrap only the total is measured (no wrappers around the targets).
For symox the native counters (_prde_native.STATS) are printed."""
import sys, time, importlib, signal
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from integrands import INTEGRANDS
B = sys.argv[1]
if B in ('symox-stock', 'symox-noprecheck'):
    from symox.integrals import _prde_native as _nat, _risch_splitfactor_native as _sf
    if B == 'symox-stock':
        for _n in ('is_log_deriv_k_t_radical_in_field', 'parametric_log_deriv_heu', 'limited_integrate'):
            setattr(_nat, _n, lambda *a, **k: None)
    B = 'symox'
    if _sf._ROOT is not None:
        _sf._ROOT['precheck'] = None
S = importlib.import_module(B)
prde = importlib.import_module(B + '.integrals.prde')
risch = importlib.import_module(B + '.integrals.risch')
NAMES = [(prde, 'limited_integrate'), (prde, 'parametric_log_deriv_heu'),
         (prde, 'is_log_deriv_k_t_radical_in_field'), (risch, 'splitfactor'),
         (prde, 'param_rischDE'), (risch, 'residue_reduce')]
stats = {}
depth = {}


def wrap(mod, name):
    f = getattr(mod, name)
    def w(*a, **k):
        d = depth.get(name, 0)
        depth[name] = d + 1
        t0 = time.perf_counter()
        try:
            return f(*a, **k)
        finally:
            depth[name] = d
            if d == 0:
                s = stats.setdefault(name, [0, 0.0, {}])
                s[0] += 1
                s[1] += time.perf_counter() - t0
                DE = a[-1] if name != 'splitfactor' else a[1]
                if name == 'is_log_deriv_k_t_radical_in_field':
                    DE = a[2]
                if name == 'parametric_log_deriv_heu':
                    DE = a[4]
                case = getattr(DE, 'case', '?')
                s[2][case] = s[2].get(case, 0) + 1
    setattr(mod, name, w)


if sys.argv[2:] != ['nowrap']:
    for m, n in NAMES:
        wrap(m, n)
def _h(*a):
    raise TimeoutError
signal.signal(signal.SIGALRM, _h)
x = S.Symbol('x')
total = 0.0
total_ok = 0.0
timed_out = 0
for e in INTEGRANDS:
    expr = S.sympify(e, locals={'x': x})
    signal.alarm(30)
    t0 = time.perf_counter()
    try:
        risch.risch_integrate(expr, x)
    except BaseException as exc:
        timed_out = timed_out + 1 if isinstance(exc, TimeoutError) else timed_out
        print('EXC', e, type(exc).__name__)
        if not isinstance(exc, TimeoutError):
            total_ok += time.perf_counter() - t0
    else:
        total_ok += time.perf_counter() - t0
    signal.alarm(0)
    print('TIME %.6f %s' % (time.perf_counter() - t0, e))
    total += time.perf_counter() - t0
print(B, 'total risch_integrate %.3f s over %d integrands' % (total, len(INTEGRANDS)))
print(B, 'without the %d timed-out integrands (30 s alarm): %.3f s' % (timed_out, total_ok))
for n, (c, t, cases) in sorted(stats.items()):
    print('%-36s calls %5d  inclusive %.3f s  cases %s' % (n, c, t, cases))
if B == 'symox':
    from symox.integrals import _prde_native
    print('STATS', dict(_prde_native.STATS))
