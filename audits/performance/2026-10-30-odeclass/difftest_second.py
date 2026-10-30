"""Differential test: native second-order matches vs the reference matchers
(2nd_linear_airy, 2nd_linear_bessel, match_2nd_hypergeometric and the
order-two power series hints of classify_ode).

usage: difftest_second.py cases.json [n_classify]
Prints one line per disagreement and a summary.
"""
import collections
import itertools
import json
import signal
import sys

from symox import *  # noqa
from symox.solvers.ode import single as s, hypergeometric as hg, _classify_native as cn
from symox.solvers.ode.single import SingleODEProblem
from symox.solvers.deutils import ode_order

globals().pop("test", None)

stats = collections.Counter()
x, y, a, b = symbols('x y a b')
f = Function('f')
fx = f(x)
SAVED = (cn.airy_match, cn.bessel_rejected, cn.match_2nd_hypergeometric, cn.power_series2)


class TO(Exception):
    pass


def alarm(*args):
    raise TO()


signal.signal(signal.SIGALRM, alarm)


def off():
    cn.airy_match = lambda solver, eq: None
    cn.bessel_rejected = lambda solver, eq: False
    cn.match_2nd_hypergeometric = lambda eq, func: cn._MISSING
    cn.power_series2 = lambda *args: cn._MISSING


def on():
    cn.airy_match, cn.bessel_rejected, cn.match_2nd_hypergeometric, cn.power_series2 = SAVED


def guarded(fn, *args, **kw):
    signal.alarm(60)
    try:
        return 'ok', fn(*args, **kw)
    except TO:
        return 'timeout', None
    except Exception as exc:
        return 'error', type(exc).__name__
    finally:
        signal.alarm(0)


def ref(fn, *args, **kw):
    off()
    try:
        return guarded(fn, *args, **kw)
    finally:
        on()


def compare(name, what, decided, native, reference):
    if not decided:
        stats['decline ' + what] += 1
        return
    if reference[0] == 'timeout':
        stats['timeout'] += 1
        return
    if native == reference:
        stats[('match ' if native[1] else 'agree ') + what] += 1
    else:
        stats['DISAGREE'] += 1
        print('DISAGREE', name, what, 'ref', reference, 'native', native, flush=True)


def airy(eq, func, x):
    sol = s.SecondLinearAiry(SingleODEProblem(eq, func, x))
    ok = sol._matches()
    return ok, (sol.rn if ok else None)


def bessel(eq, func, x):
    return s.SecondLinearBessel(SingleODEProblem(eq, func, x))._matches()


def check(name, eq, func, classify=False):
    x = func.args[0]
    if isinstance(eq, Equality):
        eq = eq.lhs - eq.rhs
    status, prob = guarded(SingleODEProblem, eq, func, x)
    if status != 'ok':
        return
    hof = prob.eq_high_order_free
    # Airy
    d = guarded(cn.airy_match, s.SecondLinearAiry(prob), hof)
    decided = d[0] == 'ok' and d[1] is not None
    compare(name, 'airy', decided, guarded(airy, eq, func, x) if decided else None, ref(airy, eq, func, x) if decided else None)
    # Bessel: the native path only rejects
    d = guarded(cn.bessel_rejected, s.SecondLinearBessel(prob), hof)
    decided = d == ('ok', True)
    compare(name, 'bessel', decided, ('ok', False) if decided else None, ref(bessel, eq, func, x) if decided else None)
    # 2nd_hypergeometric
    pre = prob.eq_preprocessed
    d = guarded(cn.match_2nd_hypergeometric, pre, func)
    decided = d[0] == 'ok' and d[1] is not cn._MISSING
    compare(name, 'hypergeometric', decided, d if decided else None, ref(hg.match_2nd_hypergeometric, pre, func) if decided else None)
    if classify:
        for x0 in (0, 1, -1, y):
            nat = guarded(classify_ode, eq, func, dict=True, x0=x0)
            refv = ref(classify_ode, eq, func, dict=True, x0=x0)
            if nat[0] != 'ok' or refv[0] != 'ok':
                compare(name, 'classify x0=%s' % x0, True, nat[0], refv[0])
                continue
            keys = lambda m: [k for k in m if k != 'order' and k != 'default' and k != 'ordered_hints']
            series = lambda m: {k: m[k] for k in m if str(k).startswith('2nd_power_series')}
            same = keys(nat[1]) == keys(refv[1]) and series(nat[1]) == series(refv[1]) and \
                nat[1].get('ordered_hints') == refv[1].get('ordered_hints')
            if same:
                stats[('match ' if series(nat[1]) else 'agree ') + 'classify'] += 1
            else:
                stats['DISAGREE'] += 1
                print('DISAGREE', name, 'classify x0=%s' % x0, 'ref', keys(refv[1]), series(refv[1]),
                      'native', keys(nat[1]), series(nat[1]), flush=True)


A = ["1", "x", "x**2", "x**2 - 1", "x*(1 - x)", "x**3", "2", "y", "(x - 1)**2", "1 + x**2",
     "exp(x)", "1.5", "x - 2", "sin(x)", "3*x**2 + 2*x", "a*x**2", "-x", "I", "0"]
B = ["0", "1", "x", "1/x", "3 - 5*x", "-2*x", "x**2", "y", "2*x - 1", "exp(x)", "1.0", "a + b*x"]
C = ["0", "1", "-x", "x", "-(2 + 3*x)", "x**2 - 4", "x**2 - Rational(1, 4)", "-6", "a*b", "-a*b",
     "y*x", "x**4", "log(x)", "1/x**2", "-x**3", "x - y", "0.0", "besselj(0, x)", "-(a + b*x)"]
R = ["0", "1", "exp(x)", "x"]


def generated():
    for i, (p, q, r) in enumerate(itertools.product(A, B, C)):
        rhs = R[i % len(R)]
        src = "(%s)*d2 + (%s)*d1 + (%s)*f(x) - (%s)" % (p, q, r, rhs)
        eq = eval(src, {'d2': fx.diff(x, 2), 'd1': fx.diff(x), 'f': f, 'x': x, 'y': y, 'a': a, 'b': b,
                        'exp': exp, 'sin': sin, 'log': log, 'I': I, 'Rational': Rational, 'besselj': besselj})
        yield src, eq


EXTRA = ["x**2*d2 + x*d1 + (x**2 - 4)*f(x)", "d2 - x*f(x)", "d2 + (2 + 3*x)*f(x)", "x*(1 - x)*d2 + (1 - 3*x)*d1 - f(x)",
         "(x - 1)**2*d2 + (x - 1)*d1 + ((x - 1)**2 - 1)*f(x)", "x**2*d2 + 3*x*d1 + (4*x**4 - 5)*f(x)",
         "d2 + x*d1 + f(x)", "x*d2 - d1 + x**3*f(x)", "d2*x**2 + x*d1 - f(x)*x - f(x)", "d3 + d2 + f(x)",
         "d1 + f(x)", "d2 + d1**2 + f(x)", "d2 + f(x)**2", "(x**2 + 1)*d2 + 2*x*d1", "d2 + f(x)/x", "d2 + d1/x + f(x)/x**2",
         "x**2*d2 - x*d1 + f(x)", "(1 - x**2)*d2 - 2*x*d1 + 6*f(x)", "d2 + sin(x)*f(x)", "d2 + exp(x)*f(x)",
         "4*x**2*d2 + 8*x*d1 + (x - 3)*f(x)", "2*x*d2 + d1 - 2*f(x)", "x*d2 + 2*d1 + x*f(x)"]


def main():
    nclass = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    cases = json.load(open(sys.argv[1]))
    import symox as M
    ns = {n: getattr(M, n) for n in dir(M) if not n.startswith('_')}
    for c in cases:
        try:
            eq = eval(c['eq'], ns)
            func = eval(c['func'], ns)
            order = ode_order(eq.lhs - eq.rhs if isinstance(eq, Equality) else eq, func)
        except Exception:
            continue
        check(c['name'], eq, func, classify=(order == 2))
    env = {'d3': fx.diff(x, 3), 'd2': fx.diff(x, 2), 'd1': fx.diff(x), 'f': f, 'x': x, 'sin': sin, 'exp': exp}
    for src in EXTRA:
        check(src, eval(src, env), fx, classify=True)
    for i, (src, eq) in enumerate(generated()):
        check(src, eq, fx, classify=(nclass and i % nclass == 0))
    tot = collections.Counter()
    for k, v in stats.items():
        tot[k.split(' ')[0]] += v
    print(dict(tot))
    print(sorted(stats.items()))


main()
