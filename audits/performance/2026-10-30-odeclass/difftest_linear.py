"""Differential test: native nth_linear_* matches vs the reference matchers.

usage: difftest_linear.py cases.json [extra_eqs.py ...]
Prints one line per disagreement and a summary.
"""
import json, sys, signal, collections
from symox import *  # noqa
from symox.solvers.ode import single as s, _classify_native as cn
from symox.solvers.ode.single import SingleODEProblem
from symox.solvers.deutils import ode_order

CLASSES = [s.NthLinearConstantCoeffHomogeneous, s.NthLinearConstantCoeffVariationOfParameters,
           s.NthLinearConstantCoeffUndeterminedCoefficients, s.NthLinearEulerEqHomogeneous,
           s.NthLinearEulerEqNonhomogeneousVariationOfParameters,
           s.NthLinearEulerEqNonhomogeneousUndeterminedCoefficients]
stats = collections.Counter()


class TO(Exception):
    pass


def alarm(*a):
    raise TO()


signal.signal(signal.SIGALRM, alarm)
SAVED = (cn.nth_linear_match, cn.linear_coefficients)


def off():
    cn.nth_linear_match = lambda solver, eq: None
    cn.linear_coefficients = lambda eq, func, order: cn._MISSING


def on():
    cn.nth_linear_match, cn.linear_coefficients = SAVED


def same(a, b):
    if a is None or b is None:
        return a is b
    return list(a.items()) == list(b.items())


def run(cls, eq, func, x, native):
    (on if native else off)()
    try:
        prob = SingleODEProblem(eq, func, x)
        sol = cls(prob)
        if native:
            src = prob.eq_preprocessed if cls in (s.NthLinearEulerEqHomogeneous, s.NthLinearEulerEqNonhomogeneousVariationOfParameters) else prob.eq_high_order_free
            if cn.nth_linear_match(sol, src) is None:
                return 'decline', None, None
        ok = sol._matches()
        ts = None
        if ok and cls is s.NthLinearConstantCoeffUndeterminedCoefficients:
            ts = set(sol.trialset)
        return ok, (sol.r if ok else None), ts
    finally:
        on()


def check(name, eq, func):
    x = func.args[0]
    if isinstance(eq, Equality):
        eq = eq.lhs - eq.rhs
    try:
        order = ode_order(eq, func)
    except Exception:
        return
    signal.alarm(60)
    try:
        prob = SingleODEProblem(eq, func, x)
        for src in (prob.eq_preprocessed, prob.eq_high_order_free):
            nat = cn.linear_coefficients(src, func, order)
            if nat is cn._MISSING:
                stats['decline coefficients'] += 1
                continue
            off()
            try:
                ref = prob.get_linear_coefficients(src, func, order)
            finally:
                on()
            if same(ref, nat):
                stats['agree coefficients'] += 1
            else:
                stats['DISAGREE'] += 1
                print('DISAGREE', name, 'coefficients', ref, nat, flush=True)
    except TO:
        stats['timeout'] += 1
    except Exception as exc:
        stats['error'] += 1
        print('ERROR', name, 'coefficients', type(exc).__name__, str(exc)[:100], flush=True)
    finally:
        signal.alarm(0)
    for cls in CLASSES:
        signal.alarm(60)
        try:
            n = run(cls, eq, func, x, True)
            if n[0] == 'decline':
                stats['decline ' + cls.hint] += 1
                continue
            r = run(cls, eq, func, x, False)
        except TO:
            stats['timeout'] += 1
            continue
        except Exception as exc:
            stats['error'] += 1
            print('ERROR', name, cls.hint, type(exc).__name__, str(exc)[:100], flush=True)
            continue
        finally:
            signal.alarm(0)
        if bool(r[0]) == bool(n[0]) and same(r[1], n[1]) and r[2] == n[2]:
            stats[('match ' if n[0] else 'agree ') + cls.hint] += 1
        else:
            stats['DISAGREE'] += 1
            print('DISAGREE', name, cls.hint, 'ref', r, 'native', n, flush=True)


def main():
    cases = json.load(open(sys.argv[1]))
    for c in cases:
        try:
            import symox as M
            ns = {n: getattr(M, n) for n in dir(M) if not n.startswith('_')}
            eq = eval(c['eq'], ns)
            func = eval(c['func'], ns)
        except Exception:
            continue
        check(c['name'], eq, func)
    for path in sys.argv[2:]:
        ns = {}
        exec(open(path).read(), ns)
        for i, (eq, func) in enumerate(ns['EQS']):
            check('%s:%d' % (path, i), eq, func)
    tot = collections.Counter()
    for k, v in stats.items():
        tot[k.split(' ')[0]] += v
    print(dict(tot))
    print(sorted(stats.items()))


main()
