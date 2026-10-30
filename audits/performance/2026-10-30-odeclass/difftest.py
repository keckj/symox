"""Differential test: native first-order matches vs the reference matcher.

usage: difftest.py cases.json [extra_eqs.py]
Prints one line per disagreement and a summary.
"""
import json, sys, signal, collections
from symox import *  # noqa
from symox.solvers.ode import single as s, _classify_native as cn
from symox.solvers.ode.single import SingleODEProblem
from symox.solvers.deutils import ode_order

CLASSES = [s.Separable, s.FirstExact, s.FirstLinear, s.Bernoulli, s.RiccatiSpecial,
           s.HomogeneousCoeffBest, s.HomogeneousCoeffSubsIndepDivDep,
           s.HomogeneousCoeffSubsDepDivIndep, s.AlmostLinear, s.LinearCoefficients,
           s.SeparableReduced]

stats = collections.Counter()


class TO(Exception):
    pass


def alarm(*a):
    raise TO()


signal.signal(signal.SIGALRM, alarm)


def reference(cls, eq, func, x):
    prob = SingleODEProblem(eq, func, x)
    saved = cn.solver_match
    cn.solver_match = lambda solver: None
    try:
        sol = cls(prob)
        ok = sol._matches()
        return ok, getattr(sol, '_wilds_match', None), [sol._wilds_match.get(w, S.Zero) for w in sol.wilds()] if getattr(sol, '_wilds_match', None) is not None else None
    finally:
        cn.solver_match = saved


def native(cls, eq, func, x):
    prob = SingleODEProblem(eq, func, x)
    sol = cls(prob)
    m = cn.solver_match(sol)
    if m is None:
        return 'decline', None, None
    if m is False:
        return False, None, None
    sol._wilds_match = m
    ok = sol._verify(func)
    return ok, m, [sol._wilds_match.get(w, S.Zero) for w in sol.wilds()]


def check(name, eq, func):
    x = func.args[0]
    if isinstance(eq, Equality):
        eq = eq.lhs - eq.rhs
    try:
        if ode_order(eq, func) != 1:
            return
    except Exception:
        return
    signal.alarm(60)
    try:
        from symox.solvers.ode._ode_helpers import _preprocess
        df = func.diff(x)
        pe = expand(_preprocess(eq, func)[0])
        col = collect(pe, df, exact=True)
        d = Wild('d', exclude=[df, func.diff(x, 2)]); e = Wild('e', exclude=[df])
        nat = cn.split_de(col, df, d, e, func.diff(x, 2))
        if nat is cn._MISSING:
            stats['decline power_series'] += 1
        else:
            ref = col.match(d + e*df)
            if ref == nat:
                stats['agree power_series'] += 1
            else:
                stats['DISAGREE'] += 1
                print('DISAGREE', name, 'power_series', ref, nat, flush=True)
    except TO:
        stats['timeout'] += 1
    except Exception as exc:
        stats['error'] += 1
        print('ERROR', name, 'power_series', type(exc).__name__, str(exc)[:100], flush=True)
    finally:
        signal.alarm(0)
    signal.alarm(60)
    try:
        from symox.solvers.ode._ode_helpers import _preprocess
        from symox.solvers.ode.riccati import match_riccati
        pe = SingleODEProblem(eq, func, x).eq_expanded
        nat = cn.riccati_match(pe, func.func, x)
        if nat is None:
            stats['decline riccati'] += 1
        else:
            ref = match_riccati(pe, func.func, x)
            if ref[0] == nat[0] and list(ref[1]) == list(nat[1]):
                stats['agree riccati'] += 1
            else:
                stats['DISAGREE'] += 1
                print('DISAGREE', name, 'riccati', ref, nat, flush=True)
    except TO:
        stats['timeout'] += 1
    except Exception as exc:
        stats['error'] += 1
        print('ERROR', name, 'riccati', type(exc).__name__, str(exc)[:100], flush=True)
    finally:
        signal.alarm(0)
    for cls in CLASSES:
        signal.alarm(60)
        try:
            n = native(cls, eq, func, x)
            if n[0] == 'decline':
                stats['decline ' + cls.hint] += 1
                signal.alarm(0)
                continue
            r = reference(cls, eq, func, x)
        except TO:
            stats['timeout'] += 1
            continue
        except Exception as exc:
            stats['error'] += 1
            print('ERROR', name, cls.hint, type(exc).__name__, str(exc)[:100], flush=True)
            continue
        finally:
            signal.alarm(0)
        rm = r[1]
        nm = n[1]
        same_dict = (rm is None and nm is None) or (rm is not None and nm is not None and dict(rm) == dict(nm))
        if r[0] == n[0] and same_dict and r[2] == n[2]:
            stats['agree ' + cls.hint] += 1
        else:
            stats['DISAGREE'] += 1
            print('DISAGREE', name, cls.hint, 'ref', r[0], rm, 'native', n[0], nm, flush=True)


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
