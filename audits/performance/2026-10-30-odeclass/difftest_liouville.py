"""Differential test: native Liouville and 2nd_nonlinear_autonomous_conserved
matches vs the reference matcher (the native hook disabled).

usage: difftest_liouville.py cases.json
Prints one line per disagreement and a summary.
"""
import collections
import itertools
import json
import signal
import sys

from symox import *  # noqa
from symox.solvers.ode import single as s, _classify_native as cn
from symox.solvers.ode.single import SingleODEProblem
from symox.solvers.deutils import ode_order

globals().pop("test", None)

CLASSES = [s.Liouville, s.SecondNonlinearAutonomousConserved]
stats = collections.Counter()
x, y, a = symbols('x y a')
f, g = Function('f'), Function('g')
fx = f(x)
SAVED = cn.solver_match


class TO(Exception):
    pass


def alarm(*args):
    raise TO()


signal.signal(signal.SIGALRM, alarm)


def run(cls, eq, func, x, native):
    sol = cls(SingleODEProblem(eq, func, x))
    if native:
        m = cn.solver_match(sol)
        if m is None:
            return 'decline'
    else:
        cn.solver_match = lambda solver: None
    try:
        ok = sol._matches()
    finally:
        cn.solver_match = SAVED
    wm = getattr(sol, '_wilds_match', None)
    return ok, (dict(wm) if wm is not None else None), (sol.wilds_match() if wm is not None else None)


def check(name, eq, func):
    if isinstance(eq, Equality):
        eq = eq.lhs - eq.rhs
    if len(getattr(func, 'args', ())) != 1:
        return
    x = func.args[0]
    try:
        if ode_order(eq, func) != 2:
            return
    except Exception:
        return
    for cls in CLASSES:
        signal.alarm(60)
        try:
            nat = run(cls, eq, func, x, True)
            if nat == 'decline':
                stats['decline ' + cls.hint] += 1
                continue
            ref = run(cls, eq, func, x, False)
            if nat == ref:
                stats[('match ' if nat[1] is not None else 'agree ') + cls.hint] += 1
            else:
                stats['DISAGREE'] += 1
                print('DISAGREE', name, cls.hint, 'ref', ref, 'native', nat, flush=True)
        except TO:
            stats['timeout'] += 1
        except Exception as exc:
            stats['error'] += 1
            print('ERROR', name, cls.hint, type(exc).__name__, str(exc)[:100], flush=True)
        finally:
            signal.alarm(0)


TERMS = ["d2", "x*d2", "fx*d2", "2*d2", "d1**2", "d1**2/fx", "exp(fx)*d1**2", "x*d1**2", "d1", "d1/x", "x*d1",
         "fx*d1", "fx", "fx**3", "sin(fx)", "1/fx**2", "x", "1", "a*d1", "g(fx)", "g(x)*d1", "d1**3", "d2**2",
         "y*d1**2", "log(fx)"]


def generated():
    env = {'fx': fx, 'd1': fx.diff(x), 'd2': fx.diff(x, 2), 'x': x, 'y': y, 'a': a,
           'sin': sin, 'exp': exp, 'log': log, 'g': g}
    for t1, t2, t3 in itertools.combinations(TERMS, 3):
        for sign in ('+', '-'):
            src = '%s + %s %s %s' % (t1, t2, sign, t3)
            try:
                yield src, eval(src, env)
            except Exception:
                continue


def main():
    import symox as M
    ns = {n: getattr(M, n) for n in dir(M) if not n.startswith('_')}
    for c in json.load(open(sys.argv[1])):
        try:
            eq = eval(c['eq'], ns)
            func = eval(c['func'], ns)
        except Exception:
            continue
        check(c['name'], eq, func)
    for src, eq in generated():
        check(src, eq, fx)
    tot = collections.Counter()
    for k, v in stats.items():
        tot[k.split(' ')[0]] += v
    print(dict(tot))
    print(sorted(stats.items()))


main()
