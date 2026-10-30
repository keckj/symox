"""Differential test: native nth_order_reducible vs the reference matcher.

usage: difftest_reducible.py cases.json
Prints one line per disagreement and a summary.
"""
import collections
import itertools
import json
import signal
import sys

from symox import *  # noqa
from symox.solvers.ode import _classify_native as cn
from symox.solvers.ode.single import NthOrderReducible, SingleODEProblem

globals().pop("test", None)

stats = collections.Counter()
x, y, a = symbols('x y a')
f, g = Function('f'), Function('g')
fx = f(x)
SAVED = cn.nth_order_reducible


class TO(Exception):
    pass


def alarm(*args):
    raise TO()


signal.signal(signal.SIGALRM, alarm)


def decide(eq, func, x):
    solver = NthOrderReducible(SingleODEProblem(eq, func, x))
    ok = solver._matches()
    return ok, (solver.smallest if ok else None)


def check(name, eq, func):
    if isinstance(eq, Equality):
        eq = eq.lhs - eq.rhs
    if len(getattr(func, 'args', ())) != 1:
        return
    x = func.args[0]
    signal.alarm(60)
    try:
        prob = SingleODEProblem(eq, func, x)
        nat = cn.nth_order_reducible(NthOrderReducible(prob))
        if nat is None:
            stats['decline'] += 1
            return
        cn.nth_order_reducible = lambda solver: None
        try:
            ref = decide(eq, func, x)
        finally:
            cn.nth_order_reducible = SAVED
        got = (bool(nat[0]), nat[1])
        if got == (bool(ref[0]), ref[1]):
            stats['match' if got[0] else 'agree'] += 1
        else:
            stats['DISAGREE'] += 1
            print('DISAGREE', name, 'ref', ref, 'native', nat, flush=True)
    except TO:
        stats['timeout'] += 1
    except Exception as exc:
        stats['error'] += 1
        print('ERROR', name, type(exc).__name__, str(exc)[:100], flush=True)
    finally:
        signal.alarm(0)


TERMS = ["0", "fx", "d1", "d2", "d3", "d1**2", "d2**2", "x*d1", "x*d2", "sin(d1)", "exp(x)*d2", "fx*d2",
         "d1*d3", "g(x)*d2", "x", "1", "d2/d1", "a*d3", "Derivative(f(y), y)", "f(2*x)", "d1.diff(y)",
         "Integral(d1, x)", "1.5*d2", "f(x)**2"]


def generated():
    env = {'fx': fx, 'd1': fx.diff(x), 'd2': fx.diff(x, 2), 'd3': fx.diff(x, 3), 'x': x, 'y': y, 'a': a,
           'sin': sin, 'exp': exp, 'g': g, 'f': f, 'Derivative': Derivative, 'Integral': Integral}
    for t1, t2, t3 in itertools.combinations(TERMS, 3):
        src = '%s + %s - (%s)' % (t1, t2, t3)
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
    print(dict(stats))


main()
