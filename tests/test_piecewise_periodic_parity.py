"""Periodic inequalities, Piecewise simplification and Piecewise powers.

The expected strings are the SymPy 1.14 outputs, except where the reference
answer is wrong; those cases assert the correct value and say why.
"""
import pytest

from symox import (S, Symbol, symbols, sin, cos, tan, pi, Abs, Interval, Eq, Ne, Le, Lt,
                   Ge, Gt, And, Or, Not, Xor, Piecewise, ITE, KroneckerDelta, solveset,
                   solve_univariate_inequality, reduce_inequalities, solve, integrate,
                   piecewise_exclusive, sqrt, Function)
from symox.functions.elementary.piecewise import piecewise_simplify
from symox.logic.boolalg import simplify_univariate

x, y, z = symbols('x y z')
r = Symbol('r', real=True)


@pytest.mark.parametrize('expr, expected', [
    (lambda: solveset(sin(x) > 0, x, S.Reals), 'Interval.open(0, pi)'),
    (lambda: solveset(cos(x) < 0, x, Interval(0, 2*pi)), 'Interval.open(pi/2, 3*pi/2)'),
    (lambda: solveset(sin(x) > S.Half, x, Interval(0, 2*pi)), 'Interval.open(pi/6, 5*pi/6)'),
    (lambda: solveset(tan(x) > 0, x, Interval(0, pi)), 'Interval.open(0, pi/2)'),
    (lambda: solveset(sin(x) > 0, x, Interval(-pi, pi)), 'Interval.open(0, pi)'),
    (lambda: solveset(sin(2*x) < 0, x, Interval(0, pi)), 'Interval.open(pi/2, pi)'),
    (lambda: solveset(sin(x)**2 > S(1)/4, x, Interval(0, pi)), 'Interval.open(pi/6, 5*pi/6)'),
    # SymPy 1.14 leaves a ConditionSet: its function_range fails on the
    # Piecewise derivative of Abs(sin(x)) (as_set of a periodic relational)
    (lambda: solveset(Abs(sin(x)) < S.Half, x, Interval(0, 2*pi)),
     'Union(Interval.Ropen(0, pi/6), Interval.open(5*pi/6, 7*pi/6), Interval.Lopen(11*pi/6, 2*pi))'),
    (lambda: solveset(cos(x) > 0, x, S.Reals),
     'Union(Interval.Ropen(0, pi/2), Interval.open(3*pi/2, 2*pi))'),
    (lambda: solveset(sin(x) > 2, x, S.Reals), 'EmptySet'),
    (lambda: solveset(sin(x) <= 1, x, S.Reals), 'Reals'),
    (lambda: solveset(sin(x) > 0, x, Interval(0, 2*pi)).as_relational(x), '(0 < x) & (x < pi)'),
    (lambda: reduce_inequalities(sin(x) > 0, x), '(0 < x) & (x < pi)'),
    (lambda: reduce_inequalities(cos(x) >= 0, x),
     '((0 <= x) & (x <= pi/2)) | ((3*pi/2 <= x) & (x < 2*pi))'),
    (lambda: reduce_inequalities(cos(x) + 2 > 0, x), 'True'),
    (lambda: solve_univariate_inequality(sin(x) > 0, x, domain=Interval(0, 2*pi)),
     '(0 < x) & (x < pi)'),
    (lambda: solve_univariate_inequality(sin(x) > 0, x, relational=False), 'Interval.open(0, pi)'),
    (lambda: solve(sin(x) > 0, x), '(0 < x) & (x < pi)'),
    (lambda: (sin(x) + 2 > 0).as_set(), 'Interval(-oo, oo)'),
    (lambda: (sin(x) + 2 < 0).as_set(), 'EmptySet'),
    (lambda: ((sin(x) + 2 > 0) & (x > 1)).as_set(), 'Interval.open(1, oo)'),
    (lambda: integrate(Piecewise((1, sin(x) > 0), (0, True)), (x, 0, 2*pi)), 'pi'),
    (lambda: integrate(Piecewise((1, sin(x) > 0), (0, True)), (x, 0, pi)), 'pi'),
    (lambda: integrate(Piecewise((sin(x), sin(x) > 0), (0, True)), (x, 0, 2*pi)), '2'),
])
def test_periodic_inequalities(expr, expected):
    assert str(expr()) == expected


@pytest.mark.parametrize('expr', [
    lambda: (sin(x) > 0).as_set(),
    lambda: (cos(x) <= 0).as_set(),
    lambda: And(sin(x) > 0, x > 0, x < 2*pi).as_set(),
    lambda: Eq(sin(x), 0).as_set(),
    lambda: Piecewise((1, sin(x) > 0), (2, True)).as_expr_set_pairs(),
])
def test_periodic_as_set_refused(expr):
    with pytest.raises(NotImplementedError) as info:
        expr()
    assert str(info.value) == '\nas_set is not implemented for relationals with periodic solutions'


def test_periodic_inequality_beyond_one_period():
    # SymPy intersects the solutions in [0, 2*pi) with the domain and returns
    # EmptySet on [2*pi, 3*pi]; the inequality is solved on the domain itself
    assert solveset(sin(x) > 0, x, Interval(2*pi, 3*pi)) == Interval.open(2*pi, 3*pi)
    # SymPy misses the endpoint 2*pi, where sin vanishes
    assert solveset(sin(x) >= 0, x, Interval(0, 2*pi)).contains(2*pi) == S.true
    # SymPy integrates one period only and returns pi
    assert integrate(Piecewise((1, sin(x) > 0), (0, True)), (x, 0, 4*pi)) == 2*pi
    assert integrate(Piecewise((1, sin(x) > 0), (0, True)), (x, -pi, 0)) == 0


def test_periodic_condition_kept_on_the_whole_line():
    # SymPy rewrites sin(x) > 0 as 0 < x < pi, which is false at x = 2*pi + 1
    pw = Piecewise((1, sin(x) > 0), (0, True))
    assert str(pw.simplify()) == 'Piecewise((1, sin(x) > 0), (0, True))'
    assert str(piecewise_simplify(pw)) == 'Piecewise((1, sin(x) > 0), (0, True))'
    assert pw.simplify().subs(x, 2*pi + 1) == 1


@pytest.mark.parametrize('expr, expected', [
    (lambda: piecewise_simplify(Piecewise((x, x < 1), (x, x < 2), (y, True))),
     'Piecewise((x, (x <= 1) | (x < 2)), (y, True))'),
    (lambda: piecewise_simplify(Piecewise((0, x < 0), (x, Eq(x, 0)), (x**2, True))),
     'Piecewise((0, x < 0), (x**2, True))'),
    (lambda: piecewise_simplify(Piecewise((1, x < 0), (1, x < 1), (2, True))),
     'Piecewise((1, (x <= 0) | (x < 1)), (2, True))'),
    (lambda: piecewise_simplify(Piecewise((x, x > 0), (-x, True))), 'Piecewise((x, x > 0), (-x, True))'),
    (lambda: Piecewise((z, And(Eq(z, 0), Eq(z + y, 0))), (1, True)).simplify(),
     'Piecewise((0, Eq(y, 0) & Eq(z, 0)), (1, True))'),
    (lambda: Piecewise((x + 1, Eq(x, 1)), (2*x, True)).simplify(), '2*x'),
    (lambda: Piecewise((Function('f')(x), Eq(Function('f')(x), 1)), (2, True)).simplify(),
     'Piecewise((1, Eq(f(x), 1)), (2, True))'),
    (lambda: (Piecewise((sin(x)**2 + cos(x)**2, x < 1), (2, True))).simplify(),
     'Piecewise((1, x < 1), (2, True))'),
    (lambda: piecewise_exclusive(Piecewise((1, x < 0), (2, x < 1), (3, True))),
     'Piecewise((1, x < 0), (2, (x >= 0) & (x < 1)), (3, x >= 1))'),
    (lambda: piecewise_exclusive(Piecewise((1, x < 0), (2, x < 1))),
     'Piecewise((1, x < 0), (2, (x >= 0) & (x < 1)), (nan, x >= 1))'),
    (lambda: piecewise_exclusive(Piecewise((1, x < 0), (2, x < 1)), skip_nan=True),
     'Piecewise((1, x < 0), (2, (x >= 0) & (x < 1)))'),
])
def test_piecewise_simplify(expr, expected):
    assert str(expr()) == expected


@pytest.mark.parametrize('expr, expected', [
    (lambda: Le(0, x).simplify(), 'x >= 0'),
    (lambda: Lt(1, x).simplify(), 'x > 1'),
    (lambda: Lt(y, x).simplify(), 'x > y'),
    (lambda: Le(2*x, 4).simplify(), 'x <= 2'),
    (lambda: Ge(0, x - 1).simplify(), 'x <= 1'),
    (lambda: Gt(-x, 2).simplify(), 'x < -2'),
    (lambda: Lt(x**2 + 2*x, 3).simplify(), 'x**2 + 2*x < 3'),
    (lambda: Lt(2*x + 4*y, 6).simplify(), 'x + 2*y < 3'),
    (lambda: Lt(2*x + 4*y, 6*x).simplify(), 'x > y'),
    (lambda: Lt(pi, 3).simplify(), 'False'),
    (lambda: Le(sin(x)**2 + cos(x)**2, 1).simplify(), 'True'),
    (lambda: Ge(r**2, 4*r).simplify(), 'r**2 - 4*r >= 0'),
    (lambda: Not(Or(x < 0, x < 1)).simplify(), 'x >= 1'),
    (lambda: And(0 <= x, y < 1).simplify(), '(x >= 0) & (y < 1)'),
    (lambda: (And(x > 0, x < 1) | (y > 1)).simplify(), '(y > 1) | ((x > 0) & (x < 1))'),
    (lambda: Xor(x > 0, y > 0).simplify(), '(x > 0) ^ (y > 0)'),
])
def test_relational_and_boolean_simplify(expr, expected):
    assert str(expr()) == expected


def test_simplify_univariate_endpoints():
    # SymPy tests the endpoints of later intervals on an earlier piece; the
    # whole condition is used here, so x = 0 stays in the result
    rv = simplify_univariate(Or(x <= 0, x >= 1))
    assert rv.subs(x, 0) == S.true and rv.subs(x, 1) == S.true
    assert rv.subs(x, S.Half) == S.false


@pytest.mark.parametrize('expr, expected', [
    (lambda: Piecewise((1, x < 0), (2, True))**2, 'Piecewise((1, x < 0), (4, True))'),
    (lambda: Piecewise((1, x < 0), (2, True))**y, 'Piecewise((1, x < 0), (2**y, True))'),
    (lambda: sqrt(Piecewise((x, x < 1), (4, True))), 'Piecewise((sqrt(x), x < 1), (2, True))'),
    (lambda: Piecewise((x, x < 1), (4, True))**-1, 'Piecewise((1/x, x < 1), (1/4, True))'),
])
def test_piecewise_power(expr, expected):
    assert str(expr()) == expected


@pytest.mark.parametrize('expr, expected', [
    (lambda: KroneckerDelta(x, y).rewrite(Piecewise), 'Piecewise((0, Ne(x, y)), (1, True))'),
    (lambda: KroneckerDelta(x, 2).rewrite(Piecewise), 'Piecewise((0, Ne(x, 2)), (1, True))'),
    (lambda: Piecewise((y > 0, x < 0), (z > 0, True)).rewrite(ITE), 'ITE(x < 0, y > 0, z > 0)'),
])
def test_piecewise_rewrites(expr, expected):
    assert str(expr()) == expected


@pytest.mark.parametrize('pw, name', [
    (lambda: Piecewise((1, x < 0), (2, True)), 'Integer'),
    (lambda: Piecewise((S.Half, x < 0), (2, True)), 'Half'),
    (lambda: Piecewise((0, x < 0), (2, True)), 'Zero'),
    (lambda: Piecewise((x + 1, x < 0), (2, True)), 'Add'),
    (lambda: Piecewise((y, x < 0), (2, True)), 'Integer'),
])
def test_ite_rewrite_type_error(pw, name):
    with pytest.raises(TypeError) as info:
        pw().rewrite(ITE)
    assert str(info.value) == '\nExpecting Boolean or bool but got `%s`' % name
