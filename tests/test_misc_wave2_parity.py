"""Parity checks against SymPy 1.14 for Min/Max ordering, polarify,
Piecewise series, relationals of non-real sides, as_expr_set_pairs and
piecewise_fold. The expected strings were produced by SymPy 1.14.0."""
import pytest

from symox import (
    Dummy, Eq, Ne, Ge, Gt, I, Integral, Interval, Le, Lt, Max, Min, O,
    Order, Piecewise, Rational, S, Sum, Symbol, Tuple, besselj, cos,
    exp, exp_polar, hyper, latex, log, meijerg, nan, oo, pi, pretty, series,
    sin, sqrt, srepr, symbols, zoo, Abs, Derivative, FiniteSet, Heaviside,
)
from symox.functions.elementary.complexes import polarify
from symox.functions.elementary.piecewise import piecewise_fold

x, y, z = symbols('x y z')
n = Symbol('n')


# -- Min / Max --------------------------------------------------------------

@pytest.mark.parametrize('build, expected', [
    (lambda: Min(x**2, x, sin(x)), 'Min(x, x**2, sin(x))'),
    (lambda: Max(x**3, cos(x), x, exp(x)), 'Max(x, x**3, exp(x), cos(x))'),
    (lambda: Max(y**2, x + y, sin(y), x*y), 'Max(y**2, x*y, x + y, sin(y))'),
    (lambda: Max(x, Min(y**2, y + 1), 3), 'Max(3, x, Min(y**2, y + 1))'),
    (lambda: Max(Min(x**2, x, sin(x)), x + 5), 'Max(x + 5, Min(x, x**2, sin(x)))'),
    (lambda: Min(x + 1, sin(x), evaluate=False), 'Min(x + 1, sin(x))'),
    (lambda: Max(x, cos(x), 2*sin(y)), 'Max(x, 2*sin(y), cos(x))'),
    (lambda: Max(Abs(x), x**2, z), 'Max(x**2, z, Abs(x))'),
    (lambda: Max(y, x), 'Max(x, y)'),
    (lambda: Min(z, x, y), 'Min(x, y, z)'),
])
def test_minmax_str_sorted(build, expected):
    e = build()
    assert str(e) == expected
    assert pretty(e) == expected


def test_minmax_args_and_srepr():
    m = Min(x**2, x, sin(x))
    assert m.args == (x, sin(x), x**2)
    assert srepr(m) == "Min(Symbol('x'), sin(Symbol('x')), Pow(Symbol('x'), Integer(2)))"
    assert Max(x**3, cos(x), x, exp(x)).args == (x, exp(x), cos(x), x**3)
    assert Max(y, x) == Max(x, y) and hash(Max(y, x)) == hash(Max(x, y))


def test_minmax_latex_sorted():
    assert latex(Min(x**2, x, sin(x))) == r'\min\left(x, x^{2}, \sin{\left(x \right)}\right)'
    assert latex(Max(x**3, cos(x), x, exp(x))) == r'\max\left(x, x^{3}, e^{x}, \cos{\left(x \right)}\right)'
    assert latex(Max(x, Min(y**2, y + 1), 3)) == r'\max\left(3, x, \min\left(y^{2}, y + 1\right)\right)'


def test_minmax_derivative():
    assert str(Max(y, x).diff(x)) == 'Heaviside(x - y)'
    assert str(Min(y, x).diff(x)) == 'Heaviside(-x + y)'
    assert str(Max(x, y, z).diff(x)) == 'Heaviside(x - Max(y, z))'
    assert str(Max(x**2, y).diff(x)) == '2*x*Heaviside(x**2 - y)'
    assert Max(y, z).diff(x) == 0


def test_exp_polar_sort_key():
    assert str(exp_polar(I*pi)*exp(x)*sin(y)) == 'exp_polar(I*pi)*exp(x)*sin(y)'
    assert str(exp_polar(2)*sin(x)) == 'exp_polar(2)*sin(x)'
    assert str(exp_polar(I)*exp(x)) == 'exp_polar(I)*exp(x)'


# -- polarify ---------------------------------------------------------------

@pytest.mark.parametrize('build, expected', [
    (lambda: polarify(meijerg([x], [y], [1], [z], x*y), subs=False),
     'meijerg(((x,), (y,)), ((1,), (z,)), x*y)'),
    (lambda: polarify(meijerg([-x], [], [1], [], -x), subs=False),
     'meijerg(((-x,), ()), ((1,), ()), x*exp_polar(I*pi))'),
    (lambda: polarify(meijerg([], [], [0], [], exp(I*pi)*x), subs=False),
     'meijerg(((), ()), ((0,), ()), x*exp_polar(I*pi))'),
    (lambda: polarify(meijerg([1], [2], [3], [4], -y), lift=True),
     'meijerg(((1,), (2,)), ((3,), (4,)), exp_polar(I*pi)*polar_lift(y))'),
    (lambda: polarify(x*meijerg([-1], [2], [], [], -y) + exp(-x), subs=False),
     'x*meijerg(((-1,), (2,)), ((), ()), y*exp_polar(I*pi)) + exp(x*exp_polar(I*pi))'),
    (lambda: polarify(hyper([-1], [2], -x), subs=False), 'hyper((-1,), (2,), -x)'),
    (lambda: polarify(hyper([-x, 2], [3], 2*y), subs=False), 'hyper((2, -x), (3,), 2*y*exp_polar(0))'),
    (lambda: polarify(hyper([Rational(1, 2)], [3], -y), subs=False), 'hyper((1/2,), (3,), -y)'),
    (lambda: polarify(Integral(exp(-x)*y, (x, 0, oo)), subs=False),
     'Integral(y*exp(x*exp_polar(I*pi)), (x, polar_lift(0), oo))'),
    (lambda: polarify(Integral(-x, (x, -1, y)), subs=False, lift=True),
     'Integral(exp_polar(I*pi)*polar_lift(x), (x, exp_polar(I*pi), polar_lift(y)))'),
    (lambda: polarify(Sum(-x**n, (n, 1, 5)), subs=False), 'Sum(x**n*exp_polar(I*pi), (n, 1, 5))'),
    (lambda: polarify(Piecewise((-x, x < -1), (2, True)), subs=False),
     'Piecewise((x*exp_polar(I*pi), x < -1), (2*exp_polar(0), True))'),
    (lambda: polarify(exp(-x)*sin(-2*y), subs=False),
     'exp_polar(I*pi)*exp(x*exp_polar(I*pi))*sin(2*y*exp_polar(0))'),
    (lambda: polarify(-x + 1, lift=True), 'polar_lift(1 - x)'),
    (lambda: polarify(Eq(-x, 2), subs=False), 'Eq(x*exp_polar(I*pi), 2*exp_polar(0))'),
    (lambda: polarify(Tuple(-x, 2), subs=False), '(x*exp_polar(I*pi), 2*exp_polar(0))'),
    (lambda: polarify(besselj(-1, -x), subs=False), 'besselj(1, x)'),
    (lambda: polarify(sqrt(-x), subs=False), '(x*exp_polar(I*pi))**(exp_polar(0)/2)'),
])
def test_polarify(build, expected):
    assert str(build()) == expected


def test_polarify_subs_meijerg():
    e, reps = polarify(meijerg([x], [y], [1], [z], x*y))
    assert str(e) == 'meijerg(((_x,), (_y,)), ((1,), (_z,)), _x*_y)'
    assert e.xreplace(reps) == meijerg([x], [y], [1], [z], x*y)


def test_hyper_orders_then_unpolarifies():
    assert str(hyper([x*exp_polar(I*pi), 2*exp_polar(0)], [3], y)) == 'hyper((2, -x), (3,), y)'


# -- Piecewise series -------------------------------------------------------

P = Piecewise((x**2, x < 0), (exp(x), True))


@pytest.mark.parametrize('build, expected', [
    (lambda: P.series(x, 0, 4), '1 + x + x**2/2 + x**3/6 + O(x**4)'),
    (lambda: series(P, x, 0, 3), '1 + x + x**2/2 + O(x**3)'),
    (lambda: Piecewise((sin(x), x > 0), (cos(x), True)).series(x, 0, 4), 'x - x**3/6 + O(x**4)'),
    (lambda: Piecewise((1/(1 - x), x < S(1)/2), (x, True)).series(x, 0, 3),
     'Piecewise((1 + x + x**2 + O(x**3), x < 1/2), (x, True))'),
    (lambda: Piecewise((x, y > 0), (x**2, True)).series(x, 0, 3), 'Piecewise((x, y > 0), (x**2, True))'),
    (lambda: Piecewise((exp(x), x < 1), (0, True)).series(x, 1, 3), '0'),
    (lambda: Piecewise((exp(x), x > 1), (sin(x), True)).nseries(x, 0, 3),
     'Piecewise((1 + x + x**2/2 + O(x**3), x > 1), (x + O(x**3), True))'),
    (lambda: Piecewise((exp(x), Eq(x, 0)), (sin(x), True)).series(x, 0, 3), 'x + O(x**3)'),
    (lambda: Piecewise((exp(x), x >= 0), (sin(x), True)).series(x, 0, 3), '1 + x + x**2/2 + O(x**3)'),
    (lambda: Piecewise((x, y < 1), (sqrt(x**2), True)).series(x, 0, 3), 'x'),
    (lambda: Piecewise((exp(x), x < 0), (cos(x), True)).series(x, 0, 3, dir='-'), '1 + x + x**2/2 + O(x**3)'),
    (lambda: Piecewise((exp(x), x < 2), (cos(x), True)).series(x, 2, 3),
     'cos(2) - (x - 2)*sin(2) - (x - 2)**2*cos(2)/2 + O((x - 2)**3, (x, 2))'),
    (lambda: Piecewise((exp(x), x < 2), (cos(x), True)).series(x, 2, 3, dir='-'),
     'exp(2) - (2 - x)*exp(2) + (2 - x)**2*exp(2)/2 + O((x - 2)**3, (x, 2))'),
    (lambda: Piecewise((1/x, x > 0), (0, True)).series(x, 0, 3), '1/x'),
    (lambda: Piecewise((exp(y*x), y > 0), (x, True)).series(x, 0, 3),
     'Piecewise((1 + x*y + x**2*y**2/2 + O(x**3), y > 0), (x, True))'),
    (lambda: Piecewise((log(1 + x), x < 1), (x, True)).series(x, 0, 4),
     'Piecewise((x - x**2/2 + x**3/3 + O(x**4), x < 1), (x, True))'),
    (lambda: Piecewise((x**2, x < 0), (exp(x), True)).nseries(x, 0, 3),
     'Piecewise((x**2, x < 0), (1 + x + x**2/2 + O(x**3), True))'),
    (lambda: Piecewise((sin(x)/x, x < 0), (1, True)).series(x), '1'),
    (lambda: series(Piecewise((1/(1 + x), x > 1), (x + x**3, True))),
     'Piecewise((1 - x + x**2 - x**3 + x**4 - x**5 + O(x**6), x > 1), (x**3 + x, True))'),
    (lambda: Piecewise((y, x < 0), (2, True)).series(x, 0, 3), '2'),
])
def test_piecewise_series(build, expected):
    assert str(build()) == expected


def test_piecewise_series_multivariate_error():
    with pytest.raises(ValueError, match='x must be given for multivariate functions.'):
        series(Piecewise((x, x > y), (1, True)))


def test_order_subs_moves_point():
    assert str(O(x**3).subs(x, x - 2)) == 'O((x - 2)**3, (x, 2))'
    assert str(O(x**3).subs(x, 2 - x)) == 'O((x - 2)**3, (x, 2))'
    assert str(O(x**3).subs(x, -x)) == 'O(x**3)'
    assert str((1 + x + O(x**3)).subs(x, x - 2)) == '-1 + x + O((x - 2)**3, (x, 2))'
    assert str(exp(-x).series(x, 0, 3).subs(x, -x)) == '1 + x + x**2/2 + O(x**3)'


# -- relationals ------------------------------------------------------------

r = Symbol('r', real=True)
c = Symbol('c', complex=True)
w = Symbol('w', real=False)


@pytest.mark.parametrize('build, expected', [
    (lambda: y < y, 'y < y'), (lambda: y <= y, 'y <= y'),
    (lambda: y > y, 'y > y'), (lambda: y >= y, 'y >= y'),
    (lambda: Lt(y, y), 'y < y'), (lambda: Le(y, y), 'y <= y'),
    (lambda: Gt(y, y), 'y > y'), (lambda: Ge(y, y), 'y >= y'),
    (lambda: r < r, 'False'), (lambda: r <= r, 'True'), (lambda: Lt(r, r), 'False'),
    (lambda: c < c, 'c < c'), (lambda: c <= c, 'c <= c'), (lambda: Lt(c, c), 'c < c'),
    (lambda: x < x + 1, 'x < x + 1'), (lambda: c < c + 1, 'c < c + 1'),
    (lambda: r < r + 1, 'True'),
    (lambda: w < w, 'w < w'), (lambda: Lt(x, x + I), 'x < x + I'),
    (lambda: Lt(x*I, x*I), 'I*x < I*x'),
    (lambda: Eq(y, y), 'True'), (lambda: Ne(y, y), 'False'),
    (lambda: oo < oo, 'False'),
])
def test_relational_self_comparison(build, expected):
    assert str(build()) == expected


@pytest.mark.parametrize('build, message', [
    (lambda: I < I, 'Invalid comparison of non-real I'),
    (lambda: Lt(I, I), 'Invalid comparison of non-real I'),
    (lambda: Le(I, I), 'Invalid comparison of non-real I'),
    (lambda: (1 + I) < (1 + I), r'Invalid comparison of non-real 1 \+ I'),
    (lambda: Symbol('i', imaginary=True) < Symbol('i', imaginary=True), 'Invalid comparison of non-real i'),
    (lambda: zoo < zoo, 'Invalid comparison of non-real zoo'),
    (lambda: Le(zoo, zoo), 'Invalid comparison of non-real zoo'),
    (lambda: nan < nan, 'Invalid NaN comparison'),
])
def test_relational_non_real_type_error(build, message):
    with pytest.raises(TypeError, match=message):
        build()


# -- as_expr_set_pairs ------------------------------------------------------

@pytest.mark.parametrize('build, expected', [
    (lambda: Piecewise((x, x < 1), (x**2, x > 3), (x + 1, True)).as_expr_set_pairs(),
     '[(x, Interval.open(-oo, 1)), (x**2, Interval.open(3, oo)), (x + 1, Interval(1, 3))]'),
    (lambda: Piecewise((x, x < 1), (x**2, x > 3), (x + 1, True)).as_expr_set_pairs(Interval(0, 5)),
     '[(x, Interval.Ropen(0, 1)), (x**2, Interval.Lopen(3, 5)), (x + 1, Interval(1, 3))]'),
    (lambda: Piecewise((1, x < 0), (2, True)).as_expr_set_pairs(),
     '[(1, Interval.open(-oo, 0)), (2, Interval(0, oo))]'),
    (lambda: Piecewise((1, (x > 0) & (x < 2)), (2, x > 1), (3, True)).as_expr_set_pairs(),
     '[(1, Interval.open(0, 2)), (2, Interval(2, oo)), (3, Interval(-oo, 0))]'),
    (lambda: Piecewise((1, Eq(x, 0)), (2, True)).as_expr_set_pairs(),
     '[(1, {0}), (2, Union(Interval.open(-oo, 0), Interval.open(0, oo)))]'),
    (lambda: Piecewise((1, x < 0), (2, x < 3)).as_expr_set_pairs(),
     '[(1, Interval.open(-oo, 0)), (2, Interval.Ropen(0, 3))]'),
    (lambda: Piecewise((1, Ne(x, 1)), (2, True)).as_expr_set_pairs(Interval(0, 3)),
     '[(1, Union(Interval.Ropen(0, 1), Interval.Lopen(1, 3))), (2, {1})]'),
    (lambda: Piecewise((1, x > 0)).as_expr_set_pairs(S.Reals), '[(1, Interval.open(0, oo))]'),
    (lambda: Piecewise((1, (x < 0) | (x > 5)), (3, True)).as_expr_set_pairs(),
     '[(1, Union(Interval.open(-oo, 0), Interval.open(5, oo))), (3, Interval(0, 5))]'),
    (lambda: Piecewise((1, Eq(x, I)), (2, True)).as_expr_set_pairs(S.Complexes), '[(2, Complexes)]'),
    (lambda: Piecewise((x, x < 1), (x**2, True)).as_expr_set_pairs(Interval(2, 3)), '[(x**2, Interval(2, 3))]'),
    (lambda: Piecewise((1, x < 0), (2, True)).as_expr_set_pairs(FiniteSet(-1, 1, 2)), '[(1, {-1}), (2, {1, 2})]'),
    (lambda: Piecewise((1, Eq(y, 0)), (2, True)).as_expr_set_pairs(),
     '[(1, {0}), (2, Union(Interval.open(-oo, 0), Interval.open(0, oo)))]'),
])
def test_as_expr_set_pairs(build, expected):
    assert str(build()) == expected


def test_as_expr_set_pairs_errors():
    with pytest.raises(NotImplementedError) as info:
        Piecewise((y, x < y), (2, True)).as_expr_set_pairs()
    assert str(info.value) == '\nmultivariate conditions are not handled.'
    with pytest.raises(ValueError) as info:
        Piecewise((1, x < 0), (2, True)).as_expr_set_pairs(S.Complexes)
    assert str(info.value) == ('\nInequalities in the complex domain are not supported. Try the real\n'
                               'domain by setting domain=S.Reals')


# -- piecewise_fold ---------------------------------------------------------

p1 = Piecewise((x, x < 1), (1, True))
p2 = Piecewise((x**2, x < 2), (0, True))
p3 = Piecewise((y, y > 0), (-y, True))
q2 = Piecewise((x**2, x < 1), (0, True))
q4 = Piecewise((x**3, x < 1), (5, x < 2), (0, True))


@pytest.mark.parametrize('build, expected', [
    (lambda: piecewise_fold(p1 + p2),
     'Piecewise((x**2 + x, (x < 1) & (x < 2)), (x, x < 1), (x**2 + 1, x < 2), (1, True))'),
    (lambda: piecewise_fold(p1 + p3),
     'Piecewise((x + y, (y > 0) & (x < 1)), (y + 1, y > 0), (x - y, x < 1), (1 - y, True))'),
    (lambda: piecewise_fold(Piecewise((1, Eq(Piecewise((x, x < 0), (1, True)), 0)), (2, True))),
     'Piecewise((1, Eq(x, 0) & (x < 0)), (2, True))'),
    (lambda: piecewise_fold(p1 + q2 + q4), 'Piecewise((x**3 + x**2 + x, x < 1), (6, x < 2), (1, True))'),
    (lambda: piecewise_fold(p1*q2*q4), 'Piecewise((x**6, x < 1), (0, True))'),
    (lambda: piecewise_fold(p1 + q2, evaluate=None), 'Piecewise((x**2 + x, x < 1), (1, True))'),
    (lambda: piecewise_fold(Piecewise((1, x < 1), (1, x < 2), (3, True))*y, evaluate=None),
     'Piecewise((y, (x < 1) | (x < 2)), (3*y, True))'),
    (lambda: piecewise_fold(Piecewise((p1, x < 0), (q2, True)), evaluate=None),
     'Piecewise((x, (x < 0) & (x < 1)), (1, x < 0), (x**2, x < 1), (0, True))'),
    (lambda: piecewise_fold(Piecewise((x, x < 1), (y, True)) + Piecewise((x, x < 1), (z, True))
                            + Piecewise((2, y < 0), (3, True))),
     'Piecewise((2*x + 2, (x < 1) & (y < 0)), (2*x + 3, x < 1), (y + z + 2, y < 0), (y + z + 3, True))'),
    (lambda: piecewise_fold(Piecewise((1, Piecewise((x < 0, y > 0), (x > 1, True))), (2, True))),
     'Piecewise((1, ((y <= 0) | (x < 0)) & ((x > 1) | (y > 0))), (2, True))'),
    (lambda: piecewise_fold(Piecewise((x, x > 0), (0, True)) + Piecewise((1, x > 0), (2, True)), evaluate=False),
     'Piecewise((x + 1, x > 0), (2, True))'),
    (lambda: piecewise_fold(log(Piecewise((x, x > 0), (1, True))) + 1), 'Piecewise((log(x) + 1, x > 0), (1, True))'),
    (lambda: piecewise_fold(Piecewise((x, x < 1)) + Piecewise((x**2, x < 1))), 'Piecewise((x**2 + x, x < 1))'),
    (lambda: piecewise_fold(Piecewise((x, x < 1), (2, x < 3)) + Piecewise((x**2, x < 1), (4, True))),
     'Piecewise((x**2 + x, x < 1), (6, x < 3))'),
    (lambda: piecewise_fold(Piecewise((Piecewise((x, y < 0), (1, True)) + Piecewise((2, y < 0), (3, True)), x < 0),
                                      (0, True))),
     'Piecewise((x + 2, (x < 0) & (y < 0)), (4, x < 0), (0, True))'),
    (lambda: piecewise_fold(Integral(p1, (x, 0, 3))), 'Piecewise((Integral(x, (x, 0, 3)), x < 1), (Integral(1, (x, 0, 3)), True))'),
    (lambda: piecewise_fold(Sum(p3, (y, 1, 3))), 'Piecewise((Sum(y, (y, 1, 3)), y > 0), (Sum(-y, (y, 1, 3)), True))'),
    (lambda: piecewise_fold(Derivative(p1, x)), 'Piecewise((Derivative(x, x), x < 1), (Derivative(1, x), True))'),
    (lambda: piecewise_fold(x), 'x'),
])
def test_piecewise_fold(build, expected):
    assert str(build()) == expected
