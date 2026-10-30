"""limgaps: exp factors merged in Gruntz, loggamma at -oo, dominant terms of
Order sums, AccumBounds series values of floor, ceiling and frac, bounded
powers of sin and cos, moved Order points, against SymPy 1.14."""
import pytest

from symox import *  # noqa: F401,F403
from symox import _core
from symox.calculus.accumulationbounds import AccumBounds

del test  # noqa: F821

x, y = symbols('x y')
p = Symbol('p', positive=True)
r = Symbol('r', real=True)


def test_gruntz_merges_cancelling_exp_factors():
    assert limit(((x**(x + 1) + (x + 1)**x)/x**(x + 1))**x, x, oo) == exp(E)
    assert limit(x**(-x - 1)*(x**(x + 1) + (x + 1)**x), x, oo) == 1


@pytest.mark.parametrize('expr, expected', [
    (loggamma(x), zoo),
    (x*loggamma(x), zoo),
    (2*loggamma(x), zoo),
    (loggamma(x)**2, zoo),
    (1/loggamma(x), 0),
    (loggamma(2*x), zoo),
    (loggamma(-x), oo),
    (loggamma(x**2), oo),
    (loggamma(-x**2), zoo),
    (loggamma(x) + 3, zoo),
])
def test_loggamma_at_minus_infinity(expr, expected):
    assert limit(expr, x, -oo) == expected


@pytest.mark.parametrize('expr', [
    loggamma(x)/x, loggamma(x) + x, loggamma(x) + loggamma(x + 1), exp(loggamma(x)),
    loggamma(x) - x*log(-x),
])
def test_loggamma_at_minus_infinity_held(expr):
    # SymPy 1.14 returns these unevaluated as well
    assert isinstance(limit(expr, x, -oo), Limit)


def test_factorial_over_x_at_minus_infinity_stays_held():
    # SymPy 1.14 gives 0 (through 0*gamma(-oo)); the quotient oscillates
    # through poles, the held value is kept pending user validation.
    assert isinstance(limit(factorial(x)/x, x, -oo), Limit)
    assert isinstance(limit(gamma(x)/x, x, -oo), Limit)


@pytest.mark.parametrize('expr, expected', [
    (lambda: O(x + exp(1/x)), lambda: O(exp(1/x))),
    (lambda: O(exp(1/x) + 1), lambda: O(exp(1/x))),
    (lambda: O(log(x) + log(log(x))), lambda: O(log(x))),
    (lambda: O(log(log(x))), lambda: O(log(log(x)))),
    (lambda: O(exp(-x) + x**-4, (x, -oo)), lambda: O(exp(-x), (x, -oo))),
    (lambda: O(2**x + factorial(x), (x, oo)), lambda: O(factorial(x), (x, oo))),
    (lambda: O(2**x + 3**x, (x, oo)), lambda: O(exp(x*log(3)), (x, oo))),
])
def test_order_dominant_terms(expr, expected):
    assert expr() == expected()


ORDER_ONE = "Order(Integer(1), Tuple(Symbol('x'), Integer(0)))"


def test_limit_of_order_elsewhere():
    assert srepr(limit(1 + O(1/x), x, 1)) == ORDER_ONE


def test_order_point_errors_and_imaginary_infinity():
    with pytest.raises(NotImplementedError):
        Order(x*y, (x, 0), (y, 1))
    with pytest.raises(NotImplementedError):
        O(x*y).subs(x, cos(x))
    with pytest.raises(NotImplementedError):
        O(x*y).subs(x, x + 1)
    assert str(O(x**2, (x, oo)).subs(x, cos(x))) == 'O(1, (x, oo*I))'
    assert str(O(x**2, (x, oo)).subs(x, sin(x))) == 'O(1, (x, -oo*I))'
    assert str(O(x**2, (x, -oo)).subs(x, cos(x))) == 'O(1, (x, -oo*I))'
    assert str(O(x**2, (x, oo)).subs(x, 2*cos(x) + 1)) == 'O(1, (x, oo*I))'
    assert O(x**2, (x, oo)).subs(x, tan(x)) == O((x - pi/2)**(-2), (x, pi/2))
    assert srepr(O(x**2, (x, 1)).subs(x, cos(x))) == ORDER_ONE


@pytest.mark.parametrize('expr, expected', [
    (lambda: frac(x).series(x, oo), lambda: AccumBounds(0, 1) + O(x**-6, (x, oo))),
    (lambda: frac(x).series(x, oo, 3), lambda: AccumBounds(0, 1) + O(x**-3, (x, oo))),
    (lambda: (2*frac(1/x)).series(x), lambda: AccumBounds(0, 2) + O(x**6)),
    (lambda: (frac(1/x) + x).series(x), lambda: AccumBounds(0, 1) + x + O(x**6)),
    (lambda: floor(1/p).series(p), lambda: AccumBounds(-1 + 1/p, 1/p)),
    (lambda: ceiling(1/p).series(p), lambda: AccumBounds(1/p, 1 + 1/p)),
    (lambda: floor(p).series(p, oo), lambda: AccumBounds(p - 1, p)),
    (lambda: floor(1/p).series(p, 0, 0), lambda: 1/p + O(1, p)),
    (lambda: floor(1/x).nseries(x, 0, 3), lambda: AccumBounds(-1, 0) + 1/x),
    (lambda: ceiling(1/x).nseries(x, 0, 3), lambda: AccumBounds(0, 1) + 1/x),
])
def test_accum_bounds_series(expr, expected):
    assert expr() == expected()


def test_accum_bounds_series_without_terms():
    assert str(frac(1/x).series(x, 0, 0)) == 'O(1)'


def test_accum_bounds_series_errors_kept():
    with pytest.raises(ValueError):
        floor(1/r).series(r)
    with pytest.raises(ValueError):
        floor(1/x).series(x)


@pytest.mark.parametrize('expr, expected', [
    (sin(x)**2, AccumBounds(0, 1)),
    (cos(x)**3 + 1, AccumBounds(0, 2)),
    (sin(x)**2 - sin(x), AccumBounds(-1, 2)),
    (sin(x)**2 + cos(x)**2, AccumBounds(0, 2)),
    (2*sin(x)**4, AccumBounds(0, 2)),
    (sin(x)**3/x + cos(x)**2, AccumBounds(0, 1)),
    (sin(2*x)**2 + cos(x), AccumBounds(-1, 2)),
])
def test_bounded_powers(expr, expected):
    assert limit(expr, x, oo) == expected
    assert _core._limit_bounded_affine(expr, x, oo, '+') is not None


@pytest.mark.parametrize('expr, expected', [
    (cos(x)*floor(x)/x, AccumBounds(-1, 1)),
    (cos(x)*floor(x)/x**2, 0),
    (sin(x)**2*floor(x)/x, AccumBounds(0, 1)),
    (-3*cos(x)*ceiling(x)/x, AccumBounds(-3, 3)),
    (cos(x)*(floor(x) + 1)/(x + 2), AccumBounds(-1, 1)),
    (cos(x)*floor(x)/(x*log(x)), 0),
])
def test_bounded_factor_times_rounding(expr, expected):
    assert limit(expr, x, oo) == expected


def test_bounded_sum_with_rounding_stays_held():
    # SymPy 1.14 leaves this unevaluated (term-wise heuristics meet oo - oo)
    assert isinstance(limit(sin(x) + floor(x) - x, x, oo), Limit)
