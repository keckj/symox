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


@pytest.mark.parametrize('expr, expected', [
    (lambda: O(exp(1/x) + x**-20), lambda: O(exp(1/x))),
    (lambda: O(exp(1/x)) + O(x), lambda: O(exp(1/x))),
    (lambda: O(exp(2/x)) + O(x**-3), lambda: O(exp(2/x))),
    (lambda: O(exp(1/x)) + O(exp(2/x)), lambda: O(exp(2/x))),
    (lambda: O(exp(-1/x)) + O(x), lambda: O(x)),
    (lambda: O(exp(1/x)) + 1, lambda: O(exp(1/x))),
    (lambda: O(x) + exp(-1/x), lambda: O(x)),
    (lambda: O(x + exp(-1/x)), lambda: O(x)),
    (lambda: O(x**2*exp(1/x) + x), lambda: O(x**2*exp(1/x))),
    (lambda: O(exp(1/x)*(x + 1)), lambda: O(exp(1/x))),
    (lambda: O(x + factorial(x), (x, oo)), lambda: O(factorial(x), (x, oo))),
    (lambda: O(x*factorial(x) + x**2, (x, oo)), lambda: O(x*factorial(x), (x, oo))),
    (lambda: O(x**2 + x**x, (x, oo)), lambda: O(exp(x*log(x)), (x, oo))),
    (lambda: O(x*2**x, (x, oo)), lambda: O(x*exp(x*log(2)), (x, oo))),
])
def test_order_sums_and_containment(expr, expected):
    assert expr() == expected()


def test_order_containment_undecided_keeps_terms():
    # limit(exp(1/x)/exp(2/x), x, 0+) is 0: not contained, both kept
    assert (O(exp(1/x)) + exp(2/x)).has(exp(2/x))


@pytest.mark.parametrize('expr, expected', [
    (lambda: O(x*sin(x), (x, oo)), lambda: O(x, (x, oo))),
    (lambda: O(sin(x), (x, oo)), lambda: O(1, (x, oo))),
    (lambda: O(sin(x)**2, (x, oo)), lambda: O(1, (x, oo))),
    (lambda: O(cos(x)/x, (x, oo)), lambda: O(1/x, (x, oo))),
    (lambda: O(tan(x), (x, oo)), lambda: O(1, (x, oo))),
    (lambda: O(sin(1/x)), lambda: O(1, x)),
    (lambda: O(log(x)*sin(1/x)), lambda: O(log(x))),
    (lambda: O(x*sin(x) + 1, (x, oo)), lambda: O(x, (x, oo))),
    (lambda: O(x**x*sin(x), (x, oo)), lambda: O(exp(x*log(x)), (x, oo))),
    (lambda: O(sin(x) + cos(x), (x, oo)), lambda: O(1, (x, oo))),
])
def test_order_bounded_trig_factors(expr, expected):
    assert expr() == expected()


def test_order_trig_factors_kept_like_reference():
    # a product of functions only raising PoleError, and an argument whose
    # value at 0 is nan, are kept by the reference too
    assert O(sin(x)*exp(x), (x, oo)).expr == exp(x)*sin(x)
    assert O(sin(1/x - 1/x**2)).expr == sin(1/x - 1/x**2)


@pytest.mark.parametrize('expr', [
    lambda: frac(x**2 + x).series(x, oo),
    lambda: frac(exp(1/x)).series(x),
    lambda: (frac(exp(1/x)) + x).series(x),
    lambda: (2*frac(exp(1/x))).series(x),
    lambda: frac(sin(1/x)).series(x),
])
def test_frac_series_with_nan_argument(expr):
    assert expr() is S.NaN


def test_floor_series_laurent_argument():
    assert floor(1/p + p).series(p) == AccumBounds(p - 1 + 1/p, p + 1/p)
    assert ceiling(1/p + p).series(p) == AccumBounds(p + 1/p, p + 1 + 1/p)


def test_bounded_limit_absorbs_symbolic_constant():
    assert limit(atan(x) + sin(x), x, oo) == AccumBounds(-1 + pi/2, 1 + pi/2)
    assert limit(sin(x) + cos(x) + E, x, oo) == AccumBounds(-2 + E, 2 + E)


_e3 = lambda: exp(x).series(x, 0, 3)  # noqa: E731


@pytest.mark.parametrize('call,expected', [
    (lambda: _e3()*O(x), 'O(x)'),
    (lambda: _e3()*O(x**2), 'O(x**2)'),
    (lambda: (_e3()*_e3()).expand(), '1 + 2*x + 2*x**2 + O(x**3)'),
    (lambda: sin(x).series(x, 0, 4)*O(x**2), 'O(x**3)'),
    (lambda: (1/x + _e3())*O(x), 'O(1)'),
    (lambda: O(x)*(x + O(x**2)), 'O(x**2)'),
    (lambda: O(x*(1 + x + O(x**2))), 'O(x)'),
    (lambda: O(x, (x, oo))*exp(1/x).series(x, oo, 3), 'O(x, (x, oo))'),
    (lambda: O(1/x, (x, oo))*(x + O(1, (x, oo))), 'O(1, (x, oo))'),
    (lambda: O(x*y)*(1 + x + O(x**2)), 'O(x**3 + x*y, x, y)'),
    (lambda: O(x*y)*(1 + x + O(x**2, x, y)), 'O(x*y, x, y)'),
    (lambda: O(x**2, x, y)*(1 + x + y + O(x**2, x, y)), 'O(x**2, x, y)'),
    (lambda: cos(x).series(x, 0, 4)*O(x**3) + x, 'x + O(x**3)'),
    (lambda: O(x)*(x + 1/x + O(x**2)), 'O(1)'),
    (lambda: O(x)*log(x)*(1 + x + O(x**2)), 'O(x*log(x))'),
    (lambda: O(x)*(2 + O(x)), 'O(x)'),
])
def test_order_absorbs_the_terms_of_products(call, expected):
    assert str(call()) == expected


def test_limit_in_a_non_symbol_variable_follows_limit_doit():
    assert limit(x, x**2, 0) == x
    assert limit(x**2, x**2, 3) == 3
    assert limit(y, x**2, 0) == y
    with pytest.raises(TypeError, match=r'^x\*\*2 should be of symbol type$'):
        limit(x**2 + 1, x**2, 0)


@pytest.mark.parametrize('direction,error,message', [
    ('-+', ValueError, "^direction must be one of '\\+', '-' or '\\+-', not -\\+$"),
    ('foo', ValueError, "^direction must be one of '\\+', '-' or '\\+-', not foo$"),
    (1, TypeError, "^direction must be of type basestring or Symbol, not <class 'int'>$"),
])
def test_limit_direction_errors_match_the_reference(direction, error, message):
    with pytest.raises(error, match=message):
        limit(1/x, x, 0, direction)
    with pytest.raises(error, match=message):
        Limit(1/x, x, 0, direction)


def test_public_gruntz_returns_accumulation_bounds_of_oscillation():
    from symox.series.gruntz import gruntz
    assert gruntz(sin(x), x, oo) == AccumBounds(-1, 1)
    assert gruntz(cos(x) + 1, x, oo) == AccumBounds(0, 2)
    assert gruntz(sin(1/x), x, 0) == AccumBounds(-1, 1)
    assert gruntz(exp(-x)*x**5, x, oo) == 0
    with pytest.raises(NotImplementedError, match='^Second argument must be a Symbol$'):
        gruntz(x, x**2, 0)
