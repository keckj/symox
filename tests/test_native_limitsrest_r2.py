"""Symbolic powers in Gruntz, nseries along cdir, limits SymPy evaluates by
heuristics, Order forms and the AccumBounds series error, against values
computed with SymPy 1.14."""
import pytest

from symox import *  # noqa: F401,F403
from symox import _core

del test  # noqa: F821

x, y = symbols('x y')
p = Symbol('p', positive=True)
r = Symbol('r', real=True)


@pytest.mark.parametrize('expr, x0, dir, expected', [
    ((1 + x**log(3))**(1/x), 0, '+', 1),
    ((1 + x**log(3))**(1/x), 0, '-', 1),
    (x*log(1 + x**(-log(3))), oo, '-', 0),
    ((1 + x**sqrt(2))**(1/x), 0, '+', 1),
    ((1 + x**(1/log(3)))**(1/x), 0, '+', oo),
    ((1 + x**log(2))**(1/x), 0, '+', oo),
    (sin(x**log(3))/x, 0, '+', 0),
    (log(1 + x**pi)/x**3, 0, '+', 0),
    ((exp(x**sqrt(3)) - 1)/x, 0, '+', 0),
    (x**log(3)/x, 0, '+', 0),
])
def test_symbolic_power_limits(expr, x0, dir, expected):
    assert limit(expr, x, x0, dir) == expected


@pytest.mark.parametrize('expr, cdir, expected', [
    (floor(x), -1, lambda: Integer(-1)),
    (ceiling(x), -1, lambda: Integer(0)),
    (frac(x), -1, lambda: x + 1),
    (sign(x), -1, lambda: Integer(-1)),
    (floor(x) + 2, -1, lambda: Integer(1)),
    (x*floor(x), -1, lambda: -x + O(x**3)),
    (x*(floor(x) + 1), -1, lambda: Integer(0)),
    (x**5*floor(x), -1, lambda: O(x**3)),
    (x*sign(x) + 1, -1, lambda: 1 - x),
    (floor(x)/x, -1, lambda: -1/x + O(x**3)),
    (floor(-x), -1, lambda: Integer(0)),
    (ceiling(-x**3), -1, lambda: Integer(1)),
    (frac(-x), -1, lambda: -x),
    (sign(x**2 - x), -1, lambda: Integer(1)),
    (atan(1/x), -1, lambda: -pi/2 - x + O(x**3)),
    (x*atan(1/x), -1, lambda: -pi*x/2 - x**2 + O(x**3)),
    (acot(x), -1, lambda: -pi/2 - x + O(x**3)),
    (cos(x)*sign(x), -1, lambda: -1 + x**2/2 + O(x**3)),
    (frac(sin(x) + 5), -1, lambda: 1 + x + O(x**3)),
    (log(-x), -1, lambda: log(x) - I*pi),
    (log(-x), I, lambda: log(x) - I*pi),
    (log(x**2 - x**3), -1, lambda: -2*I*pi + 2*log(x) - x - x**2/2 + O(x**3)),
    (log(-x + I*x**2), -1, lambda: -I*pi + log(x) - I*x + x**2/2 + O(x**3)),
    (log(-1 - I*x), -1, lambda: I*pi + I*x + x**2/2 + O(x**3)),
    (log(-1 + I*x), -1, lambda: -I*pi - I*x + x**2/2 + O(x**3)),
    (log(x*(I - 1)), -1, lambda: log(x) + log(sqrt(2)) - 5*I*pi/4),
])
def test_nseries_cdir(expr, cdir, expected):
    assert expr.nseries(x, 0, 3, cdir=cdir) == expected()


def test_nseries_cdir_not_sure_of_sign():
    with pytest.raises(NotImplementedError):
        floor(x).nseries(x, 0, 3, cdir=I)
    with pytest.raises(NotImplementedError):
        (x*floor(x)).nseries(x, 0, 3, cdir=I)


def test_heuristic_limits():
    assert limit((Abs(x + y) - Abs(x - y))/(2*x), x, 0) == sign(y)
    assert limit((Abs(x + y) - Abs(x - y))/(2*x), x, 0, '-') == sign(y)
    assert limit((Abs(x + r) - Abs(x - r))/(2*x), x, 0) == sign(r)
    assert limit(Abs(x + y), x, 0) == y*sign(y)
    assert limit(sign(x + y), x, 0) == sign(y)
    assert limit(Abs(x + y)/(x + y), x, 0) == sign(y)
    g = gamma(-oo)
    assert limit(factorial(x), x, -oo) == g
    assert limit(gamma(x), x, -oo) == g
    assert limit(factorial(2*x), x, -oo) == g
    assert limit(factorial(-x**2), x, oo) == g
    assert limit(gamma(1/x), x, 0, '-') == g
    assert limit(gamma(x) + 1, x, -oo) == 1 + g
    assert limit(1/gamma(x), x, -oo) == 1/g
    assert limit(factorial(x**2), x, -oo) == oo


def test_limit_with_order_at_another_point():
    assert limit(1 + O(x**2), x, 1) == 1
    assert limit(1 + O(x**2), x, oo) == 1
    assert limit(1 + O(x**2), x, -1) == 1
    assert limit(1 + O(x**2), x, y) == 1
    assert limit(x + O(x**2), x, 1) == 1
    assert limit(1/x + O(x), x, 2) == Rational(1, 2)
    assert limit(sin(x) + O(x**2), x, 1) == sin(1)
    assert limit(O(x**2) + 1/(x - 1), x, 1) == oo
    assert limit(1 + O(x**2), x, 0) == 1
    assert isinstance(limit(1 + O(x**2, (x, oo)), x, 0), Limit)


def test_series_abs_symbolic_lead():
    assert Abs(x + y).series(x, 0, 3) == x*sign(y) + y*sign(y)
    assert Abs(y*x + x**2).series(x, 0, 3) == x**2*sign(y) + x*y*sign(y)
    assert Abs(x + p).series(x, 0, 3) == p + x
    assert Abs(x + r).series(x, 0, 3) == r*sign(r) + x*sign(r)
    assert Abs(sin(x) + y).nseries(x, 0, 3) == y*sign(y) + x*sign(y) + O(x**3)
    assert sign(x + y).series(x, 0, 3) == sign(y)
    assert Abs(x*y).series(x, 0, 3) == x*Abs(y)
    assert sign(x**2*y).series(x, 0, 3) == sign(y)
    assert limit(Abs(x*y)/x, x, 0) == Abs(y)
    assert limit(Abs(x*y)/x, x, 0, '-') == -Abs(y)
    assert Abs(x*y).series(x, 0, 3, dir='-') == -x*Abs(y)
    assert limit(sign(x*y), x, 0, '-') == -sign(y)
    assert limit(Abs(x**2*y)/x**2, x, 0, '-') == Abs(y)


def test_order_forms():
    assert O(x**2).subs(x, cos(x)) == O((x - pi/2)**2, (x, pi/2))
    assert O(x**2).subs(x, cos(x) + 1) == O((x - pi)**4, (x, pi))
    assert O(x**2).diff(x) == O(x)
    assert str(O(x**2).diff(x, 2)) == 'O(1)'
    assert O(x).diff(x) == O(1, x)
    assert O(x**3*y, x, y).diff(x) == O(x**2*y, x, y)
    assert O(x**2).diff(y) == 0
    assert O(x**2, (x, oo)).diff(x) == O(x, (x, oo))
    assert O(x*log(x)).diff(x) == O(log(x))
    assert diff(sin(x) + O(x**3), x) == cos(x) + O(x**2)
    assert O(nan) is nan
    assert Order(log(x) + 1/log(x)) == O((log(x)**2 + 1)/log(x))
    assert Order(1 + 1/log(x)) == O((log(x) + 1)/log(x))
    assert Order(log(x) - 1/log(x)) == O((log(x) - 1)*(log(x) + 1)/log(x))
    assert Order(2*log(x) + 3/log(x)) == O((2*log(x)**2 + 3)/log(x))
    assert Order(x + 1/log(x)) == O(1/log(x))
    assert Order(log(x)**2 + log(x)) == O((log(x) + 1)*log(x))
    assert Order(x*log(x) + x) == O(x*(log(x) + 1))
    assert Order(log(x)/x + 1/x) == O((log(x) + 1)/x)
    assert Order(log(x) + 1) == O(log(x))
    assert Order(log(x)**2 + 1) == O(log(x)**2)


@pytest.mark.parametrize('expr, x0', [
    (floor(x), oo), (ceiling(x), oo), (floor(x), -oo), (floor(x)/x, oo),
    (floor(x) + x, oo), (floor(x)*sin(1/x), oo), (floor(x**2), oo),
    (floor(1/x), 0), (floor(1/x) + sin(x), 0), (floor(-1/x), 0),
    (floor(1/x + 1), 0), (2*ceiling(1/x), 0), (floor(log(x)), 0),
])
def test_series_accum_bounds_error(expr, x0):
    with pytest.raises(ValueError, match='Only real AccumulationBounds are supported'):
        expr.series(x, x0)


def test_series_accum_bounds_other_errors():
    with pytest.raises(ValueError, match='Only real AccumulationBounds'):
        floor(1/x).series(x, 0, dir='-')
    with pytest.raises(NotImplementedError):
        (1/floor(x)).series(x, oo)
    with pytest.raises(NotImplementedError):
        (floor(1/x)**2).series(x, 0)
    assert _core._series_accum_error(exp(floor(x)), x, oo, '+') is None
    assert _core._series_accum_error(floor(1/p), p, 0, '+') is None
