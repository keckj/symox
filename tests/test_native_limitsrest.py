"""Series front end (cdir, non-symbol variable, logx) and Gruntz constant
signs, against values computed with SymPy 1.14."""
import pytest

from symox import *  # noqa: F401,F403
from symox import AccumBounds, _core

del test  # noqa: F821

x, a, b, w, y = symbols('x a b w y')


@pytest.mark.parametrize('expr, expected', [
    (lambda: Abs(x).series(x, 0, cdir=-1), lambda: -x),
    (lambda: Abs(x).series(x, 0, cdir=1), lambda: x),
    (lambda: Abs(x).series(x, 0, cdir=I), lambda: -I*x),
    (lambda: log(x).series(x, 0, cdir=-1), lambda: log(-x) + I*pi),
    (lambda: log(x).series(x, 0, 3, cdir=I), lambda: log(-I*x) + I*pi/2),
    (lambda: log(x).series(x, 0, 3, cdir=-I), lambda: log(I*x) - I*pi/2),
    (lambda: sqrt(x).series(x, 0, cdir=-1), lambda: I*sqrt(-x)),
    (lambda: frac(x).series(x, cdir=-1), lambda: x + 1),
    (lambda: floor(x + 2).series(x, 0, cdir=-1), lambda: Integer(1)),
    (lambda: sign(x).series(x, 0, cdir=-1), lambda: Integer(-1)),
    (lambda: frac(sin(x) + 5).series(x, 0, 4, cdir=-1), lambda: 1 + x - x**3/6 + O(x**4)),
    (lambda: sin(x).series(x, 0, 4, cdir=I), lambda: x - x**3/6 + O(x**4)),
    (lambda: (1/x).series(x, oo, 3, cdir=-1), lambda: 1/x),
    (lambda: (1/(x + 1)).series(x, -oo, 3, cdir=1), lambda: -1/x**2 + 1/x + O(x**(-3), (x, oo))),
])
def test_series_cdir(expr, expected):
    assert expr() == expected()


def test_series_non_symbol_variable():
    r = series(exp(a*b), a*b, n=8)
    assert r.removeO() == sum((a*b)**k/factorial(k) for k in range(8))
    assert exp(x**2).series(x**2, 0, 3) == 1 + x**2 + x**4/2 + O(x**6)
    assert (x**2).series(x**2, 0, 3) == x**2


def test_series_logx_substituted_once():
    assert log(-log(x)).nseries(x, logx=-x) == log(x)
    assert log(log(x)).nseries(x, logx=-x) == log(-x)
    assert log(sin(x)).series(x, 0, 3, logx=w) == -x**2/6 + w + O(x**3)
    assert exp(x*log(x)).series(x, 0, 3, logx=w) == 1 + w*x + w**2*x**2/2 + O(x**3)
    assert (x**y).series(x, 0, 3, logx=w) == exp(w*y)


def test_series_direction_binding():
    assert _core._series_direction(Integer(0)) is None
    assert _core._series_direction(Integer(-2)) == '-'
    assert _core._series_direction(Integer(3)) == '+'
    assert _core._series_direction(1 + I) == 'I'
    with pytest.raises(NotImplementedError):
        _core._series_direction(y)


@pytest.mark.parametrize('expr, expected', [
    (2**(-x)*sqrt(4**(x + 1) + 1), 2),
    (4**x/2**(2*x), 1),
    (8**x/2**(3*x + 1), Rational(1, 2)),
    (sqrt(9**x + 3**x)/3**x, 1),
])
def test_gruntz_perfect_power_logs(expr, expected):
    # the exponent difference log(4)/log(2) - 2 is an exact zero
    assert _core._limit(expr, x, oo, '+') == expected
    assert limit(expr, x, oo) == expected


@pytest.mark.parametrize('expr', [sin(x)**log(x), sin(x)**cos(x), tan(x)**(2**sqrt(pi))])
def test_bounds_power_with_possibly_negative_base_held(expr):
    assert isinstance(limit(expr, x, oo), Limit)


def test_bounds_power_still_evaluated():
    assert limit(sin(x)**2, x, oo) == AccumBounds(0, 1)
    assert limit(frac(x)**x, x, oo) == AccumBounds(0, oo)


def test_abs_of_real_power():
    # Abs((-x/(x + 1))**x) is (x/(x + 1))**x on the positive variable
    assert limit(Abs((-x/(x + 1))**x), x, oo) == exp(-1)


def test_order_bounded_absorbs():
    assert O(x) + O(1) == O(1)
    assert O(1) + O(x) == O(1)
    assert x + O(1) + O(x) == O(1)
    assert O(1/x) + O(1) == O(1/x)
    assert O(x**2) + O(x) == O(x)
    assert O(x) + O(y) == O(y) + O(x)


def test_bounds_rational_power_times_vanishing_factor():
    assert limit(log(2 + sqrt(atan(x))*sqrt(sin(1/x))), x, 0) == log(2)
