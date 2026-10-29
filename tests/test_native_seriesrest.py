"""Limits through rounding functions, modulus squeezes and Piecewise series
at infinity (wave 22 seriesrest).

Every expected value is hard-coded from SymPy 1.14.
"""

import pytest

from symox import (
    AccumBounds, Limit, Mod, Piecewise, Rational, S, Symbol, ceiling, exp,
    floor, limit, oo, sin, symbols,
)

x, y = symbols("x y")
b = Symbol("b", positive=True)
c = Symbol("c", negative=True)
ni = Symbol("n", integer=True)


@pytest.mark.parametrize("expr, x0, d, expected", [
    # floor/ceiling/Mod at an infinite argument: one independent
    # AccumBounds per occurrence, so AB - AB does not cancel
    (floor(x + S.Half) - floor(x), oo, "-", AccumBounds(Rational(-1, 2), Rational(3, 2))),
    (ceiling(x + S.Half) - ceiling(x), oo, "-", AccumBounds(Rational(-1, 2), Rational(3, 2))),
    (floor(x + Rational(1, 3)) - floor(x), oo, "-", AccumBounds(Rational(-2, 3), Rational(4, 3))),
    (floor(x) - ceiling(x), oo, "-", AccumBounds(-2, 0)),
    (floor(x) - x, oo, "-", AccumBounds(-1, 0)),
    (floor(x) - x, -oo, "+", AccumBounds(-1, 0)),
    (ceiling(x) - x, oo, "-", AccumBounds(0, 1)),
    (2*floor(x) - 2*x, oo, "-", AccumBounds(-2, 0)),
    (Mod(x, 3), oo, "-", AccumBounds(0, 3)),
    (floor(x) + Mod(x, 2), oo, "-", oo),
    (Mod(x, 2)/x, oo, "-", S.Zero),
    (floor(x)/x, oo, "-", S.One),
    (floor(x)/x, -oo, "+", S.One),
    (floor(x**2)/x**2, oo, "-", S.One),
    (floor(x/2), oo, "-", oo),
    (x*floor(3/x)/2, 0, "+", Rational(3, 2)),
    (x*ceiling(3/x)/2, 0, "+", Rational(3, 2)),
    # finite points
    (Mod(x, 3), 6, "+", S.Zero),
    (Mod(x, 3), 6, "-", S(3)),
    (floor(x/3), 0, "-", S(-1)),
    (Mod(x, b), ni*b, "+", S.Zero),
    (Mod(x, b), ni*b, "-", b),
    (Mod(x, c), ni*c, "+", c),
    (Mod(x, c), ni*c, "-", S.Zero),
    # a power of a negative base merged into one power (issues 22181, 19739)
    ((-1)**x*2**(-x), oo, "-", S.Zero),
    ((-Rational(1, 4))**x, oo, "-", S.Zero),
    ((-2)**x*3**(-x), oo, "-", S.Zero),
    ((-1)**x*3**(-x)*2**(-x), oo, "-", S.Zero),
    ((-1)**x*exp(-x), oo, "-", S.Zero),
    ((-Rational(1, 2))**(x**2), oo, "-", S.Zero),
])
def test_limit_values(expr, x0, d, expected):
    assert limit(expr, x, x0, d) == expected


def test_floor_n_over_2():
    n = Symbol("n")
    assert limit(floor(n/2), n, oo) is oo


@pytest.mark.parametrize("expr", [
    # the reference leaves these unevaluated: an oscillating sin defeats
    # the leading term and the term-by-term heuristics meet oo - oo
    sin(x) + floor(x) - x,
    x - floor(x) + sin(x),
    sin(x)/x + floor(x) - x,
])
def test_mixed_oscillation_and_rounding_unevaluated(expr):
    assert isinstance(limit(expr, x, oo), Limit)


@pytest.mark.parametrize("expr", [
    # SymPy raises NotImplementedError (sign of log(2) - I*pi undecided);
    # these must not be squeezed to a value
    3*(-1)**x*2**(-x),
    (-1)**(2*x)*2**(-x),
    (-1)**x/x,
    (-Rational(1, 2))**x*x**5,
    (-Rational(1, 2))**x + 1,
])
def test_undecided_negative_base_not_evaluated(expr):
    try:
        r = limit(expr, x, oo)
    except NotImplementedError:
        return
    assert isinstance(r, Limit)


@pytest.mark.parametrize("pw, x0, n", [
    (Piecewise((x, y < 1), (x**5, True)), oo, 3),
    (Piecewise((x, y < 1), (x**5, True)), -oo, 3),
    (Piecewise((x**2, y > 0), (x**3, True)), oo, 4),
    (Piecewise((x, y < 1), (1/x, True)), oo, 3),
])
def test_piecewise_series_at_infinity(pw, x0, n):
    # undecidable conditions keep the Piecewise as it is
    assert pw.series(x, x0, n) == pw
