"""Certified derangement remainders, independently of SymPy's wrong zero."""

import pytest

from symox import (
    AccumBounds, E, I, Function, Limit, Rational, Symbol, exp, factorial, gamma, limit, log, oo,
    sin, sqrt, subfactorial,
)

x = Symbol("x")


@pytest.mark.parametrize("argument", [x, 2*x, x + 3, x**2, exp(x)])
def test_normalized_subfactorial(argument):
    # !a = Gamma(a + 1)/e + r(a), |r(a)| <= 1/(a + 1) for real a > -1.
    assert limit(subfactorial(argument)/factorial(argument), x, oo) == exp(-1)


@pytest.mark.parametrize("expr, expected", [
    (subfactorial(x)/gamma(x + 1), exp(-1)),
    (subfactorial(x + 1)/(x*factorial(x)), exp(-1)),
    (subfactorial(x) - factorial(x)/E, 0),
    (subfactorial(2*x)/factorial(2*x) + subfactorial(x)/factorial(x), 2/E),
    (I*subfactorial(x)/factorial(x), I/E),
    ((subfactorial(x) - factorial(x)/E)/x, 0),
    ((subfactorial(x) - factorial(x)/E)*x**Rational(1, 2), 0),
    ((subfactorial(x) - factorial(x)/E)*x/log(x), 0),
])
def test_affine_certified_remainders(expr, expected):
    assert limit(expr, x, oo) == expected


def test_finite_point_mapped_to_positive_infinity():
    expr = subfactorial(1/x)/factorial(1/x)
    assert limit(expr, x, 0, "+") == exp(-1)


@pytest.mark.parametrize("expr", [
    x*(subfactorial(x) - factorial(x)/E),
    x**2*(subfactorial(x) - factorial(x)/E),
    factorial(x)*(subfactorial(x) - factorial(x)/E),
    1/(subfactorial(x) - factorial(x)/E),
    subfactorial(x + I)/factorial(x + I),
    subfactorial(x + sqrt(-x)) - factorial(x + sqrt(-x))/E,
    subfactorial(x - sqrt(-x)) - factorial(x - sqrt(-x))/E,
    subfactorial(log(-x)) - factorial(log(-x))/E,
    subfactorial(x + log(-x)) - factorial(x + log(-x))/E,
    subfactorial(x + Function("f")(x)) - factorial(x + Function("f")(x))/E,
    subfactorial(x + Symbol("a"))/factorial(x + Symbol("a")),
])
def test_uncertified_remainders_remain_undecided(expr):
    assert isinstance(limit(expr, x, oo), Limit)


def test_oscillating_argument_does_not_become_positive_infinity():
    result = limit(subfactorial(sin(x))/factorial(sin(x)), x, oo)
    # Preserve the existing bounded-oscillation fallback.
    assert result.has(AccumBounds)


def test_integer_subsequence_independently_rules_out_zero():
    # Exact alternating-series bounds bracket exp(-1). In particular,
    # every n >= 2 has !n/n! >= 1/3, ruling out SymPy 1.14's zero limit.
    for n in (2, 3, 4, 10, 20, 40):
        ratio = subfactorial(n)/factorial(n)
        partial = sum(Rational((-1)**k, int(factorial(k))) for k in range(n + 1))
        assert ratio == partial
        assert ratio >= Rational(1, 3)
