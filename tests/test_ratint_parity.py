"""Rational integration, RootSum and supporting polynomial parity.

Every expected value below is the output of SymPy 1.14.
"""

import pytest

from symox import (
    EX,
    Lambda,
    Poly,
    QQ,
    RootSum,
    ZZ,
    expand,
    gcd,
    integrate,
    log,
    sqrt,
    srepr,
    symbols,
)
from symox.integrals.rationaltools import ratint
from symox.polys.polyerrors import MultivariatePolynomialError

x, t, a, b, y = symbols('x t a b y')

R2 = ('RootSum(_t**3*(27*a**2 + 4) - 3*_t - 1, Lambda(_t, _t*log(x + '
      '(-54*_t**2*a**2 - 8*_t**2 + 27*_t*a**2 + 4*_t + 4)/(9*a))))')


def test_ratint_symbolic_cubic_is_rootsum():
    assert str(ratint(1/(x**3 + x + a), x)) == R2


def test_integrate_symbolic_cubic_is_rootsum():
    assert str(integrate(1/(x**3 + x + a), x)) == R2


def test_ratint_symbolic_quadratic():
    assert str(ratint(1/(x**2 + a), x)) == (
        '-sqrt(-1/a)*log(-a*sqrt(-1/a) + x)/2 + sqrt(-1/a)*log(a*sqrt(-1/a) + x)/2')


def test_rootsum_requires_univariate_polynomial():
    with pytest.raises(MultivariatePolynomialError):
        RootSum(x**3 + x + a, Lambda(t, t*log(x + t)))


def test_rootsum_keeps_irreducible_quadratic():
    assert str(RootSum(x**2 + 1, Lambda(t, log(y + t)))) == (
        'RootSum(x**2 + 1, Lambda(t, log(t + y)))')


def test_rootsum_srepr():
    assert srepr(RootSum(x**3 + x + 7, Lambda(t, t*log(y + t)))) == (
        "RootSum(Add(Pow(Symbol('x'), Integer(3)), Symbol('x'), Integer(7)), "
        "Lambda(Tuple(Symbol('t')), Mul(Symbol('t'), log(Add(Symbol('t'), "
        "Symbol('y'))))), Symbol('x'))")


def test_rootsum_derivative_is_rational():
    expr = RootSum(x**3 + x + 1, Lambda(t, t*log(y + t)))
    assert str(expr.diff(y)) == '(2*y - 3)/(y**3 + y - 1)'


def test_rootsum_doit_on_rational_roots_of_quadratic():
    expr = RootSum(x**2 - 2, Lambda(t, t*log(y + t)))
    assert str(expr.doit()) == '-sqrt(2)*log(y - sqrt(2)) + sqrt(2)*log(y + sqrt(2))'


@pytest.mark.parametrize('expr, expected', [
    (-31*(x + 1/(31*y))**3, '-31*x**3 - 3*x**2/y - 3*x/(31*y**2) - 1/(961*y**3)'),
    (x + (y + 1/(3*x))**2, 'x + y**2 + 2*y/(3*x) + 1/(9*x**2)'),
])
def test_expand_distributes_rational_terms_over_denominators(expr, expected):
    assert str(expand(expr)) == expected


@pytest.mark.parametrize('expr, domain', [
    (3*t + a, EX),
    (t**2/a + 1, EX),
    (sqrt(2)*t + 1, EX),
    (2*t + 1, ZZ),
    (t/2 + 3, QQ),
    ((a + b)*t**2 + t, EX),
])
def test_non_composite_domain(expr, domain):
    assert Poly(expr, t, composite=False).domain == domain


def test_non_composite_field_domain():
    assert str(Poly(3*t + a, t, composite=False, field=True)) == "Poly(3*t + a, t, domain='EX')"


def test_polynomial_ring_ground_gcd():
    f = Poly(2*a*t + 2*a, t, domain='QQ[a]')
    g = Poly(4*a*t + 4*a, t, domain='QQ[a]')
    assert str(gcd(f, g)) == "Poly(2*a*t + 2*a, t, domain='QQ[a]')"


def test_polynomial_ring_ground_content_and_primitive():
    p = Poly(2*a*t + 4*a, t, domain='QQ[a]')
    assert str(p.content()) == '2*a'
    assert str(p.primitive()) == "(2*a, Poly(t + 2, t, domain='QQ[a]'))"
