"""Limit parity cases fixed in wave 13 (limrest).

Every expected value is hard-coded from SymPy 1.14.
"""
from symox import *
del test

import pytest

x = Symbol('x')
z = Symbol('z')


CASES = [
    # Bessel functions at infinity (Hankel asymptotic expansions)
    (lambda: limit(besselj(0, x), x, oo), S.Zero),
    (lambda: limit(besselj(1, x), x, oo), S.Zero),
    (lambda: limit(besselj(1, x), x, -oo), S.Zero),
    (lambda: limit(bessely(0, x), x, oo), S.Zero),
    (lambda: limit(bessely(1, x), x, oo), S.Zero),
    (lambda: limit(bessely(1, x), x, -oo), S.Zero),
    (lambda: limit(besseli(0, x), x, oo), oo),
    (lambda: limit(besseli(1, x), x, oo), oo),
    (lambda: limit(besseli(1, x), x, -oo), -oo),
    (lambda: limit(besseli(1, x), x, I*oo), S.Zero),
    (lambda: limit(besseli(1, x), x, -I*oo), S.Zero),
    (lambda: limit(besselk(0, x), x, oo), S.Zero),
    (lambda: limit(besselk(1, x), x, oo), S.Zero),
    (lambda: limit(besseli(1, 3*x)/(x*besseli(1, x)**3), x, oo), 2*sqrt(3)*pi/3),
    (lambda: limit(besseli(1, 3*x)/(x*besseli(1, x)**3), x, -oo), -2*sqrt(3)*pi/3),
    (lambda: limit(besseli(0, x + 1)/besseli(0, x), x, oo), E),
    # Bessel and LambertW Taylor series at zero
    (lambda: limit((besselj(z, x)/x**z).subs(z, 7), x, 0), Rational(1, 645120)),
    (lambda: limit(besselj(7, x)/x**7, x, 0), Rational(1, 645120)),
    (lambda: limit(besseli(1, x)/x, x, 0), Rational(1, 2)),
    (lambda: LambertW(x).limit(x, 0), S.Zero),
    (lambda: limit(LambertW(x)/x, x, 0), S.One),
    (lambda: limit((LambertW(x) - x)/x**2, x, 0), S.NegativeOne),
    (lambda: limit(-LambertW(-log(x))/log(x), x, 1), S.One),
    (lambda: limit(x*diff(sqrt(x)*besseli(1, sqrt(8*x)), x, x)/(sqrt(x)*besseli(1, sqrt(8*x))), x, 0), Integer(2)),
    # asech/acsch
    (lambda: limit(log(x) + asech(x), x, 0, '+'), log(2)),
    (lambda: limit(log(x) + asech(x), x, 0, '-'), log(2) + 2*I*pi),
    (lambda: limit(log(x) + asech(x), x, oo), oo),
    (lambda: limit(log(x) + acsch(x), x, 0, '+'), log(2)),
    (lambda: limit(log(x) + acsch(x), x, 0, '-'), -oo),
    (lambda: limit(log(x) + acsch(x), x, oo), oo),
    (lambda: limit(asech(x), x, 0, '+'), oo),
    (lambda: limit(acsch(x), x, 0, '+'), oo),
    (lambda: limit(asech(x), x, oo), I*pi/2),
    (lambda: limit(acsch(x), x, oo), S.Zero),
    (lambda: limit(asech(sin(x)/x), x, 0, '+-'), S.Zero),
    (lambda: limit(asech(exp(1/x)), x, 0, '+'), I*pi/2),
    (lambda: limit(asech(exp(1/x)), x, 0, '-'), oo),
    (lambda: limit(acsch(exp(1/x)), x, 0, '+'), S.Zero),
    (lambda: limit(acsch(exp(1/x)), x, 0, '-'), oo),
    # inverse functions at an infinite argument
    (lambda: limit(asin(exp(x)), x, oo, '-'), -oo*I),
    (lambda: limit(acos(exp(x)), x, oo, '-'), oo*I),
    (lambda: limit(asec(exp(x)), x, oo, '-'), pi/2),
    (lambda: limit(acsc(exp(x)), x, oo, '-'), S.Zero),
    (lambda: limit(asin(exp(1/x)), x, 0, '+'), -oo*I),
    (lambda: limit(acos(exp(1/x)), x, 0, '+'), oo*I),
    (lambda: limit(asec(sin(x)/x), x, 0, '+-'), S.Zero),
    (lambda: limit(asec(exp(1/x)), x, 0, '+'), pi/2),
    (lambda: limit(asec(exp(1/x)), x, 0, '-'), oo*I),
    (lambda: limit(acsc(exp(1/x)), x, 0, '+'), S.Zero),
    (lambda: limit(acsc(exp(1/x)), x, 0, '-'), -oo*I),
    (lambda: limit(asinh(x) - log(2*x), x, oo), S.Zero),
    (lambda: limit(acosh(x) - log(2*x), x, oo), S.Zero),
    (lambda: limit(asinh(x) + log(-2*x), x, -oo), S.Zero),
    (lambda: limit(asec(x), x, oo), pi/2),
    (lambda: limit(acsc(x), x, oo), S.Zero),
    (lambda: limit(asec(x), x, 2), pi/3),
    (lambda: limit(acos(2*x), x, oo), oo*I),
    (lambda: limit(asin(x)/log(x), x, oo), -I),
    # trigonometric pole whose location carries pi (issue 14411)
    (lambda: limit(3*sec(4*pi*x - x/3), x, 3*pi/(24*pi - 2)), -oo),
    # branch side of a sqrt sitting on its cut (issue 23752)
    (lambda: limit(sqrt(-I*x**2 + x - 3), x, 0, '+'), -sqrt(3)*I),
    (lambda: limit(sqrt(-I*x**2 + x - 3), x, 0, '-'), -sqrt(3)*I),
    (lambda: limit(sqrt(-I*x**2 + I*x - 3), x, 0, '+'), sqrt(3)*I),
    (lambda: limit(sqrt(-I*x**2 + I*x - 3), x, 0, '-'), -sqrt(3)*I),
    # directed infinity reduced to oo
    (lambda: limit(x*sin(1/x), x, I*oo), S.One),
    (lambda: limit(1/x, x, -I*oo), S.Zero),
]


@pytest.mark.parametrize('i', range(len(CASES)))
def test_limrest_case(i):
    f, expected = CASES[i]
    assert f() == expected


def test_issue_6052():
    g = hyperexpand(meijerg((), (), (1,), (0,), -x))
    assert limit(g, x, 0, '+-') == 0
    assert limit(g, x, oo) == -oo
