"""Mul.flatten combines numeric bases sharing one exponent (SymPy 1.14).

Expected values are hard-coded from SymPy 1.14.
"""
import pytest

import symox
from symox import (Float, I, Mul, Pow, Rational, S, Symbol, cos, gamma,
                   gammasimp, pi, sin, sqrt, srepr, symbols)
from symox.integrals.meijerint import meijerint_indefinite

n, m, x, y, s = symbols('n m x y s')
k = Symbol('k', integer=True)
p = Symbol('p', positive=True)
half = Rational(1, 2)


@pytest.mark.parametrize('value, expected', [
    (lambda: 2**n*3**n, '6**n'),
    (lambda: Rational(1, 4)**n*4**n, '1'),
    (lambda: 2**x*half**x, '1'),
    (lambda: 2**(x + 1)*3**(x + 1), '6**(x + 1)'),
    (lambda: sqrt(2)*sqrt(3), 'sqrt(6)'),
    (lambda: 2**n*3**m, '2**n*3**m'),
    (lambda: (-2)**n*3**n, '(-2)**n*3**n'),
    (lambda: 2**(n/2)*3**(n/2), '6**(n/2)'),
    (lambda: Float(2.5)**n*2**n, '5.0**n'),
    (lambda: Float(2.0)**n*Float(3.0)**n, '6.0**n'),
    (lambda: Float(2)**x*Float(0.5)**x, '1.0**x'),
    (lambda: 2**n*Float(3.0)**n, '6.0**n'),
    (lambda: 2**n*3**n*5**n, '30**n'),
    (lambda: 2**n*3**n*5**m*7**m, '35**m*6**n'),
    (lambda: x*2**n*3**n, '6**n*x'),
    (lambda: 3*2**n*3**n, '3*6**n'),
    (lambda: 2**k*3**k, '6**k'),
    (lambda: (-2)**k*(-3)**k, '6**k'),
    (lambda: (-1)**k*2**k, '(-2)**k'),
    (lambda: Float(-2.0)**k*3**k, '(-6.0)**k'),
    (lambda: (-3)**k*Rational(-1, 3)**k, '1'),
    (lambda: 2**I*3**I, '6**I'),
    (lambda: 2**pi*3**pi, '6**pi'),
    (lambda: 2**(2*n)*3**n, '2**(2*n)*3**n'),
    (lambda: 4**n*half**n, '4**n/2**n'),
    (lambda: 2**n*Rational(3, 2)**n, '3**n'),
    (lambda: 2**n*Rational(2, 3)**n, '(4/3)**n'),
    (lambda: Rational(2, 3)**n*Rational(3, 4)**n, '(1/2)**n'),
    (lambda: 2**n*Rational(1, 3)**n, '2**n/3**n'),
    (lambda: 2**(-n)*3**(-n), '6**(-n)'),
    (lambda: half**n*Rational(1, 3)**n, '6**(-n)'),
    (lambda: half**x*half**x, '2**(-2*x)'),
    (lambda: x*half**n, 'x/2**n'),
    (lambda: x*half**(x + 1), '2**(-x - 1)*x'),
    (lambda: 2**x*half**(x + 1), '2**x*2**(-x - 1)'),
    (lambda: 2**x*half**y, '2**x/2**y'),
    (lambda: Rational(3, 2)**x*Rational(2, 3)**x, '1'),
    (lambda: 2**n*3**n*6**(-n), '1'),
    (lambda: Pow(2, n)*Pow(3, n)*Pow(5, -n), '6**n/5**n'),
    (lambda: 2**n*3**n - 6**n, '0'),
    (lambda: 2**p*3**p, '6**p'),
    (lambda: 2**(I*x)*3**(I*x), '6**(I*x)'),
    (lambda: sqrt(2)**n*sqrt(3)**n, '6**(n/2)'),
    (lambda: 2**n*3**n*x**n, '6**n*x**n'),
    (lambda: (2**x*3**x).subs(x, y), '6**y'),
    (lambda: Mul(2**n, 3**n, evaluate=False).doit(), '6**n'),
    (lambda: Rational(1, 2)**(x + 1), '(1/2)**(x + 1)'),
])
def test_numeric_bases_with_equal_exponents(value, expected):
    assert str(value()) == expected


def test_combined_power_structure():
    assert srepr(2**n*3**n*x) == "Mul(Pow(Integer(6), Symbol('n')), Symbol('x'))"
    assert srepr(x*half**n) == "Mul(Pow(Integer(2), Mul(Integer(-1), Symbol('n'))), Symbol('x'))"
    assert srepr(half**n) == "Pow(Rational(1, 2), Symbol('n'))"
    assert (2**n*3**n).func is Pow
    assert (2**n*3**n).args == (S(6), n)


def test_gammasimp_factors_before_expanding_bases():
    e = (sqrt(pi)*(-2**s*gamma(1 - s/2)*gamma(-(s + 1)/2)*gamma((s + 3)/2)
                   + 6**s*s*gamma(-s/2)*gamma((1 - s)/2)*gamma((s + 1)/2)
                   + 6**s*gamma(1 - s/2)*gamma(-(s + 1)/2)*gamma((s + 3)/2))
         / (4*3**s*s*gamma(-s/2)*gamma((1 - s)/2)*gamma(1 - s/2)))
    assert str(gammasimp(e)) == (
        'sqrt(pi)*(2*6**s*cos(pi*s/2)*gamma(1 - s)*gamma(s/2 + 1/2) '
        '- sqrt(pi)*gamma(1 - s/2) + 6**s*sqrt(pi)*gamma(1 - s/2)/2**s)'
        '/(8*3**s*cos(pi*s/2)*gamma(1 - s)*gamma(1 - s/2))')


def test_meijerint_error_after_numeric_base_combination():
    with pytest.raises(NotImplementedError, match='Gammas partially over the strip.'):
        meijerint_indefinite(sin(x)**3/x**3, x)
