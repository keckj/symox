"""Inverse Mellin transform, Meijer G Laplace inversion and besselsimp.

The expected values are the outputs of the reference implementation.
"""
import pytest

import symox as sx
from symox import _core
from symox import (I, Rational, S, besseli, besselj, besselk, bessely, besselsimp, exp, exp_polar,
                   expand_func, gamma, inverse_laplace_transform, inverse_mellin_transform, log, oo, pi,
                   polar_lift, simplify, sqrt)


s, x, t, z = sx.symbols('s x t z')
nu = sx.Symbol('nu')
a, b = sx.symbols('a b', positive=True)
c = sx.Symbol('c')
H = sx.Heaviside
IMT = inverse_mellin_transform


@pytest.mark.parametrize('F, strip, expected', [
    (gamma(s), (0, oo), exp(-x)),
    (gamma(-s), (-oo, 0), exp(-1/x)),
    (1/(s**2 - 1), (-1, None), -x*H(1 - x)/2 - H(x - 1)/(2*x)),
    (1/(s**2 - 1), (None, 1), -x*H(1 - x)/2 - H(x - 1)/(2*x)),
    (gamma(s) + gamma(s - 1), (1, oo), (x + 1)*exp(-x)/x),
    (b**(-s/a)*sx.factorial(s/a)/s, (0, oo), exp(-b*x**a)),
    (sx.factorial(a/b + s/b)/(a + s), (-a, oo), x**a*exp(-x**b)),
    (sx.factorial(s - 1), (0, oo), exp(-x)),
    (24/s**5, (0, oo), log(x)**4*H(1 - x)),
    (pi/(s*sx.sin(pi*s)), (-1, 0), log(x + 1)),
    (pi/(s*sx.sin(pi*s/2)), (-2, 0), log(x**2 + 1)),
    (pi/(s*sx.sin(2*pi*s)), (Rational(-1, 2), 0), log(sqrt(x) + 1)),
    (pi/(s*sx.sin(pi*s)), (0, 1), log(1 + 1/x)),
    (-gamma(s + S.Half)/(sqrt(pi)*s), (Rational(-1, 2), 0), sx.erf(sqrt(x))),
    (pi/sx.cos(pi*s), (0, S.Half), sqrt(x)/(x + 1)),
    (gamma(s)*gamma(1 - s), (0, 1), 1/(x + 1)),
    (gamma(s)*gamma(c - s), (0, None), gamma(c)/(x + 1)**c),
])
def test_inverse_mellin_reference_values(F, strip, expected):
    assert IMT(F, s, x, strip) == expected


def test_inverse_mellin_conditions_and_errors():
    assert IMT(gamma(s), s, x, (0, oo), noconds=False) == (exp(-x), True)
    with pytest.raises(sx.integrals.transforms.MellinTransformStripError, match='Pole inside critical strip'):
        IMT(1/(s**2 - 1), s, x, (-2, 2))
    # no inverse: the transform stays unevaluated
    assert IMT(exp(s), s, x, (0, oo)).has(sx.InverseMellinTransform)
    assert str(IMT(gamma(s)/s, s, x, (0, oo))) == \
        'x*((log(x*exp_polar(I*pi)) - Ei(x*exp_polar(I*pi)))/x + EulerGamma/x) - log(x) - EulerGamma'


@pytest.mark.parametrize('F, expected', [
    (exp(-sqrt(s)), 'exp(-1/(4*t))*Heaviside(t)/(2*sqrt(pi)*t**(3/2))'),
    (exp(-sqrt(s))/s, '(-4*sqrt(pi)*t*erf(1/(2*sqrt(t))) + 4*sqrt(pi)*t)*Heaviside(t)/(4*sqrt(pi)*t)'),
    (log(s)/s, '(-t*log(t) - EulerGamma*t)*Heaviside(t)/t'),
    (exp(-sqrt(s))*exp(-2*s), 'exp(-1/(4*(t - 2)))*Heaviside(t - 2)/(2*sqrt(pi)*(t - 2)**(3/2))'),
    (a**b*(s + sqrt(s**2 + a**2))**(-b)/sqrt(s**2 + a**2),
     'a*a**b*a**(-b - 1)*(-I*sin(pi*b/2) + cos(pi*b/2))*Heaviside(t)*besseli(b, a*t*exp_polar(I*pi/2))'),
    (2**(-s)/s, 'Heaviside(t - log(2))'),
    (gamma(s - 1)/gamma(s + 1), '(1 - exp(-t))*exp(t)*Heaviside(t)'),
])
def test_inverse_laplace_by_meijer_inversion(F, expected):
    assert str(inverse_laplace_transform(F, s, t)) == expected


@pytest.mark.parametrize('F', [sx.atan(1/s), exp(-1/s)/sqrt(s)])
def test_inverse_laplace_without_inversion_stays_held(F):
    assert inverse_laplace_transform(F, s, t).has(sx.InverseLaplaceTransform)


@pytest.mark.parametrize('f, expected', [
    (sqrt(x + 4) - 2, '(((1/sqrt(pi), 0, meijerg(((1, 3/2), ()), ((1,), (0,)), x/4)),), True)'),
    (sqrt(x + 4) + 2,
     '(((exp_polar(I*pi)/sqrt(pi), 0, meijerg(((1, 3/2), ()), ((1,), (0,)), x*exp_polar(-2*I*pi)/4)),), True)'),
    (1/(sqrt(x + 4) + 2), '(((-exp_polar(I*pi)/(4*sqrt(pi)), 0, meijerg(((0, 1/2), ()), ((0,), (-1,)), x/4)),), True)'),
    ((sqrt(x**2 + 9) + 3)**b,
     '(((-3**b*b/(2*sqrt(pi)), 0, meijerg(((b/2 + 1/2, b/2 + 1), ()), ((0,), (b,)), x**2/9)),), True)'),
    ((x**c - 2**c)/(x - 2), '(((2**(c - 1)*sin(pi*polar_lift(c))/pi, 0, meijerg(((0, c), ()), ((0, c), ()), x/2)),), True)'),
])
def test_meijer_table_radical_and_quotient_entries(f, expected):
    assert str(tuple(_core._meijerint_debug('rewrite_single', f, x))) == expected


def test_definite_integrals_through_radical_entries():
    assert sx.integrate((sqrt(x + 4) - 2)*exp(-x), (x, 0, oo)) == sqrt(pi)*exp(4)*sx.erfc(2)/2


@pytest.mark.parametrize('e, expected', [
    (besseli(b, a*t*exp_polar(I*pi)), exp(I*pi*b)*besseli(b, a*t)),
    (besseli(b, a*t*exp_polar(I*pi/2)), exp(I*pi*b/2)*besselj(b, a*t)),
    (exp(-I*pi*b/2)*H(t)*besseli(b, a*t*exp_polar(I*pi/2)), H(t)*besselj(b, a*t)),
    (besselj(b, a*t*exp_polar(I*pi/2)), exp(I*pi*b/2)*besseli(b, a*t)),
    (besselj(nu, z*polar_lift(-1)), exp(I*pi*nu)*besselj(nu, z)),
    (besseli(nu, z*polar_lift(-I)), exp(-I*pi*nu/2)*besselj(nu, z)),
    (besseli(Rational(-1, 2), z), sqrt(2)*sx.cosh(z)/(sqrt(pi)*sqrt(z))),
    (z*besseli(0, z) + z*besseli(2, z)/2 + besseli(1, z), 3*z*besseli(0, z)/2),
    (besselj(4, z), (z**3*besselj(0, z) - 8*z**2*besselj(1, z) - 24*z*besselj(0, z) + 48*besselj(1, z))/z**3),
    (besselj(b, I*z), besselj(b, I*z)),
    (besselj(0, z) + besselj(1, z) + besselj(2, z), (z + 2)*besselj(1, z)/z),
    (besselj(3, z)*z - 4*besselj(2, z) + z*besselj(1, z), 0),
    (bessely(b, z) + bessely(b + 1, z) + bessely(b + 2, z), (2*b + z + 2)*bessely(b + 1, z)/z),
    (besselk(3, z), (z**2*besselk(1, z) + 4*z*besselk(0, z) + 8*besselk(1, z))/z**2),
])
def test_besselsimp_and_simplify(e, expected):
    assert besselsimp(e) == expected
    assert simplify(e) == expected


def test_expand_func_bessel_recurrence():
    # the unevaluated nesting of the reference recurrence
    assert str(expand_func(besselj(4, z))) == \
        'besselj(0, z) + (-6*besselj(1, z) + 6*(-4*besselj(0, z) + 8*besselj(1, z)/z)/z)/z - 2*besselj(1, z)/z'
    assert str(expand_func(besseli(Rational(7, 2), z))) == \
        'besseli(-1/2, z) + (-5*besseli(1/2, z) - 5*(-3*besseli(-1/2, z) + 3*besseli(1/2, z)/z)/z)/z - besseli(1/2, z)/z'
    assert sx.expand(expand_func(besselj(4, z)) - besselsimp(besselj(4, z))) == 0
    # a symbolic order is not reduced
    assert expand_func(besselk(nu + 3, z)) == besselk(nu + 3, z)
