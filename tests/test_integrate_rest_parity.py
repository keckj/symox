"""Remaining integration gaps closed against the reference SymPy 1.14.

The expected printed forms were produced by SymPy 1.14 on the same inputs,
except where a test documents an intentional difference.
"""
import pytest
from symox import (
    FourierTransform, Integral, Piecewise, Rational, Symbol, cos, exp,
    fourier_transform, integrate, log, oo, pi, sin, sqrt, symbols, tan, atan,
)
from symox.integrals.risch import risch_integrate

x, y, k = symbols('x y k')
a = Symbol('a')


# --- Fourier transform of a Gaussian with no assumption on its parameter ---
# The transform integral is computed under conditions, but the Gaussian
# itself is not proven integrable: the reference raises internally and
# returns the held transform with and without noconds (the value would be
# invalid for re(a) <= 0).

def test_fourier_gaussian_without_assumptions_stays_held():
    expected = FourierTransform(exp(-a*x**2), x, k)
    assert fourier_transform(exp(-a*x**2), x, k) == expected
    assert fourier_transform(exp(-a*x**2), x, k, noconds=False) == expected


def test_fourier_gaussian_with_positive_parameter():
    p = Symbol('p', positive=True)
    assert str(fourier_transform(exp(-p*x**2), x, k)) == 'sqrt(pi)*exp(-pi**2*k**2/p)/sqrt(p)'


# --- floor terms making atan(c*tan(u)) antiderivatives continuous ---

FLOOR_CASES = [
    ('integrate(1/(cos(x) + 2), x)',
     '2*sqrt(3)*(atan(sqrt(3)*tan(x/2)/3) + pi*floor((x/2 - pi/2)/pi))/3'),
    ('integrate(1/(cos(x) + 2), x, heurisch=True)',
     '2*sqrt(3)*(atan(sqrt(3)*tan(x/2)/3) + pi*floor((x/2 - pi/2)/pi))/3'),
    ('integrate(1/(2 + sin(x)), x)',
     '2*sqrt(3)*(atan(2*sqrt(3)*tan(x/2)/3 + sqrt(3)/3) + pi*floor((x/2 - pi/2)/pi))/3'),
    ('integrate(1/(1 + 3*cos(x)**2), x)',
     'atan(2*tan(x/2) - sqrt(3))/2 + atan(2*tan(x/2) + sqrt(3))/2 + pi*floor((x/2 - pi/2)/pi)'),
    ('integrate(1/(1 + cos(x)**2), x)',
     'sqrt(2)*(atan(sqrt(2)*tan(x/2) - 1) + pi*floor((x/2 - pi/2)/pi))/2'
     ' + sqrt(2)*(atan(sqrt(2)*tan(x/2) + 1) + pi*floor((x/2 - pi/2)/pi))/2'),
    ('integrate(atan(tan(x)), x)', '(atan(tan(x)) + pi*floor((x - pi/2)/pi))**2/2'),
    ('integrate(1/(cos(x) + 2), (x, 0, 2*pi))', '2*sqrt(3)*pi/3'),
    ('integrate(1/(cos(x) + 2), (x, 0, 7))', '2*sqrt(3)*atan(sqrt(3)*tan(7/2)/3)/3 + 2*sqrt(3)*pi/3'),
    ('integrate(1/(5 + 4*sin(x)), (x, 0, 2*pi))', '2*pi/3'),
    ('integrate(1/(cos(x) + 2), x, manual=True)', 'Integral(1/(cos(x) + 2), x)'),
    ('integrate(1/(cos(x) + 2), x, risch=True)', 'Integral(1/(cos(x) + 2), x)'),
    ('integrate(1/(cos(x) + 2), x, meijerg=True)', 'Integral(1/(cos(x) + 2), x)'),
]


@pytest.mark.parametrize('call, expected', FLOOR_CASES)
def test_atan_floor_terms(call, expected):
    assert str(eval(call)) == expected


def test_piecewise_floor_terms_added_once():
    # Intentional difference: the reference adds the floor terms a second
    # time to the antiderivative of a Piecewise (2*pi*floor(...)), which
    # makes it discontinuous and its definite value on (-1, 7) wrong
    # (4*sqrt(3)*pi/3 + ...). The terms are added once here.
    f = Piecewise((1/(cos(x) + 2), x > 0), (0, True))
    assert str(integrate(f, x)) == (
        'Piecewise((0, x <= 0), (2*sqrt(3)*(atan(sqrt(3)*tan(x/2)/3)'
        ' + pi*floor((x/2 - pi/2)/pi))/3 + 2*sqrt(3)*pi/3, True))')
    assert str(integrate(f, (x, -1, 7))) == '2*sqrt(3)*atan(sqrt(3)*tan(7/2)/3)/3 + 2*sqrt(3)*pi/3'


# --- RootSum of the Risch logarithmic part on exp/log towers ---

ROOTSUM_CASES = [
    ('integrate(exp(x)/(1 + exp(2*x)), x)', 'RootSum(4*_z**2 + 1, Lambda(_i, _i*log(2*_i + exp(x))))'),
    ('integrate(exp(x)/(1 + exp(2*x)), x, risch=True)', 'RootSum(4*_z**2 + 1, Lambda(_i, _i*log(2*_i + exp(x))))'),
    ('risch_integrate(exp(x)/(1 + exp(2*x)), x)', 'RootSum(4*_z**2 + 1, Lambda(_i, _i*log(2*_i + exp(x))))'),
    ('integrate(1/(x*(1 + log(x)**2)), x)', 'RootSum(4*_z**2 + 1, Lambda(_i, _i*log(2*_i + log(x))))'),
    ('integrate(1/(x*(1 + log(x)**2)), x, risch=True)', 'RootSum(4*_z**2 + 1, Lambda(_i, _i*log(2*_i + log(x))))'),
    ('integrate(atan(exp(x)), x)', 'Integral(atan(exp(x)), x)'),
    ('integrate(atan(exp(x)), x, risch=True)', 'Integral(atan(exp(x)), x)'),
    ('integrate(1/(x**2 + 1), x, risch=True)', 'atan(x)'),
    ('integrate(x/(x**4 + 1), x, risch=True)', 'atan(x**2)/2'),
    ('integrate(1/(x**3 + 1), x, risch=True)',
     'log(x + 1)/3 - log(x**2 - x + 1)/6 + sqrt(3)*atan(2*sqrt(3)*x/3 - sqrt(3)/3)/3'),
]


@pytest.mark.parametrize('call, expected', ROOTSUM_CASES)
def test_risch_root_sums(call, expected):
    assert str(eval(call)) == expected
