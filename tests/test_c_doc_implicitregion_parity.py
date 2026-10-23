"""ImplicitRegion regular points and rational parametrizations (SymPy 1.14).

Expected values are hard-coded from the reference.
"""
import pytest

from symox import Eq, Symbol, symbols, sympify
from symox.vector import ImplicitRegion
from symox.vector._ternary_reference import ternary_normal_diophantine

x, y, z, s, t = symbols('x y z s t')
T = Symbol('t', real=True)


def test_rational_parametrization_docstring():
    parabola = ImplicitRegion((x, y), y**2 - 4*x)
    assert str(parabola.rational_parametrization()) == '(4/t**2, 4/t)'
    circle = ImplicitRegion((x, y), Eq(x**2 + y**2, 4))
    assert circle.regular_point() == (0, -2)
    assert str(circle.rational_parametrization()) == '(4*t/(t**2 + 1), 4*t**2/(t**2 + 1) - 2)'
    cubic = ImplicitRegion((x, y), x**3 + x**2 - y**2)
    assert str(cubic.rational_parametrization()) == '(t**2 - 1, t*(t**2 - 1))'
    assert str(cubic.rational_parametrization(parameters=(t))) == '(t**2 - 1, t*(t**2 - 1))'
    sphere = ImplicitRegion((x, y, z), x**2 + y**2 + z**2 - 4)
    assert str(sphere.rational_parametrization(parameters=(t, s))) == \
        '(-2 + 4/(s**2 + t**2 + 1), 4*s/(s**2 + t**2 + 1), 4*t/(s**2 + t**2 + 1))'
    c = ImplicitRegion((x, y), (x - 1/2)**2 + y**2 - (1/4)**2)
    assert str(c.rational_parametrization(reg_point=(3/4, 0))) == \
        '(0.75 - 0.5/(t**2 + 1), -0.5*t/(t**2 + 1))'


@pytest.mark.parametrize('eq, point, param', [
    (x**2 + 2*y**2 - 3, (-1, -1), '((4*t + 2)/(2*t**2 + 1) - 1, t*(4*t + 2)/(2*t**2 + 1) - 1)'),
    (3*x**2 + 5*y**2 - 8, (-1, -1), '((10*t + 6)/(5*t**2 + 3) - 1, t*(10*t + 6)/(5*t**2 + 3) - 1)'),
    (x**2 + y**2 - 2, (-1, -1), '((2*t + 2)/(t**2 + 1) - 1, t*(2*t + 2)/(t**2 + 1) - 1)'),
    (2*x**2 + 3*y**2 - 5, (-1, -1), '((6*t + 4)/(3*t**2 + 2) - 1, t*(6*t + 4)/(3*t**2 + 2) - 1)'),
    (x**2 - y**2 - 3, (2, 1), '(2 + (2*t - 4)/(1 - t**2), t*(2*t - 4)/(1 - t**2) + 1)'),
    (x**2 + y**2 - 25, (0, -5), '(10*t/(t**2 + 1), 10*t**2/(t**2 + 1) - 5)'),
])
def test_conic_regular_point_and_parametrization(eq, point, param):
    region = ImplicitRegion((x, y), eq)
    assert region.regular_point() == point
    result = region.rational_parametrization()
    expected = sympify(param, locals={'t': T})
    assert all((a - b).simplify() == 0 for a, b in zip(result, expected))
    assert str(result) == param


def test_conic_without_rational_point():
    with pytest.raises(ValueError, match='Rational Point on the conic does not exist'):
        ImplicitRegion((x, y), 5*x**2 + 2*x*y + 3*y**2 + x + y - 11).regular_point()


@pytest.mark.parametrize('abc, expected', [
    ((1, 1, -1), ['(2*p*q, p**2 - q**2, p**2 + q**2)']),
    ((1, 2, -3), ['(2*p**2 - 3*q**2, -2*p**2 + 6*p*q - 3*q**2, 2*p**2 - 4*p*q + 3*q**2)']),
    ((2, 3, -5), ['(3*p**2 - 5*q**2, -3*p**2 + 10*p*q - 5*q**2, 3*p**2 - 6*p*q + 5*q**2)']),
    ((1, 1, -2), ['(p**2 - 2*q**2, -p**2 + 4*p*q - 2*q**2, p**2 - 2*p*q + 2*q**2)']),
    ((1, -12, -1), ['(12*p**2 + q**2, -2*p*q, 12*p**2 - q**2)']),
    ((2, -3, -6), ['(3*p**2 + 6*q**2, -2*p**2 - 4*p*q + 4*q**2, p**2 - 4*p*q - 2*q**2)']),
    ((1, 1, 1), ['(0, 0, 0)']),
    ((2, -7, -6), ['(0, 0, 0)']),
    ((3, -2, 6), ['(0, 0, 0)']),
    # the reference keeps an uncancelled parameter r here
    ((2, -7, 6), ['(-7*p**2 + 6*q**2, 14*p**2 - 12*p*q + 20*p*r + 12*q**2, '
                  '-7*p**2 + 28*p*q - 6*q**2 + 20*q*r)']),
])
def test_reference_ternary_normal_diophantine(abc, expected):
    # reference: sorted(map(str, diophantine(a*x**2 + b*y**2 + c*z**2)))
    assert sorted(map(str, ternary_normal_diophantine(*abc))) == expected


def test_reference_ternary_requires_positive_leading_coefficient():
    with pytest.raises(ValueError):
        ternary_normal_diophantine(-1, 1, 1)
