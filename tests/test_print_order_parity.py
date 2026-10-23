"""Add term and Mul factor order of str, latex and pretty: SymPy 1.14
``Add.as_ordered_terms`` / ``Mul.as_ordered_factors`` (expected strings
were produced by the reference SymPy)."""
import pytest

import symox
from symox import *  # noqa: F401,F403
from symox.physics.quantum import Commutator, Operator
from symox.physics.qho_1d import psi_n

# The star import brings symox.test, which pytest would collect.
del test  # noqa: F821

x, y, z, t, a, b, n, m, omega = symbols('x y z t a b n m omega')
f = Function('f')
F = Function('F')
A, B, C = Operator('A'), Operator('B'), Operator('C')

CASES = [
    ('str', 'series(gamma(x), x, 0, 3)', '1/x - EulerGamma + x*(EulerGamma**2/2 + pi**2/12) + x**2*(-EulerGamma*pi**2/12 - zeta(3)/3 - EulerGamma**3/6) + O(x**3)'),
    ('str', 'series(gamma(x), x, 0, 5)', '1/x - EulerGamma + x*(EulerGamma**2/2 + pi**2/12) + x**2*(-EulerGamma*pi**2/12 - zeta(3)/3 - EulerGamma**3/6) + x**3*(EulerGamma**4/24 + EulerGamma**2*pi**2/24 + EulerGamma*zeta(3)/3 + pi**4/160) + x**4*(-EulerGamma*pi**4/160 - pi**2*zeta(3)/36 - zeta(5)/5 - EulerGamma**2*zeta(3)/6 - EulerGamma**3*pi**2/72 - EulerGamma**5/120) + O(x**5)'),
    ('str', 'series(exp(x)*sin(y + x), x, 0, 4)', 'sin(y) + x*(sin(y) + cos(y)) + x**2*cos(y) + x**3*(-sin(y)/3 + cos(y)/3) + O(x**4)'),
    ('str', '-EulerGamma*pi**2/12 - zeta(3)/3 - EulerGamma**3/6', '-EulerGamma*pi**2/12 - zeta(3)/3 - EulerGamma**3/6'),
    ('str', 'singularityintegrate(SingularityFunction(x, 1, -1)*f(x), x)', 'f(1)*SingularityFunction(x, 1, 0)'),
    ('str', 'f(1)*SingularityFunction(x, 1, 0)', 'f(1)*SingularityFunction(x, 1, 0)'),
    ('str', 'oo - x', '-x + oo'),
    ('str', '2 - 3*x', '2 - 3*x'),
    ('str', 'pi - x', 'pi - x'),
    ('str', 'I - 2*x', '-2*x + I'),
    ('str', 'sin(x)*Max(x, y)*Min(x, y)', 'sin(x)*Max(x, y)*Min(x, y)'),
    ('str', '5*n*sin(a) + n*cosh(a) + zoo*z', '5*n*sin(a) + n*cosh(a) + zoo*z'),
    ('str', 'b**2 + zoo*t + z**x + sinh(b)', 'b**2 + zoo*t + z**x + sinh(b)'),
    ('str', 'x**2 + O(x**3) + 1/x + 3', '1/x + 3 + x**2 + O(x**3)'),
    ('str', 'exp(x) + log(x) + sin(x) + cos(x) + tan(x) + sinh(x) + re(x) + im(x) + F(x)', 'F(x) + exp(x) + log(x) + sin(x) + cos(x) + tan(x) + sinh(x) + re(x) + im(x)'),
    ('str', 'sqrt(3)*Catalan*b + EulerGamma*pi + x*y**2 + x**2*y', 'sqrt(3)*Catalan*b + x**2*y + x*y**2 + EulerGamma*pi'),
    ('str', 'Piecewise((x, x > 0), (0, True)) + Integral(f(x), x) + Derivative(f(x), x)', 'Piecewise((x, x > 0), (0, True)) + Derivative(f(x), x) + Integral(f(x), x)'),
    ('str', '(x + 1)*(y - 2)*sin(x)*exp(y)*Integral(f(y), y)', '(x + 1)*(y - 2)*exp(y)*sin(x)*Integral(f(y), y)'),
    ('terms', 'x**3 + x*y + 2*y**2 + zoo*x + 7', '[x**3, x*y, zoo*x, 2*y**2, 7]'),
    ('factors', '3*x*sin(y)*Max(x, y)*exp(z)*f(x)', '[3, x, f(x), exp(z), sin(y), Max(x, y)]'),
    ('latex', 'series(exp(x)*sin(y + x), x, 0, 4).removeO()', 'x^{3} \\left(- \\frac{\\sin{\\left(y \\right)}}{3} + \\frac{\\cos{\\left(y \\right)}}{3}\\right) + x^{2} \\cos{\\left(y \\right)} + x \\left(\\sin{\\left(y \\right)} + \\cos{\\left(y \\right)}\\right) + \\sin{\\left(y \\right)}'),
    ('str', '((A + B)**3).expand()', 'A*B*A + A*B**2 + A**2*B + A**3 + B*A*B + B*A**2 + B**2*A + B**3'),
    ('str', 'Commutator(A*B, C).expand(commutator=True)', '[A,C]*B + A*[B,C]'),
    ('str', '2*A*B - 3*x*B*A + 1', '-3*x*B*A + 1 + 2*A*B'),
    ('str', 'psi_n(0, x, m, omega)', '(m*omega)**(1/4)*exp(-m*omega*x**2/(2*hbar))/(hbar**(1/4)*pi**(1/4))'),
    ('latex', 'sin(x)*Max(x, y)*Min(x, y) + oo - x', '- x + \\sin{\\left(x \\right)} \\max\\left(x, y\\right) \\min\\left(x, y\\right) + \\infty'),
    ('pretty', '5*n*sin(a) + n*cosh(a) + zoo*z', '5*n*sin(a) + n*cosh(a) + zoo*z'),
]


def _render(kind, value):
    if kind == 'str':
        return str(value)
    if kind == 'terms':
        return str(value.as_ordered_terms())
    if kind == 'factors':
        return str(value.as_ordered_factors())
    if kind == 'latex':
        return symox.latex(value)
    return symox.pretty(value, use_unicode=False)


@pytest.mark.parametrize('kind,source,expected', CASES, ids=[c[1] for c in CASES])
def test_printing_order_matches_reference(kind, source, expected):
    assert _render(kind, eval(source)) == expected


def test_order_terms_follow_the_regular_terms():
    # Order terms sort among themselves, after the reversed regular terms.
    assert str(x**2 + O(x**3) + 1/x + 3) == '1/x + 3 + x**2 + O(x**3)'


def test_series_zeta_coefficient_is_the_zeta_function():
    # The x**2 coefficient of gamma(x) holds zeta(3), not an undefined zeta.
    c = series(gamma(x), x, 0, 3).removeO().coeff(x, 2)
    assert c == -EulerGamma*pi**2/12 - zeta(3)/3 - EulerGamma**3/6


def test_hbar_prints_as_the_reference_constant():
    from symox.physics.quantum.constants import hbar
    assert symox.latex(hbar) == r'\hbar'
    assert symox.pretty(hbar) == '\N{PLANCK CONSTANT OVER TWO PI}'
    assert symox.pretty(hbar, use_unicode=False) == 'hbar'
    assert symox.latex(psi_n(0, x, m, omega)) == (
        r'\frac{\sqrt[4]{m \omega} e^{- \frac{m \omega x^{2}}{2 \hbar}}}'
        r'{\sqrt[4]{\hbar} \sqrt[4]{\pi}}')
