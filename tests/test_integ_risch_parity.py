"""Parity of distribution integration and of the special integrands that
SymPy 1.14's Integral handles before integrating (Eq, Sum, None limits,
Order). Expected values are hard-coded from SymPy 1.14."""

import pytest

from symox import (
    Abs, CRootOf, DiracDelta, Eq, Function, Heaviside, Integral, O, Rational,
    S, Sum, cos, diff, exp, integrate, oo, sign, sin, sqrt, summation, symbols,
    Symbol, Piecewise, Ne, Poly,
)

x, y, a = symbols('x y a')
f = Function('f')


def _check_antiderivative(F, e):
    # Away from the support points the derivative must agree pointwise with
    # the integrand (zero for a pure delta), and the jumps must carry the
    # delta weights; both are tested by the definite integrals below.
    assert not F.has(Integral)


@pytest.mark.parametrize('expr, expected', [
    (DiracDelta(x**2 - 1), Heaviside(x - 1)/2 + Heaviside(x + 1)/2),
    (DiracDelta(x**3 - 2*x),
     Heaviside(x)/2 + Heaviside(x - sqrt(2))/4 + Heaviside(x + sqrt(2))/4),
    (cos(x)*DiracDelta(x**2 - 1),
     cos(1)*Heaviside(x - 1)/2 + cos(1)*Heaviside(x + 1)/2),
    (DiracDelta(x**2 - 1)*Heaviside(x), Heaviside(x - 1)/2),
    (DiracDelta(2*x**2 - 8)/x, Heaviside(x - 2)/16 - Heaviside(x + 2)/16),
    (DiracDelta(4 - x**2), Heaviside(x - 2)/4 + Heaviside(x + 2)/4),
    (DiracDelta(x**2 - 2),
     sqrt(2)*Heaviside(x - sqrt(2))/4 + sqrt(2)*Heaviside(x + sqrt(2))/4),
    (DiracDelta(x - 1)*Heaviside(x), Heaviside(x - 1)),
    (DiracDelta(x - 1)*Heaviside(x - 2), S.Zero),
    (DiracDelta(x - 1)*Heaviside(x - a), Heaviside(1 - a)*Heaviside(x - 1)),
    (DiracDelta(x - 1)*sign(x), Heaviside(x - 1)),
    (Abs(x)*DiracDelta(x + 2), 2*Heaviside(x + 2)),
    (f(x)*DiracDelta(x), f(0)*Heaviside(x)),
    (f(x)*DiracDelta(x - 1), f(1)*Heaviside(x - 1)),
])
def test_delta_indefinite(expr, expected):
    F = integrate(expr, x)
    _check_antiderivative(F, expr)
    assert F == expected


@pytest.mark.parametrize('expr, lims, expected', [
    (DiracDelta(x**2 - 1), (x, -oo, oo), S.One),
    (x**2*DiracDelta(x**2 - 4), (x, 0, oo), S.One),
    (DiracDelta(x**2 - 3*x + 2)*x, (x, 0, 3), S(3)),
    (f(x)*DiracDelta(x), (x, -oo, oo), f(0)),
    (cos(x)*DiracDelta(x**2 - 1), (x, -oo, oo), cos(1)),
])
def test_delta_definite(expr, lims, expected):
    assert integrate(expr, lims) == expected


def test_delta_nonaffine_definite_numeric():
    # sum over the simple roots r of g of h(r)/|g'(r)|
    v = integrate(exp(x)*DiracDelta(x**2 - 2), (x, -oo, oo))
    ref = (exp(sqrt(2)) + exp(-sqrt(2)))/(2*sqrt(2))
    assert abs((v - ref).evalf(30)) < S(10)**-25


def test_delta_crootof_root_kept():
    # Intentional: SymPy returns 0 for DiracDelta(x**5 - x - 1); the real
    # root CRootOf(x**5 - x - 1, 0) carries weight 1/|g'(r)|.
    r = CRootOf(x**5 - x - 1, 0)
    v = integrate(DiracDelta(x**5 - x - 1)*exp(x), (x, -oo, oo))
    ref = exp(r)/Abs(5*r**4 - 1)
    assert abs((v - ref).evalf(30)) < S(10)**-25
    assert integrate(DiracDelta(x**5 - x - 1), (x, -oo, oo)) != 0


def test_delta_held_cases():
    # A pole at a support point (SymPy: zoo) and a squared delta stay held.
    assert integrate(DiracDelta(x**2 - 1)/(x - 1), x).has(Integral)
    assert integrate(DiracDelta(x)**2, x).has(Integral)


def test_delta_derivative_definite_values():
    # Exact distributional values of derivatives of deltas.
    assert integrate(x*DiracDelta(x - 1, 1), (x, -oo, oo)) == -1
    assert integrate(x*DiracDelta(x - 1, 1)*Heaviside(x), (x, -oo, oo)) == -1
    # Intentional: SymPy gives -exp(3/2) (missing the 1/(2*|2|) scale).
    assert integrate(exp(x)*DiracDelta(2*x - 3, 1), (x, -oo, oo)) == -exp(Rational(3, 2))/4


def test_none_limits():
    assert integrate(x**2, (x, None, 1)) == Rational(1, 3)
    assert integrate(x**2, (x, 1, None)) == Rational(-1, 3)
    assert integrate('x**2', ('x', '1', None)) == Rational(-1, 3)


def test_equality_integrand():
    assert integrate(Eq(x, y), (x, 0, 1)) == Eq(S.Half, y)
    assert integrate(Eq(x, y), x) == Eq(x**2/2, x*y)


def test_order_integrand():
    assert integrate(O(x**5), x) == O(x**6)
    assert integrate(1 + x + O(x**3), x) == x + x**2/2 + O(x**4)


def test_sum_integrand():
    n = Symbol('n', integer=True)
    N = Symbol('N')
    assert integrate(Sum(x*y, (y, 1, 3)), (x, 0, 1)) == 3
    assert integrate(summation(sin(n*x), (n, 1, N)), x) == \
        Sum(Piecewise((-cos(n*x)/n, Ne(n, 0)), (0, True)), (n, 1, N))
    held = integrate(Sum(x**n, (n, 1, oo)), x)
    assert held == Integral(Sum(x**n, (n, 1, oo)), x)
    xx, nn = symbols('x n')
    for s, v in [(Sum(xx, (xx, y, nn)), y), (Sum(xx, (xx, 1, nn)), nn),
                 (Sum(xx, (xx, 1, y)), xx)]:
        with pytest.raises(ValueError, match='Limit of the sum cannot be an integration variable'):
            integrate(s, v)
