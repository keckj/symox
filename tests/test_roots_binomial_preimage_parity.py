"""Preimages of decomposed roots through a real binomial keep SymPy's form."""
from symox import Function, I, Symbol, cos, dsolve, roots, sin, sqrt, Eq

x = Symbol('x')


def test_negative_real_base_gives_imaginary_roots():
    assert roots(x**4 + 4*x**2 + 1, x) == {
        I*sqrt(sqrt(3) + 2): 1, -I*sqrt(sqrt(3) + 2): 1,
        I*sqrt(2 - sqrt(3)): 1, -I*sqrt(2 - sqrt(3)): 1}


def test_positive_real_base_gives_real_roots():
    assert roots(x**4 - 4*x**2 + 1, x) == {
        sqrt(sqrt(3) + 2): 1, -sqrt(sqrt(3) + 2): 1,
        sqrt(2 - sqrt(3)): 1, -sqrt(2 - sqrt(3)): 1}


def test_constant_coefficient_ode_uses_trigonometric_basis():
    f = Function('f')
    C1, C2, C3, C4 = (Symbol('C%d' % i) for i in range(1, 5))
    a, b = sqrt(2 - sqrt(3)), sqrt(sqrt(3) + 2)
    assert dsolve(f(x).diff(x, 4) + 4*f(x).diff(x, 2) + f(x)) == Eq(
        f(x), C1*sin(a*x) + C2*sin(b*x) + C3*cos(a*x) + C4*cos(b*x))
