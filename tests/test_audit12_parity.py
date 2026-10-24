"""Wave 12 audit (lots c1, c2, c4): power_base hint of expand.

Expected values are hard-coded from SymPy 1.14.
"""
import pytest

from symox import (
    Mul, Rational, Symbol, exp, expand, expand_power_base, log, symbols,
)


x, y, z, a, n = symbols('x y z a n')
p = Symbol('p', positive=True)
xp = Symbol('x', positive=True)
k = Symbol('k', integer=True)


@pytest.mark.parametrize('value, expected', [
    # expand docstring: power_base splits the nonnegative factors
    (lambda: ((2*y)**z).expand(power_base=True), '2**z*y**z'),
    (lambda: expand((2*y)**z), '2**z*y**z'),
    (lambda: expand((x*y)**z), '(x*y)**z'),
    (lambda: expand((p*y)**z), 'p**z*y**z'),
    (lambda: expand((-2*y)**z), '2**z*(-y)**z'),
    (lambda: expand((a*y)**k), 'a**k*y**k'),
    (lambda: expand((x*(y + z))**x, power_base=False), '(x*y + x*z)**x'),
    (lambda: ((x*(y + z))**x).expand(mul=False), '(x*(y + z))**x'),
    (lambda: expand((2*(y + z))**x), '(2*y + 2*z)**x'),
    (lambda: expand(((y + z)*x)**y), '(x*y + x*z)**y'),
    # factors of a product are visited as power nodes
    (lambda: expand_power_base(a*(a*xp)**n), 'a*a**n*x**n'),
    (lambda: expand(a*(a*xp)**n), 'a*a**n*x**n'),
    (lambda: expand((a*xp)**n*(a + 1)), 'a*a**n*x**n + a**n*x**n'),
    (lambda: expand((a*xp)**n*(a + 1), mul=False), 'a**n*x**n*(a + 1)'),
    (lambda: expand((2*xp)**n*(xp + 1)**2),
     '2**n*x**2*x**n + 2*2**n*x*x**n + 2**n*x**n'),
    (lambda: expand(exp(n*log(2*a))*(2*a)**n),
     '2**n*a**n*exp(n*log(2))*exp(n*log(a))'),
    (lambda: expand((3*a)**(n + 1)), '3*3**n*a**(n + 1)'),
    (lambda: expand((xp*(a + 1))**n*(a + 1)), 'a*(a*x + x)**n + (a*x + x)**n'),
])
def test_power_base_hint(value, expected):
    assert str(value()) == expected


def test_mul_before_power_base():
    # The reference applies mul before power_base (sorted hints).
    xs, ys, zs = symbols('x y z', positive=True)
    assert str(expand((xs*(ys + zs))**xs)) == '(x*y + x*z)**x'
    assert str(expand((xs*(ys + zs))**xs, mul=False)) == 'x**x*(y + z)**x'
    # (x*(y + z))**x is not x**x*(y + z)**x for x = -1/2, y + z = -1.
    e = expand((x*(y + z))**x)
    assert str(e) == '(x*y + x*z)**x'
    point = {x: Rational(-1, 2), y: -1, z: 0}
    assert abs(complex(e.subs(point).evalf()) - complex(((x*(y + z))**x).subs(point).evalf())) < 1e-12


def test_unsplit_products_are_kept():
    assert str(Mul(a, Mul(a**n, xp**n))) == 'a*a**n*x**n'
    assert str(expand((x*(y + z))**2)) == 'x**2*y**2 + 2*x**2*y*z + x**2*z**2'
