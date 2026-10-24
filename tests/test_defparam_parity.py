"""defparam (wave 17): distribution integrals and series as in SymPy 1.14.

Expected strings are SymPy 1.14 outputs. Held cases are those where the
SymPy value is wrong for some parameter values (documented in
docs/algorithmic-parity.md, section defparam).
"""
import mpmath
import pytest

from symox import *

del test

x = Symbol('x')
f = Function('f')
a = Symbol('a')
ar = Symbol('ar', real=True)
ap = Symbol('ap', positive=True)
an = Symbol('an', negative=True)

SAME = [
    ('integrate(DiracDelta(x**2 - ap**2), x)', 'Heaviside(-ap + x)/(2*ap) + Heaviside(ap + x)/(2*ap)'),
    ('integrate(DiracDelta(x**2 - ap), x)', 'Heaviside(-sqrt(ap) + x)/(2*sqrt(ap)) + Heaviside(sqrt(ap) + x)/(2*sqrt(ap))'),
    ('integrate(DiracDelta(x**2 - ap), (x, -oo, oo))', '1/sqrt(ap)'),
    ('integrate(DiracDelta(x**2 - ap**2)*exp(x), (x, -oo, oo))', 'exp(ap)/(2*ap) + exp(-ap)/(2*ap)'),
    ('integrate(DiracDelta(x**2 - ap**2), (x, 0, 1))', 'Heaviside(1 - ap)/(2*ap)'),
    ('integrate(DiracDelta(x**2 - an**2), (x, -oo, oo))', '-1/an'),
    ('integrate(DiracDelta((x - ap)*(x - 2*ap)), (x, -oo, oo))', '2/ap'),
    ('integrate(DiracDelta(x**2 - 3*ap*x + 2*ap**2), x)', 'Heaviside(-2*ap + x)/ap + Heaviside(-ap + x)/ap'),
    ('integrate(DiracDelta(x, 1)*f(x), x)', '-Heaviside(x)*Subs(Derivative(f(x), x), x, 0)'),
    ('integrate(DiracDelta(x - 1, 1)*f(x), x)', '-Heaviside(x - 1)*Subs(Derivative(f(x), x), x, 1)'),
    ('integrate(DiracDelta(x - 1, 1)*f(x), (x, -oo, oo))', '-Subs(Derivative(f(x), x), x, 1)'),
    ('integrate(DiracDelta(x - 1, 1)*f(x), (x, 0, 2))', '-Subs(Derivative(f(x), x), x, 1)'),
    ('integrate(DiracDelta(x, 2)*f(x), (x, -oo, oo))', 'Subs(Derivative(f(x), (x, 2)), x, 0)'),
    ('integrate(DiracDelta(x, 2)*f(x), x)', 'Heaviside(x)*Subs(Derivative(f(x), (x, 2)), x, 0)'),
    ('integrate(DiracDelta(x - ar, 1)*f(x), (x, -oo, oo))', '-Derivative(f(ar), ar)'),
    ('integrate(DiracDelta(-2*x + 1), x)', 'Heaviside(x - 1/2)/2'),
    ('integrate(DiracDelta(-2*x + 1)*x, x)', 'Heaviside(x - 1/2)/4'),
    ('integrate(DiracDelta(3*x - 1)*exp(x), x)', 'exp(1/3)*Heaviside(x - 1/3)/3'),
    ('integrate(DiracDelta(3*x - 1)*exp(x), (x, 0, 1))', 'exp(1/3)/3'),
    ('integrate(DiracDelta(ap*x - 1), x)', 'Heaviside(x - 1/ap)/ap'),
    ('integrate(DiracDelta(an*x - 1), x)', '-Heaviside(x - 1/an)/an'),
    ('integrate(DiracDelta(2*x), x)', 'Heaviside(x)/2'),
    ('integrate(DiracDelta(x/2 + 1)*x, x)', '-4*Heaviside(x + 2)'),
    ('series(Heaviside(x - 1), x)', 'O(x**6)'),
    ('series(Heaviside(x - 1)*exp(x), x)', 'O(x**6)'),
    ('series(DiracDelta(x - 1), x)', 'O(x**6)'),
    ('series(Heaviside(x + 1), x)', '1'),
    ('series(Heaviside(x - 1)/x, x)', 'O(x**6)'),
    ('series(Heaviside(x - 1) + x, x)', 'x + O(x**6)'),
    ('series(Heaviside(x-1)*sin(x), x, 0, 3)', 'O(x**3)'),
    ('integrate(Heaviside(x**2 - 1), (x, -2, 2)).n()', '2.00000000000000'),
    ('integrate(Heaviside(x**2 - 4), (x, 0, 3))', '-meijerg(((1/2, 1), ()), ((), (-1/2, 0)), 1) + 3*meijerg(((1/2, 1), ()), ((), (-1/2, 0)), 9/4)/2'),
    ('integrate(Heaviside(x**2 - 1)*exp(-x), (x, 0, oo))', '(-2*sqrt(pi)*sinh(1) + 2*sqrt(pi)*cosh(1))/(2*sqrt(pi))'),
    ('integrate(Heaviside(1 - x**2), x)', 'Piecewise((x, Abs(x**2) < 1), (x*meijerg(((1/2,), (1,)), ((0,), (-1/2,)), x**2)/2, True))'),
    ('integrate(Heaviside(x**2 - 1)*x, (x, 0, 2))', '-meijerg(((0, 1), ()), ((), (-1, 0)), 1)/2 + 2*meijerg(((0, 1), ()), ((), (-1, 0)), 4)'),
    ('series(Heaviside(x)/x, x)', '1/x'),
    ('series(Heaviside(x), x)', '1'),
    ('series(Heaviside(-x), x)', '0'),
    ('series(Heaviside(1 - x), x)', '1 + O(x**6)'),
    ('series(Heaviside(-x - 1), x)', '0'),
    ('series(DiracDelta(x + 1), x)', '0'),
    ('series(DiracDelta(x), x)', '0'),
    ('series(DiracDelta(x, 1), x)', '0'),
    ('series(DiracDelta(x - 1, 1), x)', 'O(x**6)'),
    ('series(Heaviside(x - 1)*x, x)', 'O(x**6)'),
    ('series(Heaviside(x**2 - 1), x)', 'O(x**6)'),
    ('series(Heaviside(-x**2), x)', '0'),
    ('series(Heaviside(x - 1) + 1, x)', '1 + O(x**6)'),
    ('series(Heaviside(-x)/x, x)', '0'),
    ('series(Heaviside(x - 1)/x**2 + 1/x, x)', '1/x + O(x**6)'),
    ('series(Heaviside(sin(x) - 1), x)', 'O(x**6)'),
    ('series(Heaviside(x + 1) + x**7, x)', '1 + O(x**6)'),
    ('series(Heaviside(x)*sin(x), x)', 'x - x**3/6 + x**5/120 + O(x**6)'),
    ('series(Heaviside(x - 2), x, 1)', 'O((x - 1)**6, (x, 1))'),
    ('series(Heaviside(x - 1), x, 2)', '1'),
    ('series(Heaviside(x - 1), x, 1)', '1'),
    ('series(Heaviside(x - 1)*x, x, 1, 3)', 'x'),
    ('series(Heaviside(x)*exp(x), x, 0, 4)', '1 + x + x**2/2 + x**3/6 + O(x**4)'),
    ('series(Heaviside(x + 1)*exp(x), x, 0, 4)', '1 + x + x**2/2 + x**3/6 + O(x**4)'),
    ('series(Heaviside(x - 1)*cos(x), x, 0, 4)', 'O(x**4)'),
    ("series(Heaviside(x - 1), x, 0, 4, dir='-')", '0'),
    ("series(Heaviside(x), x, 0, 4, dir='-')", '0'),
    ("series(Heaviside(-x), x, 0, 4, dir='-')", '1'),
]


@pytest.mark.parametrize('statement,expected', SAME)
def test_same_as_sympy(statement, expected):
    assert str(eval(statement)) == expected


def test_even_heaviside_on_symmetric_interval():
    # SymPy 1.14: 2 + meijerg terms that evaluate to 2 (same Meijer G form)
    r = integrate(Heaviside(1 - x**2), (x, -2, 2))
    assert r.has(meijerg)
    with mpmath.workdps(30):
        assert abs(complex(N(r, 30)) - 2) < 1e-20
    for e, lims, value in [
        (Heaviside(x**2 - 1), (x, -3, 1), 2),
        (Heaviside(1 - x**2), (x, -3, Rational(1, 2)), Rational(3, 2)),
        (x*Heaviside(x**2 - 1), (x, -3, 2), Rational(-5, 2)),
        (x**2*Heaviside(1 - x**2), (x, -3, 2), Rational(2, 3)),
    ]:
        assert abs(complex(N(integrate(e, lims))) - float(value)) < 1e-12


def test_even_heaviside_primitive_valid_for_negative_x():
    for e in [Heaviside(x**2 - 1), Heaviside(1 - x**2), x*Heaviside(x**2 - 1)]:
        F = integrate(e, x)
        for p in [-3, Rational(-1, 2), Rational(1, 2), 3]:
            assert abs(complex(N(diff(F, x).subs(x, p))) - float(e.subs(x, p))) < 1e-12


HELD = [
    # SymPy: /(2*Abs(a)), infinite at a = 0 where DiracDelta(x**2) has no simple root
    lambda: integrate(DiracDelta(x**2 - a**2), x),
    lambda: integrate(DiracDelta(x**2 - ar**2), x),
    # SymPy: -Subs(Derivative(f(x), x), x, 0), missing the 1/(2*|2|) scale (value -f'(0)/4)
    lambda: integrate(DiracDelta(2*x, 1)*f(x), (x, -oo, oo)),
    # SymPy: Heaviside(x - 1/a)/(a*Abs(a)), infinite at a = 0 where the integrand vanishes
    lambda: integrate(DiracDelta(a*x - 1)*x, x),
    # SymPy raises TypeError (divergent integral)
    lambda: integrate(Heaviside(x**2 - 1), (x, 0, oo)),
    # SymPy: a value although the integrand has a pole on the interval for some a
    lambda: integrate(1/(x + a), (x, 0, 1)),
    lambda: integrate(1/(1 + a*x), (x, 0, 1)),
    lambda: integrate(1/x, (x, a, 1)),
    lambda: integrate(1/(x**2 + a), (x, 0, oo)),
    lambda: integrate(1/(x**2 + 2*a*x + 1), (x, -oo, oo)),
]


@pytest.mark.parametrize('k', range(len(HELD)))
def test_held_where_sympy_is_wrong(k):
    assert HELD[k]().has(Integral)


def test_real_slope_delta_is_piecewise():
    # SymPy: Heaviside(x - 1/ar)/Abs(ar), infinite at ar = 0 where the integrand vanishes
    r = integrate(DiracDelta(ar*x - 1), x)
    assert r == Piecewise((Heaviside(ar*x - 1)/ar, Ne(ar, 0)), (0, True))


def test_series_exact_through_jump():
    # SymPy runs its Taylor expansion through the jump (documented correction)
    assert str(series(Heaviside(sin(x)), x, 0, 3)) == '1'
