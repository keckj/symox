"""Substitution of a power of a symbol into the bare symbol.

SymPy's Symbol._eval_subs treats a symbol y as y**1 for an old power of
y, so (y).subs(1/y, u) is 1/u; a differentiation variable changed this
way is held in Subs with a xi_i dummy. The expected strings were produced
by the reference SymPy 1.14 (manualintegrate relies on this rule to find
the substitution u = 1/y).
"""
import pytest
from symox import *

# the star import brings the library test runner, which is not a test here
del test

y, z = symbols('y z')
u = Symbol('u')
p = Symbol('p', positive=True)
n = Symbol('n', nonnegative=True)
f = Function('f')

CASES = [
    ('(y, 1/y, u)', '1/u'),
    ('(-y/(y-1), 1/y, u)', '-1/(u*(-1 + 1/u))'),
    ('(y - 1, 1/y, u)', '-1 + 1/u'),
    ('(y**3/(y+1), 1/y, u)', '1/(u**3*(1 + 1/u))'),
    ('(2*y, 1/y, u)', '2/u'),
    ('(y*z, 1/y, u)', 'z/u'),
    ('(sqrt(y), 1/y, u)', 'sqrt(1/u)'),
    ('(p, p**-2, u)', '1/sqrt(u)'),
    ('(p, p**2, u)', 'sqrt(u)'),
    ('(p, sqrt(p), u)', 'u**2'),
    ('(n, n**2, u)', 'sqrt(u)'),
    ('(n, sqrt(n), u)', 'u**2'),
    ('(y, sqrt(y), u)', 'u**2'),
    ('(y, y**Rational(1, 3), u)', 'u**3'),
    ('(y, y**-1.0, u)', 'u**(-1.0)'),
    ('(y, y**2.0, u)', 'y'),
    ('(y, y**z, u)', 'y'),
    ('(y, y**-2, u)', 'y'),
    ('(y**z, y**-z, u)', '1/u'),
    ('(sin(y) + y**2, 1/y, u)', 'sin(1/u) + u**(-2)'),
    ('(f(y), 1/y, u)', 'f(1/u)'),
    ('(y*exp(y), 1/y, u)', 'exp(1/u)/u'),
    ('(y, exp(-y), u)', 'y'),
    ('(y, (2*y)**-1, u)', 'y'),
    ('(log(y), 1/y, u)', 'log(1/u)'),
    ('(y*z, y**-2, u)', 'y*z'),
    ('(y - 1, 1/(y - 1), u)', 'y - 1'),
    ('(Integral(y, z), 1/y, u)', 'Integral(1/u, z)'),
    ('(Derivative(f(y), y), 1/y, u)', 'Subs(Derivative(f(_xi_0), _xi_0), _xi_0, 1/u)'),
    ('(Derivative(f(y), (y, 2)), 1/y, u)', 'Subs(Derivative(f(_xi_0), (_xi_0, 2)), _xi_0, 1/u)'),
    ('(Derivative(y*f(y), y), 1/y, u)', 'Subs(Derivative(_xi_0*f(_xi_0), _xi_0), _xi_0, 1/u)'),
    ('(Derivative(f(y), y), 1/y, 1/z)', 'Derivative(f(z), z)'),
    ('(Derivative(f(y, z), y), 1/y, 1/z)', 'Subs(Derivative(f(y, z), y), 1/y, 1/z)'),
    ('(y**Rational(3, 2), 1/y, u)', '(1/u)**(3/2)'),
    ('(y**-3, 1/y, u)', 'u**3'),
]


@pytest.mark.parametrize('case, expected', CASES)
def test_symbol_power_substitution_matches_the_reference(case, expected):
    expr, old, new = eval(case)
    assert str(expr.subs(old, new)) == expected


def test_manual_substitution_finds_the_reciprocal_variable():
    from symox.integrals.manualintegrate import find_substitutions
    v = Dummy('u')
    assert [(str(a), str(b), str(c.subs(v, u))) for a, b, c in
            find_substitutions(1/(y*(y - 1)), y, v)] == [('1/y', '1', '1/(u - 1)')]
    expected = 'log(y)*log(-1 + 1/y) - Integral(log(-1 + 1/y)/y, y)'
    assert str(integrate(log(y)/(y**2*(1 - 1/y)), y, manual=True)) == expected
