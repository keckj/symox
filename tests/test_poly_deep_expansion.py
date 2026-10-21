"""Deep expansion of an expression before it becomes a polynomial.

The reference expands the arguments of functions and splits exp of a sum
when it builds a polynomial, so exp(x*(x + 1)) has the generators exp(x)
and exp(x**2). cancel, factor and the Risch differential extension rely on
it: integrate(exp(-x**2/(2*(2 - I)**2)), x, risch=True) stays unevaluated
(the extension raises NotImplementedError) instead of failing in a
resultant. The expected strings were produced by SymPy 1.14.
"""
import pytest
from symox import *

# the star import brings the library test runner, which is not a test here
del test

x, y = symbols('x y')

CASES = [
    ('Poly(exp(-x**2/(2*(2 - I)**2)))', "Poly((exp(-I*x**2/25))**2*(exp(-x**2/50))**3, exp(-I*x**2/25), exp(-x**2/50), domain='ZZ')"),
    ('Poly(exp(x*(x + 1)))', "Poly((exp(x))*(exp(x**2)), exp(x), exp(x**2), domain='ZZ')"),
    ('Poly(sin(x*(x + 1)) + 1)', "Poly((sin(x**2 + x)) + 1, sin(x**2 + x), domain='ZZ')"),
    ('Poly(exp(x*(x + 1)), x)', 'EXC PolynomialError exp(x) contains an element of the set of generators.'),
    ('Poly(cos(exp(x*(y + 1))))', "Poly((cos(exp(x)*exp(x*y))), cos(exp(x)*exp(x*y)), domain='ZZ')"),
    ('Poly(x**(y*(y + 1)), x)', 'EXC PolynomialError x**(y**2 + y) contains an element of the set of generators.'),
    ('Poly(sin(x*(x + 1))**2 + x)', "Poly(x + (sin(x**2 + x))**2, x, sin(x**2 + x), domain='ZZ')"),
    ('Poly(exp(2*(x + y)))', "Poly((exp(x))**2*(exp(y))**2, exp(x), exp(y), domain='ZZ')"),
    ('Poly(log(x*(y + 1)), x)', 'EXC PolynomialError log(x*y + x) contains an element of the set of generators.'),
    ('cancel(exp(-x**2/(2*(2 - I)**2)))', 'exp(-3*x**2/50)*exp(-2*I*x**2/25)'),
    ('cancel(sin(x)/exp(x*(x + 1)))', 'exp(-x)*exp(-x**2)*sin(x)'),
    ('cancel((exp(x*(x + 1)) - exp(x**2 + x))/x)', '0'),
    ('factor(exp(-x**2/(2*(2 - I)**2)))', 'exp(-3*x**2/50)*exp(-2*I*x**2/25)'),
    ('factor(Integral(exp(-x**2/(2*(2 - I)**2)), x))', 'Integral(exp(-3*x**2/50)*exp(-2*I*x**2/25), x)'),
    ('together(1/exp(x*(x + 1)) + 1)', '(exp(x*(x + 1)) + 1)*exp(-x*(x + 1))'),
    ('gcd(exp(x*(x + 1)), exp(x)*exp(x**2))', 'exp(x)*exp(x**2)'),
    ('integrate(exp(-x**2/(2*(2 - I)**2)), x, risch=True)', 'Integral(exp(-3*x**2/50)*exp(-2*I*x**2/25), x)'),
    ('integrate(exp(-x**2/(2*(2 - I)**2)), x, heurisch=True)', 'Integral(exp(-3*x**2/50)*exp(-2*I*x**2/25), x)'),
]


@pytest.mark.parametrize('case, expected', CASES)
def test_deep_expansion_matches_the_reference(case, expected):
    try:
        result = str(eval(case))
    except Exception as error:
        result = 'EXC %s %s' % (type(error).__name__, error)
    assert result == expected


def test_heurisch_complex_erf_issue_26338():
    from symox.integrals.heurisch import heurisch
    r = Symbol('r', real=True)
    assert heurisch(exp(-x**2/(2*(2 - I)**2)), x, hints=[]) is None
    assert heurisch(exp(-r**2/(2*(2 - I)**2)), r, hints=[]) is None
