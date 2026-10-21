"""Forms of the univariate polynomial and product branches of solve.

The reference `_solve` takes the `roots` formulas (Cardano, Ferrari,
binomial) for a polynomial in one generator, falls back to `all_roots`
(radicals for linear, quadratic and binomial factors, CRootOf otherwise)
when they are incomplete, inverts the generator for every root, and solves
each factor of a product separately. These tests compare the printed
solutions with the reference.
"""
import cmath
import random

import pytest
import symox as sx
import sympy as sp


NAMES = 'x a b p q'


def build(module, source):
    namespace = dict(zip(NAMES.split(), module.symbols(NAMES)))
    for name in ('exp', 'log', 'sqrt', 'S'):
        namespace[name] = getattr(module, name)
    return eval(source, namespace), namespace['x']


def printed(solutions):
    values = []
    for solution in solutions:
        values.extend(solution.values() if isinstance(solution, dict) else [solution])
    return [str(value) for value in values]


@pytest.mark.parametrize('source,flags', [
    ('x**3 - 3*x + 1', {}),
    ('x**3 + x + 1', {}),
    ('x**3 - x - 1', {}),
    ('2*x**3 - 3*x**2 + 5', {}),
    ('x**3 - a', {}),
    ('x**3 + p*x + q', {}),
    ('x**3 - 1/a', {}),
    ('x**4 + 1', {}),
    ('x**4 - 2', {}),
    ('x**4 + x + 1', {}),
    ('x**6 + x**3 + 1', {}),
    ('x**5 - x + 1', {}),
    ('x**7 - x - 1', {}),
    ('x**3 - 3*x + 1', {'dict': True}),
    ('x**3 - 3*x + 1', {'check': False}),
    ('x**3 - 3*x + 1', {'cubics': False}),
    ('x**4 + x + 1', {'quartics': False}),
    ('(x**3 - 3*x + 1)/(x - 2)', {}),
    ('(x**3 - 2)*(x - 1)/(x - 2**(S(1)/3))', {}),
])
def test_polynomial_roots_follow_the_reference_formulas(source, flags):
    expected = sp.solve(*build(sp, source), **flags)
    actual = sx.solve(*build(sx, source), **flags)
    assert printed(actual) == printed(expected)


@pytest.mark.parametrize('source,flags', [
    ('(x**3 - 2)*(x**2 - a**2/(a + 1))', {}),
    ('(x**3 - 2)*(x**2 - a**2/(a + 1))', {'simplify': False}),
    ('x*(x**3 + x + 1)', {}),
    ('(x**2 + 1)*(x**3 - 3*x + 1)', {}),
    ('(x - 1)*(x**5 - x + 1)', {}),
    ('(x**2 - 2)*(x**5 - x + 1)', {}),
    ('(x**3 - 3*x + 1)*(x**5 - x + 1)', {}),
    ('(x**4 + 1)*(x - a)', {}),
    ('(x**2 + a)*(x**3 + a)', {}),
    ('(exp(x) - 1)*(exp(3*x) - 2)', {}),
])
def test_each_factor_of_a_product_is_solved_separately(source, flags):
    expected = sp.solve(*build(sp, source), **flags)
    actual = sx.solve(*build(sx, source), **flags)
    assert printed(actual) == printed(expected)


def test_product_factor_keeps_its_own_simplification():
    x, a = sx.symbols('x a')
    solutions = sx.solve((x**3 - 2)*(x**2 - a**2/(a + 1)), x)
    assert {'-a*sqrt(1/(a + 1))', 'a*sqrt(1/(a + 1))'} <= set(map(str, solutions))


@pytest.mark.parametrize('source', [
    'exp(x)**3 - 2', 'exp(3*x + 1) - 2', 'exp(3*x/2) - 2', 'exp(3*a*x) - 2',
    'exp(x + 1) - 2', 'exp(a*x + 1) - 2',
    '2**x - 8', '2**(x/2) - 3', '3**(2*x + 1) - 27', '2**(2*x) - 3*2**x + 2',
])
def test_generator_is_inverted_for_every_root(source):
    expected = sp.solve(*build(sp, source))
    actual = sx.solve(*build(sx, source))
    assert printed(actual) == printed(expected)


random.seed(3)
PARAMETERS = {name: complex(random.uniform(0.3, 1.7), random.uniform(-0.4, 0.4)) for name in 'abpq'}


def numeric(value):
    value = sp.sympify(str(value), locals={n: sp.Symbol(n) for n in NAMES.split()})
    value = value.subs({sp.Symbol(k): sp.Float(v.real, 30) + sp.I*sp.Float(v.imag, 30) for k, v in PARAMETERS.items()})
    return complex(value.evalf(30))


@pytest.mark.parametrize('source', [
    'exp(2*x) - 2', 'exp(-2*x) - 2', 'exp(-3*x) - 2', 'exp(2*x + 1) - 2',
    'exp(2*a*x) - 3', 'exp(x) - exp(2*x)/3 - 1', 'exp(3*x) - 3*exp(x) + 1',
    '2**(2*x) - 8', '2**(3*x) - 8', '2**(-2*x) - 8', '2**(3*x/2) - 3',
    'a**(2*x) - 2',
])
def test_power_of_exponential_generator_gives_every_branch(source):
    # The reference takes exp(p*x/q) as a power of exp(x/q); the forms of
    # the logarithms differ (see docs/algorithmic-parity.md), the values not.
    expected = [numeric(v) for v in printed(sp.solve(*build(sp, source)))]
    actual = [numeric(v) for v in printed(sx.solve(*build(sx, source)))]
    close = lambda u, w: abs(u - w) <= 1e-12*(1 + abs(u))
    assert len(actual) == len(expected)
    assert all(any(close(u, w) for w in actual) for u in expected)
    assert all(any(close(u, w) for w in expected) for u in actual)
    f, x = build(sp, source)
    for value in actual:
        residual = complex(f.subs({sp.Symbol(k): v for k, v in PARAMETERS.items()}).subs(x, value).evalf(30))
        assert cmath.isclose(residual, 0, abs_tol=1e-12)


def test_product_solutions_are_checked_against_simple_denominators_only():
    # The reference solves each factor of a product and only rejects the
    # roots of a denominator, so x*log(x) keeps the limit point 0.
    x = sx.Symbol('x')
    assert sx.solve(x*sx.log(x), x) == [0, 1]
    assert sx.solve((x - 1)*sx.log(x - 1), x) == [1, 2]


def test_incomplete_symbolic_quintic_returns_the_partial_roots():
    # With incomplete=True (the default) the reference returns the partial
    # result of roots(), which is empty here; incomplete=False still raises
    # because no root was proven absent.
    x, a = sx.symbols('x a')
    assert sx.solve(x**5 + a*x + 1, x) == []
    with pytest.raises(NotImplementedError):
        sx.solve(x**5 + a*x + 1, x, incomplete=False)


def test_cubic_in_a_power_generator_gives_three_solutions():
    # The reference prints the third solution in a form whose logarithm has
    # an exactly negative argument, which its evalf places at -I*pi; the
    # native value is on the principal branch (+I*pi). Both are solutions.
    f, x = build(sp, '2**(3*x) - 3*2**x + 1')
    actual = sx.solve(*build(sx, '2**(3*x) - 3*2**x + 1'))
    assert len(actual) == len(sp.solve(f, x)) == 3
    for value in printed(actual):
        assert abs(complex(f.subs(x, sp.sympify(value)).evalf(30))) < 1e-25


@pytest.mark.parametrize('source', [
    '2*c - 2*log(x) + log(y - 1) - log(y + 1)',
    'c - log(x) + log(y - 1)/2 - log(y + 1)/2',
])
def test_logarithm_generator_is_inverted_root_by_root(source):
    # A logarithm is not a polynomial generator of the reference, which
    # inverts it by _tsolve; substituting a symbolic inverse would print
    # (-x**2*exp(-2*c) - 1)/(x**2*exp(-2*c) - 1) instead.
    names = 'x y c'

    def solved(module):
        x, y, c = module.symbols(names)
        return module.solve(eval(source, {'x': x, 'y': y, 'c': c, 'log': module.log}), y)

    assert printed(solved(sx)) == printed(solved(sp))
