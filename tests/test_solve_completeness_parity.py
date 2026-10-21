"""Completeness of solve against the reference: products, perfect power
bases and partial polynomial roots."""

import cmath

import pytest
import sympy as sp

import symox as sx


NAMES = 'x a b y'


def build(module, source):
    namespace = dict(zip(NAMES.split(), module.symbols(NAMES)))
    for name in ('exp', 'log', 'sqrt', 'sin', 'pi', 'S', 'Rational'):
        namespace[name] = getattr(module, name)
    return eval(source, namespace), namespace['x']


def printed(solutions):
    values = []
    for solution in solutions:
        values.extend(solution.values() if isinstance(solution, dict) else [solution])
    return [str(value) for value in values]


def numeric(value):
    return complex(sp.sympify(value).evalf(30))


@pytest.mark.parametrize('source', [
    'x*log(x)', '(x - 1)*log(x - 1)', 'x**2*log(x)', 'x*log(x + 1)',
    'x*log(x)/(x - 1)', '(x - 2)*log(x)*sin(x)', 'x*exp(1/x)',
])
def test_product_keeps_the_roots_of_every_factor(source):
    expected = sp.solve(*build(sp, source))
    actual = sx.solve(*build(sx, source))
    assert printed(actual) == printed(expected)


@pytest.mark.parametrize('source', [
    '8**x - 4', '9**x - 3', '16**x - 2', '27**x - 9', '100**x - 10',
    '4**x - 2', '4**x - 3', '4**x - 4', '4**x - a', '4**x - y',
    '4**(x + 1) - 8', '4**(x**2) - 8', 'Rational(1, 9)**x - 3',
    'Rational(1, 4)**x - 2', 'Rational(1, 4)**x - 4', '36**x - 6',
    '4**x - 2**x', '4**x + 2**x - 6',
    'exp(2*x) - 2', 'exp(-2*x) - 2', 'exp(2*x + 1) - 2', 'exp(2*a*x) - 3',
])
def test_perfect_power_base_gives_the_reference_solutions(source):
    expected = sp.solve(*build(sp, source))
    actual = sx.solve(*build(sx, source))
    assert printed(actual) == printed(expected)


@pytest.mark.parametrize('source', [
    '4**x - 8', 'Rational(1, 4)**x - 8', '4**(-x) - 8', '4**(2*x) - 8',
    '4**x - 16', '4**x + 8', '8**x - 2', '9**x - 27', '2**(2*x) - 8',
    'Rational(1, 4)**x - 16',
])
def test_perfect_power_base_gives_every_branch(source):
    # The reference prints some logarithms combined by its simplification
    # step, e.g. (log(8)/2 + I*pi)/log(2) for 3/2 + I*pi/log(2) (see
    # docs/algorithmic-parity.md); the values and their order agree.
    f, x = build(sp, source)
    expected = [numeric(v) for v in printed(sp.solve(f, x))]
    solutions = printed(sx.solve(*build(sx, source)))
    actual = [numeric(v) for v in solutions]
    assert len(actual) == len(expected)
    for u, w in zip(actual, expected):
        assert abs(u - w) <= 1e-12*(1 + abs(w))
    for value in solutions:
        residual = complex(f.subs(x, sp.sympify(value)).evalf(30))
        assert cmath.isclose(residual, 0, abs_tol=1e-20)


def test_required_perfect_power_forms():
    x = sx.Symbol('x')
    assert sx.solve(4**x - 8, x) == [sx.Rational(3, 2), sx.Rational(3, 2) + sx.I*sx.pi/sx.log(2)]
    assert sx.solve(sx.Rational(1, 4)**x - 8, x) == [
        -sx.Rational(3, 2), -sx.Rational(3, 2) + sx.I*sx.pi/sx.log(2)]


@pytest.mark.parametrize('source', [
    'x**5 + a*x + 1', 'x**5 + a*x + b', '(x - 1)*(x**5 + a*x + 1)',
    'x**5 + sqrt(2)*x + 1', 'x**5 + pi*x + 1',
])
def test_incomplete_polynomial_returns_the_partial_roots(source):
    expected = sp.solve(*build(sp, source))
    actual = sx.solve(*build(sx, source))
    assert printed(actual) == printed(expected)


@pytest.mark.parametrize('source', ['x**5 + a*x + 1', 'x**5 + sqrt(2)*x + 1'])
def test_complete_polynomial_request_raises_when_roots_are_missing(source):
    with pytest.raises(NotImplementedError):
        sp.solve(*build(sp, source), incomplete=False)
    with pytest.raises(NotImplementedError):
        sx.solve(*build(sx, source), incomplete=False)


def test_rational_quintic_is_still_solved_exactly():
    # Partial roots never replace the exact CRootOf answer over QQ.
    expected = sp.solve(*build(sp, 'x**5 - x + 1'))
    actual = sx.solve(*build(sx, 'x**5 - x + 1'))
    assert printed(actual) == printed(expected)
