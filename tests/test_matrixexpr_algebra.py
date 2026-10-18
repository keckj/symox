"""Ordered extraction, local analysis and substitutions on matrix DAGs."""
from itertools import islice

import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize('index', range(18))
def test_operation_labels(index):
    x, y = sx.symbols('x y')
    n = sx.Symbol('n', integer=True, positive=True)
    a, b = sx.MatrixSymbol('A', n+1, n+1), sx.MatrixSymbol('B', n+1, n+1)
    cases = [a, a.T, a.I, a**2, -a, -2*a, a-b, -a-b, 2*x*a+b,
             (x/y)*a, 2*x/y*a, sx.Rational(-2, 3)*a,
             sx.exp(x)*a, sx.sin(x)**2*a, sx.sqrt(x)*a,
             a*b, (x**(-2)+sx.cos(x))*a,
             sx.Function('f')(x)*a]
    e = cases[index]
    assert e.count_ops() == sp.count_ops(e.to_sympy())
    assert e.count_ops(visual=True).to_sympy() == sp.count_ops(e.to_sympy(), visual=True)


@pytest.mark.parametrize('power', [-3, -1, 0, 1, 2, 5])
def test_leading_and_coefficient_extraction(power):
    x = sx.Symbol('x', positive=True)
    a, b = sx.MatrixSymbol('A', 2, 2), sx.MatrixSymbol('B', 2, 2)
    e = (2*x**power)*a
    coefficient, exponent = e.as_coeff_exponent(x)
    assert exponent == power and coefficient == 2*a
    assert (coefficient*x**exponent).doit() == e
    e = sx.exp(x)*x**power*a + x**(power+1)*b
    coefficient, exponent = e.leadterm(x)
    assert exponent == power and coefficient == a
    assert e.as_leading_term(x) == x**power*a


def test_ordered_extraction_and_collect():
    x = sx.Symbol('x')
    a, b, c = [sx.MatrixSymbol(name, 2, 2) for name in 'ABC']
    assert (a*b*c).extract_multiplicatively(a*b) == c
    assert (a*b*c).extract_multiplicatively(b*a) is None
    assert (2*x*a).extract_multiplicatively(x) == 2*a
    assert (2*a+b).extract_additively(a) == a+b
    assert (a+b).extract_additively(2*a) is None
    assert (a+b).extract_additively(c) is None
    assert a.extract_additively(a) == 0
    assert (a*b).as_coefficient(a) == b
    expression = (x*x+x)*a+2*x*b
    grouped = expression.collect(x, evaluate=False)
    assert grouped == {x*x: a, x: a+2*b}
    assert expression.collect(x).expand() == expression.expand()
    with pytest.raises(NotImplementedError):
        expression.collect(a)


def test_branch_factor_and_algebraic_type_convention():
    a = sx.MatrixSymbol('A', 2, 2)
    e = sx.exp_polar(4*sx.pi*sx.I)*a
    assert e.extract_branch_factor() == (a, sx.Integer(2))
    assert a.as_poly(sx.Symbol('x')) is None


@pytest.mark.parametrize('kind', ['exp', 'sin', 'rational', 'sqrt', 'log', 'essential'])
@pytest.mark.parametrize('point', [0, 1])
def test_meromorphic_coefficients(kind, point):
    x = sx.Symbol('x')
    a = sx.MatrixSymbol('A', 2, 2)
    scalar = {'exp': sx.exp(x), 'sin': sx.sin(x), 'rational': 1/x,
              'sqrt': sx.sqrt(x), 'log': sx.log(x), 'essential': sx.exp(1/x)}[kind]
    assert (scalar*a).is_meromorphic(x, point) == scalar.to_sympy().is_meromorphic(x.to_sympy(), point)


def test_hypergeometric_matrix_ratio_requires_common_scalar_ratio():
    k = sx.Symbol('k', integer=True, positive=True)
    a, b = sx.MatrixSymbol('A', 2, 2), sx.MatrixSymbol('B', 2, 2)
    assert (sx.factorial(k)*a).is_hypergeometric(k)
    assert (sx.factorial(k)*a+3*sx.factorial(k)*b).is_hypergeometric(k)
    assert not (a+sx.factorial(k)*b).is_hypergeometric(k)
    assert not (sx.sin(k)*a).is_hypergeometric(k)


def test_lazy_series_and_matrix_polynomial_inverse():
    x = sx.Symbol('x')
    a, b = sx.MatrixSymbol('A', 2, 2), sx.MatrixSymbol('B', 2, 2)
    assert list((a+x*x*b).lseries(x)) == [a, x*x*b]
    assert list(islice((sx.exp(x)*a).lseries(x), 4)) == [a, x*a, x*x/2*a, x**3/6*a]
    i, j = sx.symbols('i j', integer=True)
    matrix = sx.FunctionMatrix(2, 2, sx.Lambda((i, j), (x+1)/2+(x-1)/2*(-1)**(i+j)))
    assert matrix.as_explicit() == sx.Matrix([[x, 1], [1, x]])
    inverse = matrix.invert(x*x+1, x)
    product = matrix.as_explicit()*inverse
    for row in range(2):
        for col in range(2):
            assert sx.rem(product[row, col], x*x+1, x) == int(row == col)


def test_matrix_rounding_and_recursive_application():
    i, j = sx.symbols('i j', integer=True)
    e = sx.FunctionMatrix(2, 2, sx.Lambda((i, j), sx.Rational(1, 3)*(i+j)))
    assert e.round(2) == e.as_explicit().applyfunc(lambda entry: entry.round(2))
    a = sx.MatrixSymbol('A', 2, 2)
    with pytest.raises(TypeError):
        a.round()
    x, y = sx.symbols('x y')
    held = sx.MatMul(sx.Lambda(x, x+1), a)
    assert held.rcall(y).doit() == (y+1)*a


def test_matrix_and_scalar_replacements():
    x = sx.Symbol('x')
    a, b, c = [sx.MatrixSymbol(name, 2, 2) for name in 'ABC']
    e = sx.sin(x)*a*b
    result, mapping = e.replace(a, c, map=True)
    assert result.doit() == sx.sin(x)*c*b and mapping == {a: c}
    assert e.replace(sx.sin, sx.cos) == sx.cos(x)*a*b
    w = sx.Wild('w')
    assert e.replace(sx.sin(w), sx.cos(w)) == sx.cos(x)*a*b
    result = e.replace(sx.MatrixSymbol, lambda name, n, m: sx.MatrixSymbol(name.lower(), n, m))
    assert result == sx.sin(x)*sx.MatrixSymbol('a', 2, 2)*sx.MatrixSymbol('b', 2, 2)
    assert e.replace(lambda value: value == a, lambda value: c).doit() == sx.sin(x)*c*b
