"""Array references encoded as Symbol names by the sym/pyneqsys ecosystem."""
import numpy as np
import pytest
import sympy as sp
import symox as sx
from symox.printing.numpy import NumPyPrinter


@pytest.mark.parametrize('name', ['x[..., 0]', 'x[..., 1]', 'x[:, -1]', 'x[0][1]', 'x[1:][..., 0]'])
def test_array_reference_symbols_execute_like_the_reference_printer(name):
    from sympy.printing.numpy import NumPyPrinter as ReferencePrinter
    x = np.array([[.2, .4], [.5, .8], [1.1, 1.3]], dtype=np.float64)
    a, b = sx.Symbol(name), sp.Symbol(name)
    expression = sx.sin(a)+sx.Rational(1, 3)*a**2
    reference = sp.sin(b)+sp.Rational(1, 3)*b**2
    source, expected = NumPyPrinter().doprint(expression), ReferencePrinter().doprint(reference)
    np.testing.assert_allclose(eval(source, {'numpy': np, 'x': x}),
                               eval(expected, {'numpy': np, 'x': x}), rtol=1e-14)


def test_indexed_symbols_keep_printer_subclass_dispatch():
    class ConsumerPrinter(NumPyPrinter):
        def _print_Symbol(self, expression):
            return super()._print_Symbol(expression)
    a, b = sx.Symbol('x[..., 0]'), sx.Symbol('x[..., 1]')
    source = ConsumerPrinter().doprint(a**2+b)
    inputs = np.array([[1., 2.], [3., 4.]])
    np.testing.assert_array_equal(eval(source, {'numpy': np, 'x': inputs}), [3., 13.])


@pytest.mark.parametrize('name', ['x[missing()', 'x[call()]', 'x.__class__[0]', 'lambda', 'x[0]; print(1)'])
def test_non_reference_fragments_keep_existing_validation(name):
    with pytest.raises(ValueError, match='invalid Python identifier'):
        NumPyPrinter().doprint(sx.Symbol(name))


@pytest.mark.parametrize('constructor', ['Matrix', 'ImmutableMatrix', 'SparseMatrix', 'ImmutableSparseMatrix'])
@pytest.mark.parametrize('shape', [(2, 2), (0, 3), (3, 0)])
def test_numpy_matrix_protocol_preserves_shape_and_entries(constructor, shape):
    values = list(range(shape[0]*shape[1]))
    matrix = getattr(sx, constructor)(*shape, values)
    reference = getattr(sp, constructor)(*shape, values)
    actual, expected = np.asarray(matrix), np.asarray(reference)
    assert actual.shape == expected.shape == shape
    assert actual.dtype == expected.dtype == object
    assert actual.tolist() == expected.tolist()
    np.testing.assert_array_equal(np.asarray(matrix, dtype=float), np.asarray(reference, dtype=float))
    with pytest.raises(TypeError, match='copy=False'):
        matrix.__array__(copy=False)


def test_numpy_matrix_protocol_copies_and_keeps_exact_symbolic_entries():
    x = sx.Symbol('x')
    matrix = sx.Matrix([[x, sx.Rational(1, 3)]])
    converted = np.asarray(matrix)
    assert converted[0, 0] == x
    assert converted[0, 1] == sx.Rational(1, 3)
    converted[0, 0] = sx.Integer(7)
    assert matrix[0, 0] == x


@pytest.mark.parametrize('factory', [lambda s, x, y: (x+1)**(-y),
    lambda s, x, y: (2*x)**(y+1), lambda s, x, y: (-x)**s.Rational(2, 3),
    lambda s, x, y: s.Pow(s.Pow(x, y, evaluate=False), 2, evaluate=False)])
def test_power_printer_preserves_base_and_exponent_precedence(factory):
    from sympy.printing.numpy import NumPyPrinter as ReferencePrinter
    x, y = sx.symbols('x y')
    rx, ry = sp.symbols('x y')
    expression = factory(sx, x, y)
    reference = factory(sp, rx, ry)
    values = {'numpy': np, 'x': np.array([.3, 1.2, 2.5]), 'y': np.array([2., 3., 4.])}
    with np.errstate(invalid='ignore'):
        actual = eval(NumPyPrinter().doprint(expression), values)
        expected = eval(ReferencePrinter().doprint(reference), values)
    np.testing.assert_allclose(actual, expected, rtol=1e-14, equal_nan=True)
