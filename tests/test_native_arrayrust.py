"""Native array-element derivatives, matrix-expression argument checks and
the native matrix-line derivative kernel (crates/symox-core
matrixexpr/lines.rs). Expected values are hard-coded from SymPy 1.14."""
import pytest

import symox._matrix_derivatives as md
from symox import (Array, Determinant, Identity, Inverse, KroneckerDelta, Matrix, MatrixSymbol,
                   OneMatrix, Rational, SympifyError, Trace, cos, diff, sin, symbols)
from symox import _core
from symox.matrices.expressions.hadamard import HadamardProduct
from symox.tensor.array.expressions import ArrayElement, ArraySymbol

i, j, k, l, a, x_ = symbols('i j k l a x')
AS = ArraySymbol('A', (3, 2))
BS = ArraySymbol('B', (3, 2))


def kind(value):
    return type(value).__name__


# -- ArrayElement derivatives (SymPy ArrayElement._eval_derivative) ---------

@pytest.mark.parametrize('make,expected', [
    (lambda: AS[i, j].diff(AS[k, l]), 'KroneckerDelta(i, k)*KroneckerDelta(j, l)'),
    (lambda: AS[i, j].diff(AS[i, j]), '1'),
    (lambda: AS[i, j].diff(BS[k, l]), '0'),
    (lambda: AS[i, j].diff(x_), '0'),
    (lambda: (x_*AS[i, j]).diff(AS[k, l]), 'x*KroneckerDelta(i, k)*KroneckerDelta(j, l)'),
    (lambda: (AS[i, j]**2).diff(AS[i, j]), '2*A[i, j]'),
    (lambda: sin(AS[0, 1]).diff(AS[0, 1]), 'cos(A[0, 1])'),
    (lambda: diff(AS[0, 1]*AS[1, 1], AS[1, 1]), 'A[0, 1]'),
    (lambda: AS[0, 1].diff(AS[1, 1]), '0'),
    (lambda: AS[0, 1].diff(BS[0, 1]), '0'),
    (lambda: (x_*AS[0, 1]).diff(AS[0, 1]), 'x'),
    (lambda: AS[i, j].diff(AS[k, l], 2), '0'),
    (lambda: ArrayElement(AS, (i, j)).diff(ArrayElement(AS, (k, l))),
     'KroneckerDelta(i, k)*KroneckerDelta(j, l)'),
    (lambda: (AS[i, j]*BS[k, l]).diff(BS[k, l]), 'A[i, j]'),
])
def test_array_element_derivative(make, expected):
    assert str(make()) == expected


def test_array_element_derivative_values():
    # Formerly 0 (the chain rule over the name and the indices).
    d = (x_*AS[i, j]).diff(AS[k, l])
    assert d == x_*KroneckerDelta(i, k)*KroneckerDelta(j, l)
    assert d.subs({k: i, l: j}) == x_
    assert Matrix([AS[0, 0], x_]).diff(AS[0, 0]) == Matrix([[1], [0]])
    # SymPy prints the KroneckerDelta factors before the ArrayElement
    # (A[i, j]*KroneckerDelta(...) here: the sort order of the element node
    # is a recorded gap); the values are SymPy's.
    assert (AS[i, j]**2).diff(AS[k, l]) == 2*KroneckerDelta(i, k)*KroneckerDelta(j, l)*AS[i, j]._symox_()
    assert (AS[i, j]*AS[k, l]).diff(AS[0, 1]) == (
        KroneckerDelta(0, i)*KroneckerDelta(1, j)*AS[k, l]._symox_()
        + KroneckerDelta(0, k)*KroneckerDelta(1, l)*AS[i, j]._symox_())
    assert cos(AS[i, j]).diff(AS[k, l]) == -sin(AS[i, j]._symox_())*KroneckerDelta(i, k)*KroneckerDelta(j, l)


# -- argument checks --------------------------------------------------------

n_ = symbols('n')
MA = MatrixSymbol('A', n_, n_)
MB = MatrixSymbol('B', n_, n_)


@pytest.mark.parametrize('dims,shown', [
    ((-1, 2), '-1'), ((2.0, 2), '2.00000000000000'), ((2j, 2), '2.0*I'),
    ((2, -1), '-1'), ((2, 2.0), '2.00000000000000'), ((2, 2j), '2.0*I'),
    ((symbols('m', integer=False), 2), 'm'), ((symbols('m', negative=True), 2), 'm'),
])
def test_matrix_symbol_dimension_check(dims, shown):
    with pytest.raises(ValueError, match='The dimension specification %s should be a nonnegative integer.'
                       % shown.replace('*', r'\*').replace('.', r'\.')):
        MatrixSymbol('A', *dims)


def test_matrix_symbol_string_dimension():
    with pytest.raises(SympifyError):
        MatrixSymbol(1, 2, 'A')


def test_hadamard_scalar_argument_and_transpose():
    for scalar in (n_, 1):
        with pytest.raises(TypeError, match='Mix of Matrix and Scalar symbols'):
            HadamardProduct(MA, scalar)
    assert HadamardProduct(MA, MB, MA).args == (MA, MB, MA)
    assert HadamardProduct(MA, MB, MA).T.args == (MA.T, MB.T, MA.T)
    assert str(HadamardProduct(MA, MB, MA).doit()) == 'HadamardProduct(HadamardPower(A, 2), B)'


@pytest.mark.parametrize('make', [
    lambda: Array([[1, 2], [3, 4]]) - Matrix([[1, 2], [3, 4]]),
    lambda: Array([[1, 2], [3, 4]]) + Matrix([[1, 2], [3, 4]]),
    lambda: Matrix([[1, 2], [3, 4]]) - Array([[1, 2], [3, 4]]),
    lambda: Matrix([[1, 2], [3, 4]]) + Array([[1, 2], [3, 4]]),
])
def test_array_matrix_addition_is_refused(make):
    with pytest.raises(TypeError):
        make()


# -- matrix-line derivatives --------------------------------------------------

kk = symbols('k')
X = MatrixSymbol('X', kk, kk)
A = MatrixSymbol('A', kk, kk)
B = MatrixSymbol('B', kk, kk)
xv = MatrixSymbol('x', kk, 1)
yv = MatrixSymbol('y', kk, 1)

# (expression, variable, SymPy 1.14 type, SymPy 1.14 str)
LINES = [
    (lambda: xv.T*xv, xv, 'MatMul', '2*x'),
    (lambda: xv.T*A*xv, xv, 'MatAdd', 'A*x + A.T*x'),
    (lambda: xv.T*A*yv, xv, 'MatMul', 'A*y'),
    (lambda: yv.T*A*xv, xv, 'MatMul', 'A.T*y'),
    (lambda: a*xv, xv, 'MatMul', 'a*I'),
    (lambda: A*xv, xv, 'Transpose', 'A.T'),
    (lambda: xv.T*A, xv, 'MatrixSymbol', 'A'),
    (lambda: Trace(X), X, 'Identity', 'I'),
    (lambda: Trace(A*X), X, 'Transpose', 'A.T'),
    (lambda: Trace(X*A), X, 'Transpose', 'A.T'),
    (lambda: Trace(A*X*B), X, 'MatMul', 'A.T*B.T'),
    (lambda: Trace(X.T), X, 'Identity', 'I'),
    (lambda: Trace(A*X.T*B), X, 'MatMul', 'B*A'),
    (lambda: Trace(X*X), X, 'MatMul', '2*X.T'),
    (lambda: Trace(X**2), X, 'MatMul', '2*X.T'),
    (lambda: Trace(X**3), X, 'MatMul', '3*X.T**2'),
    (lambda: Trace(A*X**2*B), X, 'MatAdd', 'A.T*B.T*X.T + X.T*A.T*B.T'),
    (lambda: Trace(X.T*A*X), X, 'MatAdd', 'A*X + A.T*X'),
    (lambda: Trace(X*A*X.T), X, 'MatAdd', 'X*A + X*A.T'),
    (lambda: Trace(X.T*X), X, 'MatMul', '2*X'),
    (lambda: Trace(Inverse(X)), X, 'MatMul', '-X.T**(-2)'),
    (lambda: Trace(X**-1), X, 'MatMul', '-X.T**(-2)'),
    (lambda: Trace(A*X**-1*B), X, 'MatMul', '-X.T**(-1)*A.T*B.T*X.T**(-1)'),
    (lambda: Trace(X)**2, X, 'MatMul', '(2*Trace(X))*I'),
    (lambda: Trace(X)*Trace(A*X), X, 'MatAdd', 'Trace(A*X)*I + Trace(X)*A.T'),
    (lambda: 2*Trace(X), X, 'MatMul', '2*I'),
    (lambda: Trace(X)**3, X, 'MatMul', '(3*Trace(X)**2)*I'),
    (lambda: Trace(X)*a, X, 'MatMul', 'a*I'),
    (lambda: xv.T*X*xv, X, 'MatMul', 'x*x.T'),
    (lambda: xv.T*X.T*yv, X, 'MatMul', 'y*x.T'),
    (lambda: xv.T*X**2*yv, X, 'MatAdd', 'x*y.T*X.T + X.T*x*y.T'),
    (lambda: xv.T*X**-1*yv, X, 'MatMul', '-X.T**(-1)*x*y.T*X.T**(-1)'),
    (lambda: xv.T*Inverse(X)*yv, X, 'MatMul', '-X.T**(-1)*x*y.T*X.T**(-1)'),
    (lambda: xv.T*A*X*B*yv, X, 'MatMul', 'A.T*x*y.T*B.T'),
    (lambda: xv.T*xv*a, xv, 'MatMul', '(2*a)*x'),
    (lambda: Trace(X*Identity(kk)), X, 'Identity', 'I'),
    (lambda: Trace(OneMatrix(kk, kk)*X), X, 'OneMatrix', '1'),
    (lambda: xv.T*(A + B)*xv, xv, 'MatAdd', '(A + B)*x + (A.T + B.T)*x'),
    (lambda: (A + B)*xv, xv, 'MatAdd', 'A.T + B.T'),
    (lambda: (xv + yv).T*(xv + yv), xv, 'MatAdd', '2*x + 2*y'),
    (lambda: Trace(Identity(kk)*X*X), X, 'MatMul', '2*X.T'),
    (lambda: Trace(-X), X, 'MatMul', '-I'),
    (lambda: Trace(3*A*X), X, 'MatMul', '3*A.T'),
    (lambda: xv.T*yv, xv, 'MatrixSymbol', 'y'),
    (lambda: yv.T*xv, xv, 'MatrixSymbol', 'y'),
    (lambda: xv.T*Identity(kk)*xv, xv, 'MatMul', '2*x'),
    (lambda: Trace(X + A), X, 'Identity', 'I'),
    (lambda: Trace((X + A)*(X + B)), X, 'MatAdd', 'A.T + B.T + 2*X.T'),
    (lambda: Trace(A*X*B*X), X, 'MatAdd', 'A.T*X.T*B.T + B.T*X.T*A.T'),
    (lambda: Trace(A*X.T*B*X), X, 'MatAdd', 'B*X*A + B.T*X*A.T'),
    (lambda: Trace(X.T*X*X), X, 'MatAdd', 'X**2 + X*X.T + X.T*X'),
    (lambda: Trace(X)*Determinant(A), X, 'MatMul', 'Determinant(A)*I'),
    (lambda: xv.T*X*X.T*xv, X, 'MatMul', '2*x*x.T*X'),
]


@pytest.mark.parametrize('make,var,typ,expected', LINES)
def test_matrix_line_derivative_matches_sympy(make, var, typ, expected):
    expr = make()
    native = _core._matrix_lines_derivative(expr, var)
    assert native is not None
    assert kind(native) == typ and str(native) == expected
    d = expr.diff(var)
    assert kind(d) == typ and str(d) == expected and d == native


class _PythonOnly:
    """The `_core` module seen by the Python algorithm, without the kernel."""

    def __getattr__(self, name):
        if name == '_matrix_lines_derivative':
            return lambda expr, x: None
        return getattr(_core, name)


@pytest.mark.parametrize('make,var,typ,expected', LINES)
def test_matrix_line_derivative_matches_python_algorithm(make, var, typ, expected, monkeypatch):
    expr = make()
    native = _core._matrix_lines_derivative(expr, var)
    monkeypatch.setattr(md, '_core', _PythonOnly())
    assert md.matrix_derivative(expr, var) == native


@pytest.mark.parametrize('make,var', [
    (lambda: A*X, X),                       # rank four (array expression)
    (lambda: X.T*X, X),
    (lambda: Inverse(X), X),
    (lambda: Trace(X**symbols('p')), X),    # symbolic power
    (lambda: Determinant(X), X),
    (lambda: sin(Trace(X)), X),
    (lambda: Trace(HadamardProduct(X, A)), X),   # element-wise: array path
    (lambda: (xv.T*xv)**Rational(1, 2), xv),      # power of a 1x1 matrix
])
def test_matrix_line_kernel_declines(make, var):
    assert _core._matrix_lines_derivative(make(), var) is None


def test_declined_power_of_1x1_keeps_python_result():
    d = ((xv.T*xv)**Rational(1, 2)).diff(xv)
    assert kind(d) == 'MatMul' and str(d) == 'x*(x.T*x)**(-1/2)'
