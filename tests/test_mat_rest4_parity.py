"""mat-rest4 parity: matrix derivatives by a MatrixSymbol, relational
matrix sides and immutable results of entrywise operations.

Expected values are hard-coded from SymPy 1.14.
"""
import pytest

from symox import (Eq, Identity, Inverse, Matrix, MatrixSymbol, N, SparseMatrix,
                   ImmutableSparseMatrix, Symbol, Trace, ZeroMatrix, cancel, expand,
                   factor, simplify, sqrt, srepr)

k = Symbol('k')
X = MatrixSymbol('X', k, k)
Y = MatrixSymbol('Y', k, k)
A = MatrixSymbol('A', k, k)
x = MatrixSymbol('x', k, 1)
y = MatrixSymbol('y', k, 1)
a = MatrixSymbol('a', k, 1)
b = MatrixSymbol('b', k, 1)


@pytest.mark.parametrize('expr, expected', [
    (lambda: Trace(X*X).diff(X), '2*X.T'),
    (lambda: Trace(X**2).diff(X), '2*X.T'),
    (lambda: Trace(X).diff(X), 'I'),
    (lambda: (x.T*X*x).diff(x), 'X*x + X.T*x'),
    (lambda: Trace(X.T*Y).diff(X), 'Y'),
    (lambda: Trace(A*X*Y).diff(X), 'A.T*Y.T'),
    (lambda: (x.T*y).diff(x), 'y'),
    (lambda: (x.T*x).diff(x), '2*x'),
    (lambda: (a.T*X*b).diff(X), 'a*b.T'),
    (lambda: Trace(Inverse(X)).diff(X), '-X.T**(-2)'),
    (lambda: Inverse(X).diff(X), 'PermuteDims(ArrayTensorProduct(-X.T**(-1), X**(-1)), (3)(1 2))'),
    (lambda: (X**2).diff(X), 'ArrayAdd(PermuteDims(ArrayTensorProduct(I, X), (3)(1 2)), '
                             'PermuteDims(ArrayTensorProduct(X.T, I), (3)(1 2)))'),
    (lambda: X.T.diff(X), 'PermuteDims(ArrayTensorProduct(I, I), (1 2 3))'),
    (lambda: (X.T*Y).diff(X), 'PermuteDims(ArrayTensorProduct(Y, I), (1 2 3))'),
    (lambda: (Y*X.T).diff(X), 'PermuteDims(ArrayTensorProduct(I, Y.T), (1 2 3))'),
    (lambda: (X.T*X).diff(X), 'ArrayAdd(PermuteDims(ArrayTensorProduct(X, I), (1 2 3)), '
                              'PermuteDims(ArrayTensorProduct(X, I), (3)(1 2)))'),
    (lambda: (x*x.T*x).diff(x), 'KroneckerProduct(I, x.T*x) + 2*x*x.T'),
    (lambda: sqrt(x.T*y).diff(x), '(1/2)*y*(x.T*y)**(-1/2)'),
    (lambda: sqrt(x.T*x).diff(x), 'x*(x.T*x)**(-1/2)'),
])
def test_matrix_derivatives(expr, expected):
    assert str(expr()) == expected


def test_scalar_free_of_variable_gives_zero_matrix():
    assert k.diff(A) == ZeroMatrix(k, k)


def test_scalar_matrix_function_coefficient_prints_its_head():
    assert str(Trace(A)*Identity(k)) == 'Trace(A)*I'


def test_relational_matrix_sides_are_dense_matrices():
    e = Eq(Matrix([[1, 2], [3, 4]]), Matrix([[1, 2], [3, 5]]), evaluate=False)
    assert type(e.lhs).__name__ == 'ImmutableDenseMatrix'
    assert type(e.rhs).__name__ == 'ImmutableDenseMatrix'
    assert srepr(e.lhs) == 'ImmutableDenseMatrix([[Integer(1), Integer(2)], [Integer(3), Integer(4)]])'
    z = Symbol('z')
    assert Eq(z + 1, 2, evaluate=False).lhs == z + 1


@pytest.mark.parametrize('function', [simplify, expand, factor, cancel, N])
def test_entrywise_operations_return_immutable_matrices(function):
    z = Symbol('z')
    assert type(function(SparseMatrix([[z, 0], [0, z + z]]))).__name__ == 'ImmutableSparseMatrix'
    assert type(function(Matrix([[z, 0], [0, z + z]]))).__name__ == 'ImmutableDenseMatrix'
    assert type(function(ImmutableSparseMatrix([[z, 0], [0, z]]))).__name__ == 'ImmutableSparseMatrix'
