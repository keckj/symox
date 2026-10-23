"""Matrix class hierarchy, matrix Eq/Ne and sparse CSR helpers (SymPy 1.14).

Expected values were computed with SymPy 1.14.
"""
import pytest

from symox import *
del test
from symox.matrices import MatrixBase, MutableDenseMatrix, ImmutableDenseMatrix
from symox.matrices.dense import DenseMatrix
from symox.matrices.sparsetools import _csrtodok, _doktocsr

x, y = symbols('x y')
M = ImmutableMatrix([[1, 2], [3, 4]])


def test_dense_class_names():
    assert Matrix is MutableDenseMatrix
    assert ImmutableMatrix is ImmutableDenseMatrix
    assert type(Matrix([1])).__name__ == 'MutableDenseMatrix'
    assert type(ImmutableMatrix([1])).__name__ == 'ImmutableDenseMatrix'
    assert type(Matrix([1]).as_immutable()) is ImmutableDenseMatrix
    assert type(ImmutableMatrix([1]).as_mutable()) is MutableDenseMatrix
    assert type(eye(2)) is MutableDenseMatrix
    assert type(ImmutableMatrix.eye(2)) is ImmutableDenseMatrix


def test_immutable_results_keep_their_class():
    assert type(M + M) is ImmutableDenseMatrix
    assert type(M * Matrix([[1, 0], [0, 1]])) is ImmutableDenseMatrix
    assert type(M.T) is ImmutableDenseMatrix
    assert type(M.inv()) is ImmutableDenseMatrix


def test_dense_isinstance():
    assert isinstance(M, DenseMatrix)
    assert not isinstance(M, MutableDenseMatrix)
    assert not isinstance(Matrix([1]), ImmutableDenseMatrix)
    assert isinstance(M, Basic)
    assert not isinstance(Matrix([1]), Basic)
    assert isinstance(M, MatrixBase) and isinstance(Matrix([1]), MatrixBase)
    assert issubclass(ImmutableDenseMatrix, DenseMatrix)
    assert issubclass(MutableDenseMatrix, DenseMatrix)
    assert not issubclass(MutableDenseMatrix, ImmutableDenseMatrix)
    assert hash(ImmutableMatrix([1, 2])) == hash(ImmutableMatrix([1, 2]))
    with pytest.raises(TypeError):
        M[0, 0] = 5


def test_special_constructors_cls_option():
    assert type(eye(2, cls=ImmutableMatrix)) is ImmutableDenseMatrix
    assert type(zeros(2, cls=Matrix)) is MutableDenseMatrix
    assert type(Matrix.eye(2, cls=ImmutableMatrix)) is ImmutableDenseMatrix
    s = ones(2, 3, cls=SparseMatrix)
    assert type(s).__name__ == 'MutableSparseMatrix'
    assert s == ones(2, 3)
    e = Matrix._eval_eye(2, 3)
    assert type(e) is MutableDenseMatrix
    assert e == Matrix([[1, 0, 0], [0, 1, 0]])
    z = ImmutableMatrix._eval_zeros(1, 2)
    assert type(z) is ImmutableDenseMatrix
    assert z == Matrix([[0, 0]])
    assert Matrix._eval_ones(2, 1) == Matrix([[1], [1]])


def test_matrix_equality_relations():
    assert Eq(M, ImmutableMatrix([[1, 2], [3, 4]])) is S.true
    assert Eq(M, ImmutableMatrix([[1, 2], [3, 5]])) is S.false
    assert Eq(Matrix([[1, 2], [3, 4]]), M) is S.true
    assert Eq(M, ImmutableMatrix([1, 2, 3, 4])) is S.false
    assert Eq(M, x) is S.false
    assert Ne(M, M) is S.false
    assert Ne(M, ImmutableMatrix([1, 2])) is S.true
    e = Eq(ImmutableMatrix([x]), ImmutableMatrix([y]))
    assert str(e) == 'Eq(Matrix([[x]]), Matrix([[y]]))'
    assert e.subs(x, y) is S.true
    assert Eq(MatrixSymbol('B', 3, 3), M) is S.false
    A = MatrixSymbol('A', 2, 2)
    assert str(Eq(A, M)) == 'Eq(A, Matrix([\n[1, 2],\n[3, 4]]))'


def test_sparse_csr_helpers():
    a = SparseMatrix([[1, 2, 0, 0], [0, 3, 9, 0], [0, 1, 4, 0]])
    assert _doktocsr(a) == [[1, 2, 3, 9, 1, 4], [0, 1, 1, 2, 1, 2], [0, 2, 4, 6], [3, 4]]
    f = SparseMatrix(7, 8, {(2, 3): 5, (4, 5): 12})
    assert _doktocsr(f) == [[5, 12], [3, 5], [0, 0, 0, 1, 1, 2, 2, 2], [7, 8]]
    m = _csrtodok([[5, 7, 5], [2, 1, 3], [0, 1, 1, 3], [3, 4]])
    assert type(m).__name__ == 'MutableSparseMatrix'
    assert m == SparseMatrix(3, 4, {(0, 2): 5, (2, 1): 7, (2, 3): 5})


def test_trace_arg():
    X = MatrixSymbol('X', 2, 2)
    assert Trace(X).arg == X
    assert Q.positive(x).arg == x
