"""Matrix features of the c3 audit, with expected values from SymPy 1.14."""
import pytest

from symox import (BlockMatrix, Float, ImmutableMatrix, KroneckerProduct, Matrix,
                   MatPow, MatrixSymbol, NonInvertibleMatrixError, S, ShapeError,
                   SparseMatrix, block_collapse, eye, ones, symbols)
from symox.abc import x, y
from symox.matrices.expressions import combine_kronecker
from symox.utilities.misc import func_name

m, n = symbols('m n')
A = MatrixSymbol('A', n, n)
B = MatrixSymbol('B', n, m)
C = MatrixSymbol('C', m, n)
D = MatrixSymbol('D', m, m)
X = BlockMatrix([[A, B], [C, D]])


def test_block_lu_decomposition():
    L, U = X.LUdecomposition()
    assert L.blocks[1, 0] == C*MatPow(A, -S.Half)
    assert L.blocks[1, 1] == (D - C*A.I*B)**S.Half
    assert U.blocks[0, 1] == MatPow(A, -S.Half)*B
    assert block_collapse(L*U) == X
    with pytest.raises(ShapeError, match='Block LU decomposition is supported only'):
        BlockMatrix([[A]]).LUdecomposition()


def test_generalized_schur():
    assert X.schur() == D - C*A.I*B
    assert X.schur('D') == A - B*D.I*C
    assert X.schur('B', generalized=True) == C - D*(B.T*B).inv()*B.T*A
    assert X.schur('C', generalized=True) == B - A*(C.T*C).inv()*C.T*D
    with pytest.raises(ShapeError, match='Schur Complement can only be calculated'):
        BlockMatrix([[A]]).schur()


def test_matpow_inverse():
    assert repr((A**2).I) == 'A**(-2)'
    assert (A**S.Half).inv() == MatPow(A, -S.Half)
    assert A**S.Half*(A**S.Half).I == MatPow(A, 0).doit()


def test_kronecker_shape_predicates():
    k, l = symbols('k l', integer=True)
    P = MatrixSymbol('P', k, l)
    Q = MatrixSymbol('Q', l, k)
    assert KroneckerProduct(P, Q).has_matching_shape(KroneckerProduct(Q, P)) is True
    assert KroneckerProduct(P, Q).has_matching_shape(KroneckerProduct(P, Q)) is False
    assert KroneckerProduct(P, Q).has_matching_shape(P) is False
    E = MatrixSymbol('E', k, k)
    F = MatrixSymbol('F', l, l)
    G = MatrixSymbol('G', k, k)
    H = MatrixSymbol('H', l, l)
    assert KroneckerProduct(E, F).structurally_equal(KroneckerProduct(G, H)) is True
    assert KroneckerProduct(E, F).structurally_equal(KroneckerProduct(H, G)) is False
    assert KroneckerProduct(E, F).structurally_equal(G) is False
    assert X.structurally_equal(BlockMatrix([[A, B], [C, D]]))


def test_combine_kronecker():
    k, l = symbols('k l', integer=True)
    P = MatrixSymbol('P', k, l)
    Q = MatrixSymbol('Q', l, k)
    E = MatrixSymbol('E', l, l)
    F = MatrixSymbol('F', k, k)
    assert combine_kronecker(KroneckerProduct(P, Q)*KroneckerProduct(Q, P)) == KroneckerProduct(P*Q, Q*P)
    assert combine_kronecker(KroneckerProduct(E, F)**k) == KroneckerProduct(E**k, F**k)
    # SymPy gives KroneckerProduct(P + Q.T, Q + P.T), which is not equal.
    s = KroneckerProduct(P, Q) + KroneckerProduct(Q.T, P.T)
    assert combine_kronecker(s) == s
    # Non-integer powers are not distributed.
    h = KroneckerProduct(E, F)**S.Half
    assert combine_kronecker(h) == h


def test_matrix_has_class():
    a = Matrix(((1, x), (0.2, 3)))
    b = SparseMatrix(((1, x), (0.2, 3)))
    for mat in (a, b):
        assert mat.has(x) is True
        assert mat.has(y) is False
        assert mat.has(Float) is True
    assert Matrix([[1, x]]).has(Float) is False


def test_diag_cls_name():
    assert func_name(Matrix.diag(1, cls=ImmutableMatrix)) == 'ImmutableDenseMatrix'
    assert func_name(Matrix.diag(1)) == 'MutableDenseMatrix'
    assert func_name(ImmutableMatrix([1])) == 'ImmutableDenseMatrix'


def test_copyin_scalar_key():
    M = Matrix([[0, 1], [2, 3], [4, 5]])
    for make in (eye, lambda k: SparseMatrix(eye(k))):
        I = make(3)
        I[:3, :2] = M
        assert I.tolist() == [[0, 1, 0], [2, 3, 0], [4, 5, 1]]
        I[0, 1] = M
        assert I.tolist() == [[0, 0, 1], [2, 2, 3], [4, 4, 5]]
        I = make(3)
        I[(0, 1)] = [[7], [8]]
        assert I.tolist() == [[1, 7, 0], [0, 8, 0], [0, 0, 1]]
        I = make(3)
        I[1] = M
        assert I.tolist() == [[1, 0, 1], [0, 2, 3], [0, 4, 5]]
        I = make(3)
        I[0, 0] = [1, 2]
        assert I.tolist() == [[1, 0, 0], [2, 1, 0], [0, 0, 1]]
        with pytest.raises(ShapeError):
            make(3)[2, 2] = M


def test_sparse_creation():
    with pytest.raises(ValueError, match=r'The length of the flat list \(2\) does not match the specified size \(2 \* 2\)\.'):
        SparseMatrix(2, 2, [1, 2])
    assert SparseMatrix(2, 2, [[1, 2]]).tolist() == [[1, 2], [0, 0]]
    with pytest.raises(ValueError, match=r'The location \(0, 2\) is out of the designated range\[0, 1\]x\[0, 1\]'):
        SparseMatrix(2, 2, [[1, 2, 3]])
    assert SparseMatrix(None, [[1, 2, 3]]).tolist() == [[1, 2, 3]]
    with pytest.raises(ValueError, match=r'There is a collision at \(1, 1\) for 2 and 1\.'):
        SparseMatrix(3, 3, {(0, 0): ones(2), (1, 1): 2})
    with pytest.raises(ValueError, match='Pass rows=None and no cols for autosizing.'):
        SparseMatrix(2, None, [1])
    assert SparseMatrix(2, 2, lambda i, j: i + j).tolist() == [[0, 1], [1, 2]]
    assert SparseMatrix(None, {(1, 1): 3}).shape == (2, 2)
    assert SparseMatrix(2, 2, {(0, 0): [[1, 2], [3, 4]]}).tolist() == [[1, 2], [3, 4]]
    assert SparseMatrix(2, 3, range(6)).tolist() == [[0, 1, 2], [3, 4, 5]]
