"""Matrix leftovers of mat-rest; expected values come from SymPy 1.14."""
import pytest

from symox import *
from symox.matrices.normalforms import smith_normal_form, hermite_normal_form
from symox.polys.domains import ZZ, QQ, EXRAW
from symox.geometry import Point

del test

a, b, n, x, y = symbols('a b n x y')


def test_non_integer_matrix_powers():
    assert Matrix([[33, 24], [48, 57]])**S.Half == Matrix([[5, 2], [4, 7]])
    assert Matrix([[2, 1], [1, 2]])**Rational(1, 2) == Matrix([
        [S.Half + sqrt(3)/2, -S.Half + sqrt(3)/2],
        [-S.Half + sqrt(3)/2, S.Half + sqrt(3)/2]])
    assert Matrix([[1, a], [0, 1]])**n == Matrix([[1, a*n], [0, 1]])
    assert Matrix([[b, a], [0, b]])**n == Matrix([[b**n, a*b**(n - 1)*n], [0, b**n]])
    assert Matrix([[2, 1], [1, 2]])**10.0 == Matrix([[29525, 29524], [29524, 29525]])


def test_singular_matrix_powers():
    M = Matrix([[0, 1], [0, 0]])
    with pytest.raises(NonInvertibleMatrixError):
        M**Rational(21, 10)
    p = M**n
    assert p.is_MatPow
    assert p.exp == n


def test_jordan_block_helpers():
    J = Matrix([[2, 1], [0, 2]])
    assert J._eval_matrix_exp_jblock() == Matrix([[exp(2), exp(2)], [0, exp(2)]])
    assert J._eval_matrix_log_jblock() == Matrix([[log(2), S.Half], [0, log(2)]])


def test_float_mod_by_integer_is_exact_zero():
    assert type(Float(10.0) % 1).__name__ == 'Zero'
    assert Float(10.5) % 1 == Float(0.5)
    assert (Float(10.0) % Float(2)).is_Float
    assert Mod(Float(10.0), 3) == Float(1.0)


def test_normal_forms():
    A = Matrix([[12, 6, 4], [3, 9, 6], [2, 16, 14]])
    assert smith_normal_form(A, domain=ZZ) == Matrix([[1, 0, 0], [0, 10, 0], [0, 0, 30]])
    assert hermite_normal_form(A) == Matrix([[10, 0, 2], [0, 15, 3], [0, 0, 2]])


def test_diff_by_matrix_and_array():
    A = Matrix([[x, y], [y, x]])
    d = A.diff(A)
    assert type(d).__name__ == 'MutableDenseNDimArray'
    assert d.tolist() == [[[[1, 0], [0, 1]], [[0, 1], [1, 0]]],
                          [[[0, 1], [1, 0]], [[1, 0], [0, 1]]]]
    assert type(ImmutableMatrix(A).diff(A)).__name__ == 'ImmutableDenseNDimArray'
    assert diff(x**2*y, [[x, y]]).tolist() == [2*x*y, x**2]
    assert diff(x**2*y, ([x, y], 2)).tolist() == [[2*y, 2*x], [2*x, 0]]


def test_mod_of_matrix():
    assert Mod(Matrix([-1, 1, 2]), 3) == Matrix([2, 1, 2])
    assert Mod(Matrix([10, 20, -15]), 4) == Matrix([2, 0, 1])
    assert Mod(Matrix([x, 1]), 3) == Matrix([Mod(x, 3), 1])


def test_match_matrix_symbol():
    A = MatrixSymbol('A', 3, 3)
    B = MatrixSymbol('B', 3, 3)
    w = Wild('w')
    v = Wild('v', exclude=[A])
    assert A.match(w) == {w: A}
    assert (A*B).match(w) == {w: A*B}
    assert A.match(v) is None
    assert B.match(v) == {v: B}


def test_epsilon_eq():
    assert Float(1.).epsilon_eq(Float(1.0000000001)) is S.false
    assert Float(1.).epsilon_eq(1) is S.true


def test_point_rows():
    m = Matrix([Point(0, 0) - Point(1, 2), Point(1, -sqrt(3))])
    assert m == Matrix([[-1, -2], [1, -sqrt(3)]])


def test_deferred_vector_not_iterable():
    assert not iterable(DeferredVector('X'))
    with pytest.raises(TypeError):
        Matrix(DeferredVector('V'))


def test_rep_domains():
    from symox.polys.matrices import DM
    for Mat in (Matrix, ImmutableMatrix, SparseMatrix, ImmutableSparseMatrix):
        M = Mat([[1, 2], [3, 4]])
        assert M._rep == DM([[1, 2], [3, 4]], ZZ).to_sparse()
        assert (M/2)._rep.domain == QQ
        assert Mat([[x, 1]])._rep.domain == EXRAW
        assert Mat._fromrep(DM([[1, 2], [3, 4]], ZZ)) == Matrix([[1, 2], [3, 4]])


def test_pinv_ed_symbolic():
    c, d = symbols('c d')
    A = Matrix([[a, b], [c, d]])
    P = A.pinv(method='ED')
    vals = {a: 1, b: 2, c: 3, d: 5}
    assert all(abs(N(e.subs(vals))) < 1e-12 for e in P - A.inv())
    B = Matrix([[1, 2], [3, 4], [5, 6]])
    assert all(abs(e.evalf()) < 1e-12 for e in B.pinv(method='ED') - B.pinv())


def test_expand_hints_on_matrix():
    M = Matrix([[x*(x + 1), 2], [log(x*y), 1]])
    assert expand_mul(M) == Matrix([[x**2 + x, 2], [log(x*y), 1]])
    assert expand_log(M, force=True) == Matrix([[x*(x + 1), 2], [log(x) + log(y), 1]])


def test_pinv_eigendecomposition():
    B = Matrix([[1, 2], [2, 4]])
    assert B.pinv(method='ED') == Matrix([[Rational(1, 25), Rational(2, 25)],
                                          [Rational(2, 25), Rational(4, 25)]])
    A = Matrix([[a, b], [b, a]])
    assert simplify(A.pinv(method='ED') - A.inv()) == zeros(2, 2)


def test_matrix_rep():
    from symox.polys.matrices import DM
    assert Matrix([[1, 2], [3, 4]])._rep == DM([[1, 2], [3, 4]], ZZ).to_sparse()
    assert (Matrix([[1, 2], [3, 4]])/2)._rep.domain == QQ
    assert Matrix([[x, 1]])._rep.domain == EXRAW
    assert SparseMatrix([[1, 2]])._rep.domain == ZZ
    M = Matrix._fromrep(DM([[1, 2], [3, 4]], ZZ).to_sparse())
    assert M == Matrix([[1, 2], [3, 4]])


def test_scalar_diff_by_matrix_stays_matrix():
    res = x.diff(Matrix([[x, y]]))
    assert isinstance(res, ImmutableDenseMatrix)
    assert res == Matrix([[1, 0]])
    res = diff(x**3, Matrix([[x, y]]))
    assert isinstance(res, ImmutableDenseMatrix)
    assert res == Matrix([[3*x**2, 0]])
    assert diff(x, SparseMatrix([[x, y]])) == Matrix([[1, 0]])
    assert diff(x*y, Array([x, y])).tolist() == [y, x]
