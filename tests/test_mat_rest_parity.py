"""Matrix leftovers of SymPy 1.14 parity (task mat-rest).

Expected values are hard-coded from SymPy 1.14.
"""
import pytest

from symox import (
    CRootOf, I, ImmutableMatrix, Matrix, MatrixError,
    NonInvertibleMatrixError, NonSquareMatrixError, PolynomialError,
    Rational, Symbol, adjoint, cos, eye, ones, sin, solve, sqrt, symbols,
    zeros,
)

x = Symbol('x')
y = Symbol('y')

M5 = Matrix([[3, 0, 0, 0, -3], [0, -3, -3, 0, 3], [0, 3, 0, 3, 0],
             [0, 0, 3, 0, 3], [3, 0, 0, 3, 0]])
# SymPy prints the generator of the first cached PurePoly (x or _lambda);
# compare with the generator normalized.
ROOTS = ['3*CRootOf(x**5 - x**3 + 1, %d)' % i for i in range(5)]


def _norm(e):
    return str(e).replace('_lambda', 'x')


def test_eigenvals_crootof_dict_and_list():
    ev = M5.eigenvals()
    assert sorted(_norm(k) for k in ev) == ROOTS
    assert set(ev.values()) == {1}
    assert [_norm(e) for e in M5.eigenvals(multiple=True)] == ROOTS


def test_eigenvals_mixed_blocks():
    # charpoly lambda**5 + lambda + 1 = (l**2 + l + 1)*(l**3 - l**2 + 1)
    q = Matrix([[0, 0, 0, 0, -1], [1, 0, 0, 0, -1], [0, 1, 0, 0, 0],
                [0, 0, 1, 0, 0], [0, 0, 0, 1, 0]])
    ev = q.eigenvals()
    assert sum(ev.values()) == 5
    assert ev[Rational(-1, 2) - sqrt(3)*I/2] == 1
    assert not any(e.has(CRootOf) for e in ev)


def test_eigenvects_crootof():
    res = M5.eigenvects()
    assert [_norm(r[0]) for r in res] == ROOTS
    for val, mult, vects in res:
        assert mult == 1 and len(vects) == 1
        r = val.args[1] if val.args[0] == 3 else val.args[0]
        v = vects[0]
        assert list(v)[4] == 1
        assert (v[0] - (r**3 + r**4)).expand() == 0
        assert (v[1] - (r**3 + 1 - r)).expand() == 0
        assert (v[2] - (-r**4 + r**2 - r**3)).expand() == 0
        assert (v[3] - (-r**4 + r - r**3)).expand() == 0
        # M*v == val*v numerically
        w = (M5*v - val*v).evalf(30)
        assert all(abs(c) < 1e-20 for c in w)


def test_jordan_form_crootof_error():
    with pytest.raises(MatrixError) as e:
        M5.jordan_form()
    assert str(e.value) == ("Jordan normal form is not implemented if the "
                            "matrix have eigenvalues in CRootOf form")


def test_jordan_form_floats():
    P, J = Matrix([[1.0, 1], [0, 1.0]]).jordan_form()
    assert str(P) == 'Matrix([[1.00000000000000, 0], [0, 1.00000000000000]])'
    assert str(J) == ('Matrix([[1.00000000000000, 1.00000000000000], '
                      '[0, 1.00000000000000]])')
    res = Matrix([[1.5, 2], [0.5, 1]]).jordan_form()
    assert isinstance(res, list)
    P, J = res
    assert str(P) == ('Matrix([[-1.56155281280883, 2.56155281280883], '
                      '[1.00000000000000, 1.00000000000000]])')
    assert str(J) == 'Matrix([[0.219223593595585, 0], [0, 2.28077640640442]])'
    J = Matrix([[2.5, 1, 0], [0, 2.5, 0], [0, 0, 1.5]]).jordan_form(
        calc_transform=False)
    assert str(J) == ('Matrix([[1.50000000000000, 0, 0], [0, 2.50000000000000,'
                      ' 1.00000000000000], [0, 0, 2.50000000000000]])')


def test_bool_mask_indexing():
    M = Matrix(3, 3, range(9))
    assert M[[True, False, True], :] == Matrix([[0, 1, 2], [6, 7, 8]])
    assert M.extract([True, False, True], [False, True, True]) == \
        Matrix([[1, 2], [7, 8]])


def test_minor_cofactor_index_errors():
    M = Matrix(3, 3, range(9))
    msg = ('`i` and `j` must satisfy 0 <= i < ``M.rows`` (3)and 0 <= j < '
           '``M.cols`` (3).')
    with pytest.raises(ValueError) as e:
        M.minor(3, 0)
    assert str(e.value) == msg
    with pytest.raises(ValueError) as e:
        M.cofactor(0, 5)
    assert str(e.value) == msg
    assert M.minor(-1, -1) == -3
    assert M.minor_submatrix(-1, 0) == Matrix([[1, 2], [4, 5]])
    with pytest.raises(NonSquareMatrixError):
        Matrix(2, 3, range(6)).cofactor(0, 0)


def test_creation_dimension_errors():
    with pytest.raises(TypeError) as e:
        Matrix(2, 2)
    assert str(e.value) == ('\nData type not understood; expecting list of '
                            'lists or lists of values.')
    for f, msg in [(lambda: eye(3.0), '3.0 is not an integer'),
                   (lambda: ones(2, 3.0), '3.0 is not an integer'),
                   (lambda: Matrix(2.0, 2, [1, 2, 3, 4]),
                    '2.0 is not an integer'),
                   (lambda: zeros(-1, 2), 'Cannot create a -1 x 2 matrix. '
                    'Both dimensions must be positive')]:
        with pytest.raises(ValueError) as e:
            f()
        assert str(e.value) == msg


def test_is_identity_and_hessenberg():
    assert eye(3).is_Identity is True
    assert Matrix([[1, 0], [0, 2]]).is_Identity is False
    assert Matrix(2, 3, [1, 0, 0, 0, 1, 0]).is_Identity is False
    assert Matrix([[1, 2, 3], [4, 5, 6], [0, 7, 8]]).is_upper_hessenberg
    assert not Matrix([[1, 2, 3], [4, 5, 6], [1, 7, 8]]).is_upper_hessenberg
    assert Matrix([[1, 2, 0], [4, 5, 6], [1, 7, 8]]).is_lower_hessenberg


def test_solve_with_matrix_of_symbols():
    a, b = symbols('a b')
    X = Matrix([a, b])
    assert solve(Matrix([[1, 2], [3, 4]])*X - Matrix([1, 1]), X) == \
        {a: -1, b: 1}


def test_singular_symbolic_inverse_methods():
    m = Matrix([[-1, -1, 0], [x, 1, 1], [1, x, -1]])
    for f in (lambda: m.inv(method='LU'), lambda: m.inverse_LU(),
              lambda: m.inv(method='BLOCK')):
        with pytest.raises(NonInvertibleMatrixError) as e:
            f()
        assert str(e.value) == 'Matrix det == 0; not invertible.'
    A = Matrix([[x, 1], [1, x]])
    expected = Matrix([[-x/(1 - x**2), x**2/(1 - x**2) + 1],
                       [1/(1 - x**2), -x/(1 - x**2)]])
    assert A.inv(method='LU') == expected
    assert A.inv(method='BLOCK') == expected
    assert Matrix([[1, 2], [3, 4]]).inv(method='BLOCK') == \
        Matrix([[-2, 1], [Rational(3, 2), Rational(-1, 2)]])
    assert ImmutableMatrix([[1, x], [0, 1]]).inv(method='LU') == \
        Matrix([[1, -x], [0, 1]])


def test_noncommutative_inverse_and_adjoint():
    A, B, C, D = symbols('A B C D', commutative=False)
    M = Matrix([[A, B], [C, D]])
    with pytest.raises(PolynomialError) as e:
        M.inv()
    assert str(e.value) == 'non-commutative expressions are not supported'
    assert str(M.inv(method='LU')) == (
        'Matrix([[A**(-1)*(1 + B*(-C*A**(-1)*B + D)**(-1)*C*A**(-1)), '
        '-A**(-1)*B*(-C*A**(-1)*B + D)**(-1)], [-(-C*A**(-1)*B + D)**(-1)'
        '*C*A**(-1), (-C*A**(-1)*B + D)**(-1)]])')
    from symox import conjugate, transpose
    assert adjoint(sin(A)) == transpose(sin(conjugate(A)))
    assert adjoint(cos(A)*B) == adjoint(B)*transpose(cos(conjugate(A)))
    assert adjoint(sin(x)) == sin(x).conjugate()


def test_dotprodsimp_context_manager():
    from symox.matrices import dotprodsimp as dps_public
    from symox.matrices.utilities import (
        dotprodsimp, _get_intermediate_simp_bool)
    assert dps_public is dotprodsimp
    with dotprodsimp(True):
        assert _get_intermediate_simp_bool(False) is True
    with dotprodsimp(False):
        assert _get_intermediate_simp_bool(True) is False
    assert _get_intermediate_simp_bool(False) is False
    with dotprodsimp(True):
        assert Matrix([[x, 1], [1, x]])**2 == Matrix([[x**2 + 1, 2*x],
                                                      [2*x, x**2 + 1]])


def test_hermitian_fallback_and_dirac():
    assert Matrix([[1, 1], [I, 1]]).is_hermitian is False
    assert Matrix([[1, 1 - I], [1 + I, 2]]).is_hermitian is True
    assert Matrix([[Symbol('a'), 1], [1, 2]]).is_hermitian is None
    assert Matrix([[Symbol('a', real=True), I], [-I, 2]]).is_hermitian is True
    assert Matrix(4, 1, [1, I, 2, 3]).D == Matrix([[1, -I, -2, -3]])
    with pytest.raises(AttributeError):
        Matrix(3, 1, [1, 2, 3]).D


def test_sparse_list_keys():
    from symox import SparseMatrix
    S = SparseMatrix(3, 3, range(9))
    assert S[[0, 0, 2], 1] == Matrix([[1], [1], [7]])
    assert S[[True, False, True], [2, 0]] == Matrix([[2, 0], [8, 6]])
    assert S[1, [2, 1]] == Matrix([[5, 4]])


def test_find_reasonable_pivot_naive():
    from symox.matrices.determinant import _find_reasonable_pivot_naive
    assert _find_reasonable_pivot_naive([0, x, 3]) == (0, 0, True, [])
    assert _find_reasonable_pivot_naive([0, 0]) == (0, 0, True, [])


def test_ragged_creation_and_diag():
    from symox import diag
    for data, msg in [([[1, 2], [3]], 'mismatched dimensions'),
                      ([1, 2, [3]], 'expecting list of lists')]:
        with pytest.raises(ValueError) as e:
            Matrix(data)
        assert str(e.value) == msg
    with pytest.raises(ValueError) as e:
        diag([[1, 2], [3]])
    assert str(e.value) == 'mismatched dimensions'
    assert diag([1, 2, [3]], strict=False) == Matrix([[1], [2], [3]])
    assert diag([1, [2, 3]], strict=False) == Matrix([[1, 0], [2, 3]])
    assert diag([[1, 2], [3]], strict=False) == Matrix([[1, 2], [3, 0]])
    assert Matrix([[]]).shape == (1, 0)
    assert Matrix([[], []]).shape == (2, 0)
