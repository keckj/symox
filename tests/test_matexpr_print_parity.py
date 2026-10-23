"""String forms of matrix expressions, with expected values from SymPy 1.14."""
import pytest

from symox import (Adjoint, BlockDiagMatrix, BlockMatrix, HadamardPower,
                   HadamardProduct, I, Identity, ImmutableMatrix, Inverse,
                   KroneckerProduct, MatPow, Matrix, MatrixSymbol, Rational,
                   Transpose, ZeroMatrix, sin, sstr, symbols)
from symox import _matrixexpr_methods

n, x = symbols('n x')
A, B, C, D = [MatrixSymbol(s, n, n) for s in 'ABCD']
X = MatrixSymbol('X', 3, 3)
Y = MatrixSymbol('Y', 2, 2)

CASES = [
    (lambda: Inverse(A*B), '(A*B)**(-1)'),
    (lambda: Inverse(A + B), '(A + B)**(-1)'),
    (lambda: Inverse(A**2), '(A**2)**(-1)'),
    (lambda: Inverse(A.T), 'A.T**(-1)'),
    (lambda: Inverse(Inverse(A)), '(A**(-1))**(-1)'),
    (lambda: Inverse(A**Rational(1, 2)), '(A**(1/2))**(-1)'),
    (lambda: Inverse(Identity(n) + A), '(I + A)**(-1)'),
    (lambda: Inverse(-A + B), '(-A + B)**(-1)'),
    (lambda: Inverse(KroneckerProduct(A, B)), '(KroneckerProduct(A, B))**(-1)'),
    (lambda: Inverse(HadamardPower(A, 2)), '(HadamardPower(A, 2))**(-1)'),
    (lambda: (A.T*B).I, 'B**(-1)*A.T**(-1)'),
    (lambda: C.T*C.I, 'C.T*C**(-1)'),
    (lambda: MatPow(A, Rational(1, 2)), 'A**(1/2)'),
    (lambda: MatPow(A, Rational(-1, 2)), 'A**(-1/2)'),
    (lambda: A**Rational(3, 2), 'A**(3/2)'),
    (lambda: (A + B)**Rational(1, 2), '(A + B)**(1/2)'),
    (lambda: MatPow(A, x), 'A**x'),
    (lambda: MatPow(A, 2*x), 'A**(2*x)'),
    (lambda: MatPow(A, x**2), 'A**(x**2)'),
    (lambda: MatPow(A, sin(x)), 'A**sin(x)'),
    (lambda: MatPow(A, -2), 'A**(-2)'),
    (lambda: MatPow(A, 0.5), 'A**0.5'),
    (lambda: MatPow(A, -0.5), 'A**(-0.5)'),
    (lambda: MatPow(A, I), 'A**I'),
    (lambda: (A*B)**2, '(A*B)**2'),
    (lambda: (A**2).T, 'A.T**2'),
    (lambda: Transpose(Inverse(A)), '(A**(-1)).T'),
    (lambda: Transpose(HadamardProduct(A, B)), '(HadamardProduct(A, B)).T'),
    (lambda: Adjoint(A*B), 'Adjoint(A*B)'),
    (lambda: Adjoint(Inverse(A)), 'Adjoint(A**(-1))'),
    (lambda: A**Rational(-1, 2)*B, 'A**(-1/2)*B'),
    (lambda: (A*B.T + C.I)**Rational(1, 3), '(C**(-1) + A*B.T)**(1/3)'),
    (lambda: Rational(1, 2)*A, '(1/2)*A'),
    (lambda: x*A, 'x*A'),
    (lambda: (x + 1)*A, '(x + 1)*A'),
    (lambda: -x*A, '(-x)*A'),
    (lambda: -(x + 1)*A, '(-x - 1)*A'),
    (lambda: -I*A, '-I*A'),
    (lambda: -2*I*A, '-(2*I)*A'),
    (lambda: (1 + I)*A, '(1 + I)*A'),
    (lambda: -A*B, '-A*B'),
    (lambda: A*HadamardProduct(B, C), 'A*(HadamardProduct(B, C))'),
    (lambda: KroneckerProduct(A, B)*C, '(KroneckerProduct(A, B))*C'),
    (lambda: A + B - C, 'A + B - C'),
    (lambda: A - B*D.I*C, 'A - B*D**(-1)*C'),
    (lambda: X[0:2, 1:3], 'X[:2, 1:]'),
    (lambda: X[::2, 1:], 'X[::2, 1:]'),
    (lambda: ImmutableMatrix([[1, 2], [3, 4]])*Y, 'Matrix([\n[1, 2],\n[3, 4]])*Y'),
    (lambda: BlockMatrix([[A, B], [C, D]]), 'Matrix([\n[A, B],\n[C, D]])'),
    (lambda: BlockMatrix([[A]]), 'Matrix([[A]])'),
    (lambda: BlockMatrix([[A, B, C]]), 'Matrix([[A, B, C]])'),
    (lambda: BlockMatrix([[A], [B]]), 'Matrix([\n[A],\n[B]])'),
    (lambda: BlockMatrix([[A.I, B*C], [C + D, D**2]]),
     'Matrix([\n[A**(-1),  B*C],\n[  C + D, D**2]])'),
    (lambda: BlockMatrix([[A, ZeroMatrix(n, n)], [Identity(n), D]]), 'Matrix([\n[A, 0],\n[I, D]])'),
    (lambda: BlockMatrix([[A, B], [C, D]]).T, 'Matrix([\n[A.T, C.T],\n[B.T, D.T]])'),
    (lambda: Inverse(BlockMatrix([[A, B], [C, D]])), 'Matrix([\n[A, B],\n[C, D]])**(-1)'),
    (lambda: BlockDiagMatrix(A, B), 'Matrix([\n[A, 0],\n[0, B]])'),
    (lambda: BlockDiagMatrix(A, B, C), 'Matrix([\n[A, 0, 0],\n[0, B, 0],\n[0, 0, C]])'),
    (lambda: BlockMatrix([[A, B]])*BlockMatrix([[C], [D]]), 'Matrix([[A, B]])*Matrix([\n[C],\n[D]])'),
]

# The native printer keeps the internal order of MatAdd terms; every other
# case matches SymPy without the Python printer.
NATIVE_ORDER = {'A + B - C': '-C + A + B', 'A - B*D**(-1)*C': '-B*D**(-1)*C + A'}


@pytest.mark.parametrize('build, expected', CASES)
def test_str_repr_sstr(build, expected):
    value = build()
    assert str(value) == expected
    assert repr(value) == expected
    assert sstr(value) == expected


@pytest.mark.parametrize('build, expected', CASES)
def test_native_str(build, expected):
    value = build()
    native = _matrixexpr_methods._native_str
    assert native(value) == NATIVE_ORDER.get(expected, expected)


def test_sstr_explicit_matrices():
    assert sstr(Matrix([[1, 22], [3, 4]])) == 'Matrix([\n[1, 22],\n[3,  4]])'
    assert sstr(ImmutableMatrix([[1, 2]])) == 'Matrix([[1, 2]])'
    assert sstr(Matrix(0, 2, [])) == 'Matrix(0, 2, [])'


def test_block_schur_and_lu_doctests():
    m, n = symbols('m n')
    A = MatrixSymbol('A', n, n)
    B = MatrixSymbol('B', n, m)
    C = MatrixSymbol('C', m, n)
    D = MatrixSymbol('D', m, m)
    X = BlockMatrix([[A, B], [C, D]])
    expected = {
        ('A', False): '-C*A**(-1)*B + D',
        ('A', True): '-C*A**(-1)*B + D',
        ('B', False): 'C - D*B**(-1)*A',
        ('B', True): 'C - D*(B.T*B)**(-1)*B.T*A',
        ('C', False): '-A*C**(-1)*D + B',
        ('C', True): '-A*(C.T*C)**(-1)*C.T*D + B',
        ('D', False): 'A - B*D**(-1)*C',
        ('D', True): 'A - B*D**(-1)*C',
    }
    for (pivot, generalized), text in expected.items():
        assert str(X.schur(pivot, generalized=generalized)) == text
    L, U = X.LUdecomposition()
    assert str(L) == ('Matrix([\n[   A**(1/2),                         0],\n'
                      '[C*A**(-1/2), (-C*A**(-1)*B + D)**(1/2)]])')
    assert str(U) == ('Matrix([\n[A**(1/2),               A**(-1/2)*B],\n'
                      '[       0, (-C*A**(-1)*B + D)**(1/2)]])')


def test_block_diag_lu_doctest():
    m, n = symbols('m n')
    X = BlockDiagMatrix(MatrixSymbol('A', n, n), MatrixSymbol('D', m, m))
    L, U = X.LUdecomposition()
    expected = 'Matrix([\n[A**(1/2),        0],\n[       0, D**(1/2)]])'
    assert str(L) == expected
    assert str(U) == expected
