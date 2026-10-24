"""mat-rest5 parity: adjoint/transpose of matrices holding quantum operators,
of Abs/sign/powers of a noncommutative symbol, and nominal classes of the
matrix expressions.

Expected values are hard-coded from SymPy 1.14.
"""
import pickle

import pytest

from symox import *
del test
from symox.physics.quantum import Dagger, Operator

A = Symbol('A', commutative=False)
a = Operator('a')
X = MatrixSymbol('X', 2, 2)
Y = MatrixSymbol('Y', 2, 2)


@pytest.mark.parametrize('src, expected', [
    ("Matrix([[a, 1]]).adjoint()", "Matrix([[Dagger(a)], [1]])"),
    ("Matrix([[a, 1]]).H", "Matrix([[Dagger(a)], [1]])"),
    ("Matrix([[a, 1]]).T", "Matrix([[a], [1]])"),
    ("Matrix([[a, 1]]).conjugate()", "Matrix([[conjugate(a), 1]])"),
    ("Dagger(Matrix([[a, 1]]))", "Matrix([[Dagger(a)], [1]])"),
    ("Matrix([[a, 1]]).adjoint()[0] == Dagger(a)", "True"),
    ("Matrix([[a, 1]])*Matrix([[a], [1]])", "Matrix([[1 + a**2]])"),
    ("adjoint(Abs(A))", "Abs(A)"),
    ("transpose(Abs(A))", "Abs(A)"),
    ("conjugate(Abs(A))", "Abs(A)"),
    ("adjoint(sign(A))", "sign(conjugate(A))"),
    ("transpose(sign(A))", "sign(A)"),
    ("transpose(sqrt(A))", "transpose(sqrt(A))"),
    ("adjoint(A**2)", "Dagger(A)**2"),
    ("transpose(A**2)", "transpose(A)**2"),
    ("adjoint(A*a)", "Dagger(a)*Dagger(A)"),
    ("adjoint(a)", "Dagger(a)"),
])
def test_nc_adjoint_transpose(src, expected):
    assert str(eval(src)) == expected


@pytest.mark.parametrize('src, name, module', [
    ("X*Y", "MatMul", "matmul"),
    ("X+Y", "MatAdd", "matadd"),
    ("X**2", "MatPow", "matpow"),
    ("X.T", "Transpose", "transpose"),
    ("X.I", "Inverse", "inverse"),
    ("Adjoint(X)", "Adjoint", "adjoint"),
    ("Inverse(X)", "Inverse", "inverse"),
    ("Transpose(X)", "Transpose", "transpose"),
    ("BlockMatrix([[X, Y]])", "BlockMatrix", "blockmatrix"),
    ("BlockDiagMatrix(X, Y)", "BlockDiagMatrix", "blockmatrix"),
    ("Identity(2)", "Identity", "special"),
    ("ZeroMatrix(2, 2)", "ZeroMatrix", "special"),
    ("OneMatrix(2, 2)", "OneMatrix", "special"),
    ("HadamardProduct(X, Y)", "HadamardProduct", "hadamard"),
    ("KroneckerProduct(X, Y)", "KroneckerProduct", "kronecker"),
    ("DiagonalMatrix(X)", "DiagonalMatrix", "diagonal"),
    ("X[0:1, :]", "MatrixSlice", "slice"),
    ("MatPow(X, 3)", "MatPow", "matpow"),
    ("MatMul(X, Y, evaluate=False)", "MatMul", "matmul"),
    ("HadamardPower(X, 2)", "HadamardPower", "hadamard"),
    ("DiagMatrix(MatrixSymbol('v', 2, 1))", "DiagMatrix", "diagonal"),
])
def test_matrix_expression_classes(src, name, module):
    v = eval(src)
    t = type(v)
    assert t.__name__ == name
    assert t.__module__ == 'symox.matrices.expressions.' + module
    assert t is v.func
    assert isinstance(v, MatrixExpr)
    assert pickle.loads(pickle.dumps(v)) == v


def test_matrix_expression_class_methods():
    # class_key/fromiter keep working on the classes and on instances
    assert MatrixExpr.class_key() == MatrixExpr.class_key()
    assert MatMul.fromiter([X, Y]) == X*Y
    assert (X*Y).fromiter([X, Y]) == X*Y


def test_poly_entries_warn():
    from symox.matrices.normalforms import invariant_factors
    from symox.testing.pytest import warns_deprecated_sympy
    x = Symbol('x')
    with warns_deprecated_sympy():
        m = Matrix([[Poly(x-1), Poly(1, x), Poly(-1, x)],
                    [0, Poly(x), Poly(-1, x)],
                    [Poly(0, x), Poly(-1, x), Poly(x)]])
    assert invariant_factors(m, domain=QQ[x]) == (1, x - 1, x**2 - 1)


def test_rawmatrix_warns():
    from symox.polys.solvers import RawMatrix
    from symox.testing.pytest import warns_deprecated_sympy
    with warns_deprecated_sympy():
        RawMatrix([[1, 2], [3, 4]])


k = Symbol('k')
_A = MatrixSymbol('A', k, k)
_X = MatrixSymbol('X', k, k)
_Y = MatrixSymbol('Y', k, k)


@pytest.mark.parametrize('expr, expected', [
    (lambda: (Trace(_A)*_A).diff(_A),
     "ArrayAdd(ArrayTensorProduct(I, A), PermuteDims(ArrayTensorProduct(Trace(A)*I, I), (3)(1 2)))"),
    (lambda: (Trace(_A)*_X).diff(_A), "ArrayTensorProduct(I, X)"),
    (lambda: (Trace(_A)*_X*_Y).diff(_A), "ArrayTensorProduct(I, X*Y)"),
    (lambda: (2*Trace(_A)*_X).diff(_A), "ArrayTensorProduct(2*I, X)"),
    (lambda: (Trace(_A*_X)*_Y).diff(_A), "ArrayTensorProduct(X.T, Y)"),
    (lambda: ((Trace(_A)*_X).T).diff(_A), "ArrayTensorProduct(I, X.T)"),
    (lambda: (Trace(_A)*_X + _A).diff(_A),
     "ArrayAdd(ArrayTensorProduct(I, X), PermuteDims(ArrayTensorProduct(I, I), (3)(1 2)))"),
    (lambda: (Trace(_A)*Trace(_X)*_Y).diff(_A), "ArrayTensorProduct(Trace(X)*I, Y)"),
    (lambda: (Trace(_A)*_A*_X).diff(_A),
     "ArrayAdd(ArrayTensorProduct(I, A*X), PermuteDims(ArrayTensorProduct(Trace(A)*I, X), (3)(1 2)))"),
    (lambda: (_X*Trace(_A)).diff(_A), "ArrayTensorProduct(I, X)"),
])
def test_trace_times_matrix_derivative(expr, expected):
    assert str(expr()) == expected
