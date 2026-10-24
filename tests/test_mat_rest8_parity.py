"""mat-rest8 (wave 14): element-wise matrix derivatives, held matrix
Derivative, DiagMatrix.doit and product entries.

Expected values are SymPy 1.14 outputs, except where a comment says SymPy is
wrong or raises.
"""
import pytest

from symox import *
del test
from symox.matrices.expressions.hadamard import (HadamardPower, HadamardProduct,
                                                 hadamard_power, hadamard_product)
from symox.matrices.expressions.matexpr import MatrixElement

i, j, k, m, n = symbols("i j k m n")
X = MatrixSymbol("X", k, k)
Y = MatrixSymbol("Y", k, k)
x = MatrixSymbol("x", k, 1)
y = MatrixSymbol("y", k, 1)
A = MatrixSymbol("A", k, k)
B = MatrixSymbol("B", k, k)
a = MatrixSymbol("a", k, 1)
b = MatrixSymbol("b", k, 1)


@pytest.mark.parametrize("make, expected", [
    (lambda: x.applyfunc(tan).diff(x),
     "DiagMatrix(Lambda(_d, tan(_d)**2 + 1).(x))"),
    (lambda: (A*x.applyfunc(exp)).diff(x),
     "DiagMatrix(Lambda(_d, exp(_d)).(x))*A.T"),
    (lambda: (x.T*A*x + k*y.applyfunc(sin).T*x).diff(x),
     "k*Lambda(_d, sin(_d)).(y) + A*x + A.T*x"),
    (lambda: (x.applyfunc(sin).T*y).diff(x),
     "ArrayDiagonal(Lambda(_d, cos(_d)).(x)*y.T, (0, 1))"),
    (lambda: (a.T*X*b).applyfunc(sin).diff(X),
     "a*Lambda(_d, cos(_d)).(a.T*X*b)*b.T"),
    (lambda: (a.T*X.applyfunc(sin)*b).diff(X),
     "DiagMatrix(a)*Lambda(_d, cos(_d)).(X)*DiagMatrix(b)"),
    (lambda: (a.T*(A*X*B).applyfunc(sin)*b).diff(X),
     "A.T*DiagMatrix(a)*Lambda(_d, cos(_d)).(A*X*B)*DiagMatrix(b)*B.T"),
    (lambda: (a.T*(A*X*b).applyfunc(sin)*b.T).diff(X),
     "ArrayTensorProduct(A.T*DiagMatrix(a)*Lambda(_d, cos(_d)).(A*X*b)*b.T, b)"),
    (lambda: (a.T*A*X.applyfunc(sin)*B*b).diff(X),
     "HadamardProduct(A.T*a*b.T*B.T, Lambda(_d, cos(_d)).(X))"),
    (lambda: (a.T*(A*X.applyfunc(sin)*B).applyfunc(log)*b).diff(X),
     "HadamardProduct(A.T*DiagMatrix(a)*Lambda(_d, 1/_d).(A*Lambda(_d, sin(_d)).(X)*B)"
     "*DiagMatrix(b)*B.T, Lambda(_d, cos(_d)).(X))"),
    (lambda: (a.T*X.applyfunc(sin).applyfunc(log)*b).diff(X),
     "DiagMatrix(a)*Lambda(_d, cos(_d)/sin(_d)).(X)*DiagMatrix(b)"),
    (lambda: (a.T*X.applyfunc(exp)*b).diff(X),
     "DiagMatrix(a)*Lambda(_d, exp(_d)).(X)*DiagMatrix(b)"),
    (lambda: X.applyfunc(sin).diff(X),
     "ArrayDiagonal(ArrayTensorProduct(Lambda(_d, cos(_d)).(X), I, I), (0, 3), (1, 5))"),
    (lambda: Trace(X.applyfunc(sin)).diff(X),
     "ArrayDiagonal(Lambda(_d, cos(_d)).(X), (0, 1))"),
    (lambda: (x.applyfunc(sin).T*x.applyfunc(cos)).diff(x),
     "ArrayAdd(ArrayDiagonal(Lambda(_d, cos(_d)).(x)*Lambda(_d, cos(_d)).(x.T), (0, 1)), "
     "ArrayDiagonal(Lambda(_d, sin(_d)).(x)*Lambda(_d, -sin(_d)).(x.T), (0, 1)))"),
])
def test_applyfunc_derivatives(make, expected):
    assert str(make()) == expected


@pytest.mark.parametrize("make, expected", [
    (lambda: hadamard_product(a, x, b).diff(x), "DiagMatrix(HadamardProduct(a, b))"),
    (lambda: (a.T*hadamard_product(A, X, B)*b).diff(X), "HadamardProduct(a*b.T, A, B)"),
    (lambda: hadamard_power(x, 2).diff(x), "2*DiagMatrix(x)"),
    (lambda: hadamard_power(x.T, 2).diff(x).doit(), "2*DiagMatrix(x)"),
    (lambda: hadamard_power(x, S.Half).diff(x), "(1/2)*DiagMatrix(HadamardPower(x, -1/2))"),
    (lambda: hadamard_power(a.T*X*b, 2).diff(X), "2*a*a.T*X*b*b.T"),
    (lambda: hadamard_product(x, x).diff(x), "2*DiagMatrix(x)"),
    (lambda: hadamard_product(X, A).diff(X),
     "ArrayDiagonal(ArrayTensorProduct(I, I, A), (1, 4), (3, 5))"),
    (lambda: hadamard_power(X, 2).diff(X),
     "ArrayAdd(ArrayDiagonal(ArrayTensorProduct(I, I, X), (1, 4), (3, 5)), "
     "ArrayDiagonal(ArrayTensorProduct(I, I, X), (1, 4), (3, 5)))"),
    (lambda: (a.T*hadamard_power(X, 3)*b).diff(X), "3*(HadamardProduct(a*b.T, HadamardPower(X, 2)))"),
    (lambda: Trace(hadamard_product(X, A)).diff(X), "ArrayDiagonal(A, (0, 1))"),
    (lambda: hadamard_product(A*x, x).diff(x), "DiagMatrix(A*x) + A.T*DiagMatrix(x)"),
    (lambda: hadamard_product(i*OneMatrix(k, 1), x, y).diff(i), "HadamardProduct(x, y)"),
    (lambda: hadamard_power(x, 2).diff(i), "0"),
])
def test_hadamard_derivatives(make, expected):
    assert str(make()) == expected


def test_hadamard_derivative_error():
    # SymPy: ValueError (diagonalizing indices of different dimensions).
    with pytest.raises(ValueError):
        (x.T*hadamard_product(x, y)).diff(x)


@pytest.mark.parametrize("make, expected", [
    (lambda: HadamardProduct(X.applyfunc(sin), A*B).doit(),
     "HadamardProduct(A*B, Lambda(_d, sin(_d)).(X))"),
    (lambda: hadamard_product(hadamard_power(X, 2), a*a.T),
     "HadamardProduct(a*a.T, HadamardPower(X, 2))"),
    (lambda: hadamard_product(A, hadamard_power(X, 2)), "HadamardProduct(HadamardPower(X, 2), A)"),
    (lambda: hadamard_product(X, X.applyfunc(sin)), "HadamardProduct(Lambda(_d, sin(_d)).(X), X)"),
    (lambda: (i**2*x).applyfunc(sin).diff(i),
     "HadamardProduct((2*i)*x, Lambda(_d, cos(_d)).(i**2*x))"),
])
def test_hadamard_canonical_order(make, expected):
    assert str(make()) == expected


@pytest.mark.parametrize("make, expected", [
    (lambda: DiagMatrix(S.Half*hadamard_power(x, -S.Half)).doit(),
     "(1/2)*DiagMatrix(HadamardPower(x, -1/2))"),
    (lambda: DiagMatrix(2*A*x).doit(), "2*DiagMatrix(A*x)"),
    (lambda: DiagMatrix(i*x).doit(), "i*DiagMatrix(x)"),
    (lambda: DiagMatrix((i*x).T).doit(), "i*DiagMatrix(x)"),
    (lambda: DiagMatrix(x.T).doit(), "DiagMatrix(x)"),
])
def test_diagmatrix_doit(make, expected):
    assert str(make()) == expected


def test_held_matrix_derivative():
    d = Derivative(X, X)
    assert str(d) == "Derivative(X, X)"
    assert isinstance(d, Derivative) and d.is_Derivative
    assert d.expr == X and d.variables == (X,)
    assert latex(d) == r"\frac{d}{d X} X"
    d2 = Derivative(X*A, X)
    assert str(d2) == "Derivative(X*A, X)"
    assert latex(d2) == r"\frac{\partial}{\partial X} X A"
    assert d2.free_symbols == {X, A}
    assert str(Derivative(X*A, X, Y)) == "Derivative(X*A, X, Y)"
    assert str(Derivative(x.T*A*x, x)) == "Derivative(x.T*A*x, x)"
    assert str(Derivative(i*X, i)) == "Derivative(i*X, i)"
    assert str(Derivative(X, (X, 2))) == "Derivative(X, (X, 2))"


def test_held_matrix_derivative_doit():
    # SymPy's doit gives 0 for Derivative(X, X) (wrong), raises TypeError for
    # Derivative(X*A, X) and NotImplementedError for Derivative(x.T*A*x, x).
    # symox gives the derivative (X.diff(X)).
    assert Derivative(X, X).doit() == X.diff(X)
    assert str(Derivative(X, X).doit()) == "PermuteDims(ArrayTensorProduct(I, I), (3)(1 2))"
    assert Derivative(X*A, X).doit() == (X*A).diff(X)
    assert str(Derivative(x.T*A*x, x).doit()) == "A*x + A.T*x"
    assert Derivative(X, X, evaluate=True) == X.diff(X)
    # Unsupported higher derivatives stay held.
    assert str(Derivative(X, (X, 2)).doit()) == "Derivative(X, (X, 2))"
    # SymPy values.
    assert str(Derivative(i*X, i).doit()) == "X"
    assert str(Derivative(i**2*X, (i, 2)).doit()) == "2*X"
    assert str(Derivative(i*X, i, j).doit()) == "0"
    assert str(Derivative(X, i).doit()) == "0"
    assert str(Derivative(X, Y).doit()) == "0"


def test_product_entry_partial_doit():
    C = MatrixSymbol("C", k, 2)
    D = MatrixSymbol("D", 2, k)
    E = MatrixSymbol("E", k, k)
    assert str((C*D*E)[i, j]) == (
        "Sum((C[i, 0]*D[0, _i_2] + C[i, 1]*D[1, _i_2])*E[_i_2, j], (_i_2, 0, k - 1))")
    assert str((E*C*D)[i, j]) == (
        "Sum(C[_i_1, _i_2]*D[_i_2, j]*E[i, _i_1], (_i_1, 0, k - 1), (_i_2, 0, 1))")
    assert str((C*D)[i, j]) == "C[i, 0]*D[0, j] + C[i, 1]*D[1, j]"
    assert str((D*E*C)[0, 1]) == (
        "Sum(C[_i_2, 1]*D[0, _i_1]*E[_i_1, _i_2], (_i_1, 0, k - 1), (_i_2, 0, k - 1))")


def test_product_entry_constructor_held():
    A2 = MatrixSymbol("A", 2, n)
    B2 = MatrixSymbol("B", n, 2)
    e = MatrixElement(A2*B2, 0, 0).subs(n, 2)
    assert str(e) == "(A*B)[0, 0]"
    assert str(e.doit()) == "A[0, 0]*B[0, 0] + A[0, 1]*B[1, 0]"
    F = MatrixSymbol("F", 2, 2)
    G = MatrixSymbol("G", 2, 2)
    assert str(MatrixElement(F*G, 0, 0)) == "(F*G)[0, 0]"
    assert str((F*G)[0, 0]) == "F[0, 0]*G[0, 0] + F[0, 1]*G[1, 0]"
