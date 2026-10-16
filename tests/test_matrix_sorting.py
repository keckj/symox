"""Presentation ordering must include the typed matrix parent and its shape."""
import pytest
import sympy as sp
import symox as sx


def normalized(key):
    if isinstance(key, tuple):
        return tuple(normalized(item) for item in key)
    return key.to_sympy() if isinstance(key, sx.Expr) else key


@pytest.mark.parametrize("operation", [
    lambda s, A, B: A,
    lambda s, A, B: A[0, 1],
    lambda s, A, B: s.Trace(A),
    lambda s, A, B: s.Determinant(A),
    lambda s, A, B: 2*s.Trace(A),
    lambda s, A, B: s.Trace(A)**3,
    lambda s, A, B: -s.Determinant(A)/2,
    lambda s, A, B: s.MatAdd(B, A, evaluate=False),
    lambda s, A, B: s.MatAdd(A, 2*B, evaluate=False),
    lambda s, A, B: s.MatAdd(-B, A, evaluate=False),
    lambda s, A, B: s.MatAdd(A, s.Symbol('x')*B, evaluate=False),
    lambda s, A, B: s.MatAdd(B*A, A*B, evaluate=False),
    lambda s, A, B: s.MatAdd(s.Symbol('x')*s.Symbol('y')*B, s.Symbol('y')*A.T, evaluate=False),
    lambda s, A, B: s.MatAdd(2*s.Symbol('x')*B, s.Symbol('x')*A, evaluate=False),
    lambda s, A, B: s.MatMul(2, B, A, evaluate=False),
    lambda s, A, B: s.MatPow(A, 3, evaluate=False),
    lambda s, A, B: s.Inverse(A),
    lambda s, A, B: s.Transpose(A),
    lambda s, A, B: s.Adjoint(A),
    lambda s, A, B: s.conjugate(A),
    lambda s, A, B: s.Identity(2),
    lambda s, A, B: s.ZeroMatrix(2, 3),
    lambda s, A, B: s.OneMatrix(2, 3),
    lambda s, A, B: s.HadamardProduct(B, A),
    lambda s, A, B: s.HadamardPower(A, 2),
    lambda s, A, B: s.HadamardPower(2, A),
    lambda s, A, B: s.HadamardPower(A, B),
    lambda s, A, B: s.KroneckerProduct(B, A),
    lambda s, A, B: s.MatrixExpr(s.ImmutableMatrix([[1, 2], [3, 4]])) if s is sx else s.ImmutableMatrix([[1, 2], [3, 4]]),
])
def test_matrix_sort_keys_match_sympy(operation):
    actual = operation(sx, sx.MatrixSymbol('A', 2, 2), sx.MatrixSymbol('B', 2, 2))
    expected = operation(sp, sp.MatrixSymbol('A', 2, 2), sp.MatrixSymbol('B', 2, 2))
    assert normalized(sx.default_sort_key(actual)) == sp.default_sort_key(expected)


@pytest.mark.parametrize("operation", [
    lambda A: A,
    lambda A: A[0, 0],
    lambda A: sx.Trace(A),
    lambda A: sx.Determinant(A),
])
def test_parent_and_dimensions_distinguish_sort_keys(operation):
    values = [operation(sx.MatrixSymbol(name, n, n)) for name, n in [('A', 2), ('B', 2), ('A', 3)]]
    keys = [sx.default_sort_key(value) for value in values]
    assert all(keys[i] != keys[j] for i in range(3) for j in range(i))


def test_scalar_matrix_printing_uses_the_same_parent_order():
    A, B = [sx.MatrixSymbol(name, 2, 2) for name in ['A', 'B']]
    a, b = [sp.MatrixSymbol(name, 2, 2) for name in ['A', 'B']]
    for expression, reference in [
        (sx.Trace(B)+sx.Trace(A), sp.Trace(b)+sp.Trace(a)),
        (sx.Determinant(B)+sx.Determinant(A), sp.Determinant(b)+sp.Determinant(a)),
        (B[0, 0]+A[0, 0], b[0, 0]+a[0, 0]),
        (2*sx.Trace(B)+sx.Trace(A), 2*sp.Trace(b)+sp.Trace(a)),
        (sx.Trace(B)*sx.Trace(A), sp.Trace(b)*sp.Trace(a)),
    ]:
        assert str(expression) == str(reference)
