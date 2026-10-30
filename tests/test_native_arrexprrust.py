"""Native rank-four matrix derivatives (crates/symox-core matrixexpr/lines.rs,
``matrix_lines_rank_four``): SymPy's array form of a matrix derived by a
matrix. Expected strings are hard-coded from SymPy 1.14; every case is also
checked against the Python ``_rank_four`` path that ran before the port."""
import pytest

import symox._matrix_derivatives as md
from symox import Identity, Inverse, MatrixSymbol, Trace, symbols
from symox import _core
from symox.tensor.array.expressions import ArrayAdd, ArrayTensorProduct, PermuteDims

A, B, C, X, R = (MatrixSymbol(s, 3, 3) for s in 'ABCXR')
x = MatrixSymbol('x', 3, 1)
a = MatrixSymbol('a', 3, 1)
P = MatrixSymbol('P', 2, 3)
Q = MatrixSymbol('Q', 3, 2)
W = MatrixSymbol('W', 2, 3)
n, t = symbols('n t')
M = MatrixSymbol('M', n, n)
N = MatrixSymbol('N', n, n)

# (case, SymPy 1.14 str)
SYMPY = [
    ('X.diff(X)',
     'PermuteDims(ArrayTensorProduct(I, I), (3)(1 2))'),
    ('X.T.diff(X)',
     'PermuteDims(ArrayTensorProduct(I, I), (1 2 3))'),
    ('(A*X).diff(X)',
     'PermuteDims(ArrayTensorProduct(A.T, I), (3)(1 2))'),
    ('(X*A).diff(X)',
     'PermuteDims(ArrayTensorProduct(I, A), (3)(1 2))'),
    ('(A*X*B).diff(X)',
     'PermuteDims(ArrayTensorProduct(A.T, B), (3)(1 2))'),
    ('(A*X.T*B).diff(X)',
     'PermuteDims(ArrayTensorProduct(B, A.T), (1 2 3))'),
    ('(X*X).diff(X)',
     'ArrayAdd(PermuteDims(ArrayTensorProduct(I, X), (3)(1 2)), PermuteDims(ArrayTensorProduct(X.T, I), (3)(1 2)))'),
    ('(X**2).diff(X)',
     'ArrayAdd(PermuteDims(ArrayTensorProduct(I, X), (3)(1 2)), PermuteDims(ArrayTensorProduct(X.T, I), (3)(1 2)))'),
    ('(X.T*X).diff(X)',
     'ArrayAdd(PermuteDims(ArrayTensorProduct(X, I), (1 2 3)), PermuteDims(ArrayTensorProduct(X, I), (3)(1 2)))'),
    ('(X*X.T).diff(X)',
     'ArrayAdd(PermuteDims(ArrayTensorProduct(I, X.T), (3)(1 2)), PermuteDims(ArrayTensorProduct(I, X.T), (1 2 3)))'),
    ('(A*X+X*B).diff(X)',
     'ArrayAdd(PermuteDims(ArrayTensorProduct(A.T, I), (3)(1 2)), PermuteDims(ArrayTensorProduct(I, B), (3)(1 2)))'),
    ('(X+A).diff(X)',
     'PermuteDims(ArrayTensorProduct(I, I), (3)(1 2))'),
    ('(2*A*X).diff(X)',
     'PermuteDims(ArrayTensorProduct(2*A.T, I), (3)(1 2))'),
    ('(t*X).diff(X)',
     'PermuteDims(ArrayTensorProduct(t*I, I), (3)(1 2))'),
    ('(Trace(X)*A).diff(X)',
     'ArrayTensorProduct(I, A)'),
    ('(Trace(A*X)*B).diff(X)',
     'ArrayTensorProduct(A.T, B)'),
    ('(X**-1).diff(X)',
     'PermuteDims(ArrayTensorProduct(-X.T**(-1), X**(-1)), (3)(1 2))'),
    ('(A*X**-1*B).diff(X)',
     'PermuteDims(ArrayTensorProduct(-X.T**(-1)*A.T, X**(-1)*B), (3)(1 2))'),
    ('(X*A*X).diff(X)',
     'ArrayAdd(PermuteDims(ArrayTensorProduct(I, A*X), (3)(1 2)), PermuteDims(ArrayTensorProduct(A.T*X.T, I), (3)(1 2)))'),
    ('(P*W.T).diff(W)',
     'PermuteDims(ArrayTensorProduct(I, P.T), (1 2 3))'),
    ('(Q*W).diff(W)',
     'PermuteDims(ArrayTensorProduct(Q.T, I), (3)(1 2))'),
    ('(W*Q).diff(W)',
     'PermuteDims(ArrayTensorProduct(I, Q), (3)(1 2))'),
    ('(M*N).diff(N)',
     'PermuteDims(ArrayTensorProduct(M.T, I), (3)(1 2))'),
    ('(N.T*M*N).diff(N)',
     'ArrayAdd(PermuteDims(ArrayTensorProduct(M*N, I), (1 2 3)), PermuteDims(ArrayTensorProduct(M.T*N, I), (3)(1 2)))'),
    ('(x*x.T).diff(x)',
     'ArrayAdd(PermuteDims(ArrayTensorProduct(I, x.T), (3)(1 2)), PermuteDims(ArrayTensorProduct(I, x.T), (1 2 3)))'),
    ('(Identity(3)*X*A).diff(X)',
     'PermuteDims(ArrayTensorProduct(I, A), (3)(1 2))'),
    ('(X.T*A.T).diff(X)',
     'PermuteDims(ArrayTensorProduct(A.T, I), (1 2 3))'),
    ('((A*X).T).diff(X)',
     'PermuteDims(ArrayTensorProduct(A.T, I), (1 2 3))'),
    ('(A*X*B + C).diff(X)',
     'PermuteDims(ArrayTensorProduct(A.T, B), (3)(1 2))'),
    ('(Inverse(X)*A).diff(X)',
     'PermuteDims(ArrayTensorProduct(-X.T**(-1), X**(-1)*A), (3)(1 2))'),
    ('(A*Inverse(X.T)).diff(X)',
     'PermuteDims(ArrayTensorProduct(X.T**(-1), -X**(-1)*A.T), (1 2 3))'),
    ('(t**2*X*A).diff(X)',
     'PermuteDims(ArrayTensorProduct(t**2*I, A), (3)(1 2))'),
    ('(A*B*X*C).diff(X)',
     'PermuteDims(ArrayTensorProduct(B.T*A.T, C), (3)(1 2))'),
    ('(X.T*A*X).diff(X)',
     'ArrayAdd(PermuteDims(ArrayTensorProduct(A*X, I), (1 2 3)), PermuteDims(ArrayTensorProduct(A.T*X, I), (3)(1 2)))'),
    ('(W.T*W).diff(W)',
     'ArrayAdd(PermuteDims(ArrayTensorProduct(W, I), (1 2 3)), PermuteDims(ArrayTensorProduct(W, I), (3)(1 2)))'),
    ('(W*W.T).diff(W)',
     'ArrayAdd(PermuteDims(ArrayTensorProduct(I, W.T), (3)(1 2)), PermuteDims(ArrayTensorProduct(I, W.T), (1 2 3)))'),
    ('(P.T*W).diff(W)',
     'PermuteDims(ArrayTensorProduct(P, I), (3)(1 2))'),
    ('(x*a.T).diff(x)',
     'PermuteDims(ArrayTensorProduct(I, a.T), (3)(1 2))'),
    ('(a*x.T).diff(x)',
     'PermuteDims(ArrayTensorProduct(I, a.T), (1 2 3))'),
    ('(X*x*a.T).diff(X)',
     'PermuteDims(ArrayTensorProduct(I, x*a.T), (3)(1 2))'),
    ('(-X).diff(X)',
     'PermuteDims(ArrayTensorProduct(-I, I), (3)(1 2))'),
    ('(X - A).diff(X)',
     'PermuteDims(ArrayTensorProduct(I, I), (3)(1 2))'),
    ('(Trace(X)*X).diff(X)',
     'ArrayAdd(ArrayTensorProduct(I, X), PermuteDims(ArrayTensorProduct(Trace(X)*I, I), (3)(1 2)))'),
    ('(A*X*Trace(B*X)).diff(X)',
     'ArrayAdd(ArrayTensorProduct(B.T, A*X), PermuteDims(ArrayTensorProduct(Trace(B*X)*A.T, I), (3)(1 2)))'),
    ('((X*A).T*B).diff(X)',
     'PermuteDims(ArrayTensorProduct(B, A), (1 2 3))'),
    ('(X*A).T.diff(X)',
     'PermuteDims(ArrayTensorProduct(I, A), (1 2 3))'),
    ('(M*N*M).diff(N)',
     'PermuteDims(ArrayTensorProduct(M.T, M), (3)(1 2))'),
    ('(M.T*N.T).diff(N)',
     'PermuteDims(ArrayTensorProduct(I, M), (1 2 3))'),
    ('(N**2).diff(N)',
     'ArrayAdd(PermuteDims(ArrayTensorProduct(I, N), (3)(1 2)), PermuteDims(ArrayTensorProduct(N.T, I), (3)(1 2)))'),
    ('(N**-1).diff(N)',
     'PermuteDims(ArrayTensorProduct(-N.T**(-1), N**(-1)), (3)(1 2))'),
    ('(R*X).diff(X)',
     'PermuteDims(ArrayTensorProduct(R.T, I), (3)(1 2))'),
]

# symox and SymPy 1.14 differ (term order of products of three or more
# factors, place of a scalar coefficient, SymPy RecursionError for X**-2):
# recorded gaps; the native and Python paths must still agree.
OTHER = [
    '(X**3).diff(X)',
    '(X**4).diff(X)',
    '(A*X - X.T*B).diff(X)',
    '(X**-2).diff(X)',
    '(Trace(X)**2*A).diff(X)',
    '(Trace(X*A)*X).diff(X)',
    '(X*A*X*B*X).diff(X)',
    '(A*X**3).diff(X)',
    '(3*X.T + 2*X).diff(X)',
]


def _python_path(expr, x, monkeypatch):
    with monkeypatch.context() as m:
        m.setattr(_core, '_matrix_lines_rank_four', lambda e, v: None)
        return md.matrix_derivative(expr, x)


@pytest.mark.parametrize('case,expected', SYMPY)
def test_rank_four_sympy(case, expected, monkeypatch):
    expr, var = eval(case.rsplit('.diff(', 1)[0]), eval(case.rsplit('.diff(', 1)[1][:-1])
    assert _core._matrix_lines_rank_four(expr, var) is not None
    result = eval(case)
    assert str(result) == expected
    assert result == _python_path(expr, var, monkeypatch)
    assert str(_python_path(expr, var, monkeypatch)) == expected


@pytest.mark.parametrize('case', OTHER)
def test_rank_four_same_as_python(case, monkeypatch):
    expr, var = eval(case.rsplit('.diff(', 1)[0]), eval(case.rsplit('.diff(', 1)[1][:-1])
    result = eval(case)
    assert str(result) == str(_python_path(expr, var, monkeypatch))
    assert type(result) is type(_python_path(expr, var, monkeypatch))


def test_rank_four_types():
    d = (A*X).diff(X)
    assert type(d) is PermuteDims and type(d.expr) is ArrayTensorProduct
    assert d.permutation.array_form == [0, 2, 1, 3]
    assert d.shape == (3, 3, 3, 3)
    assert type((X*X).diff(X)) is ArrayAdd
    assert type((Trace(X)*A).diff(X)) is ArrayTensorProduct
    # Explicit entries agree with the element-wise derivative.
    e = d.as_explicit()
    assert e[0, 1, 1, 1] == A[1, 0] and e[0, 1, 0, 1] == A[0, 0] and e[0, 1, 0, 0] == 0


@pytest.mark.parametrize('expr,var', [
    (lambda: (A*X).T*X, lambda: X),          # 1 of 2 lines held: same as Python
    (lambda: A*X + B, lambda: X),            # a constant addend: _rank_four refuses
    (lambda: Trace(A*X), lambda: X),         # rank two
    (lambda: (x.T*x), lambda: x),            # rank zero
    (lambda: A*X, lambda: X.T),              # not a MatrixSymbol variable
    (lambda: A*X**t, lambda: X),             # symbolic power
])
def test_rank_four_declines_like_python(expr, var, monkeypatch):
    e, v = expr(), var()
    native = _core._matrix_lines_rank_four(e, v)
    try:
        python = _python_path(e, v, monkeypatch)
    except Exception as exc:
        python = exc
    if native is None:
        return
    assert not isinstance(python, Exception)
    assert str(e.diff(v)) == str(python)


def test_rank_four_rejects_non_matrix_arguments():
    assert _core._matrix_lines_rank_four(A*X, 3) is None
    assert _core._matrix_lines_rank_four([1], X) is None
