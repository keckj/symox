"""Wave 41: Kronecker, matrix expression, trace and comprehension models are
computed natively, with SymPy 1.14's results and errors."""
import pytest
import sympy as sp
import symox as sx
from symox import _core


def _spy(monkeypatch, name):
    calls = []
    original = getattr(_core, name)
    def spy(*args, **kwargs):
        calls.append(args)
        return original(*args, **kwargs)
    monkeypatch.setattr(_core, name, spy)
    return calls


def test_combine_kronecker_is_native(monkeypatch):
    from symox.matrices.expressions.kronecker import combine_kronecker
    calls = _spy(monkeypatch, '_combine_kronecker')
    A, B, C, D = (sx.MatrixSymbol(name, 2, 2) for name in 'ABCD')
    assert combine_kronecker(sx.KroneckerProduct(A, B)*sx.KroneckerProduct(C, D)) == sx.KroneckerProduct(A*C, B*D)
    assert combine_kronecker(sx.KroneckerProduct(A, B)**3) == sx.KroneckerProduct(A**3, B**3)
    assert calls


def test_kronecker_models_match_sympy():
    A, B = sx.MatrixSymbol('A', 2, 3), sx.MatrixSymbol('B', 3, 2)
    W, Z = sx.MatrixSymbol('W', 3, 3), sx.MatrixSymbol('Z', 3, 3)
    assert (sx.KroneckerProduct(W, Z)*sx.KroneckerProduct(W.I, Z.I)).is_Identity
    assert sx.KroneckerProduct(sx.Identity(2), sx.Identity(3)) == sx.Identity(6)
    assert sx.KroneckerProduct(sx.eye(2), sx.eye(3)) == sx.eye(6)
    X, Y = sx.MatrixSymbol('X', 2, 2), sx.MatrixSymbol('Y', 2, 2)
    expanded = sx.KroneckerProduct(X + Y, X + Y).expand(kroneckerproduct=True)
    assert str(expanded) == str(sp.KroneckerProduct(*[sp.MatrixSymbol('X', 2, 2) + sp.MatrixSymbol('Y', 2, 2)]*2)
                                .expand(kroneckerproduct=True))
    i, j, m, n, o, p = sx.symbols('i j m n o p')
    entry = sx.KroneckerProduct(sx.MatrixSymbol('A', n, m), sx.MatrixSymbol('B', o, p))._entry(i, j)
    assert str(entry) == 'A[Mod(floor(i/o), n), Mod(floor(j/p), m)]*B[Mod(i, o), Mod(j, p)]'
    assert str(sx.KroneckerProduct(A, B).conjugate()) == str(
        sp.KroneckerProduct(sp.MatrixSymbol('A', 2, 3), sp.MatrixSymbol('B', 3, 2)).conjugate())


def test_matrix_expression_errors_match_sympy():
    n = sx.Symbol('n')
    A, B = sx.MatrixSymbol('A', n, n), sx.MatrixSymbol('B', 2, 2)
    with pytest.raises(ValueError, match='symbolic shape'):A.as_explicit()
    with pytest.raises(sx.NonSquareMatrixError, match='Power of non-square matrix'):
        sx.MatrixSymbol('C', 3, 2)**2
    with pytest.raises(sx.NonSquareMatrixError, match='Inverse of non-square matrix'):
        sx.MatrixSymbol('C', 3, 2).I
    with pytest.raises(NotImplementedError):2/B
    assert sx.Eq(sx.MatrixSymbol('P', 3, 1), sx.MatrixSymbol('Q', 2, 1)) == sx.false
    assert sx.Ne(sx.MatrixSymbol('P', 3, 1), sx.MatrixSymbol('Q', 2, 1)) == sx.true
    # Only a proven dimension mismatch decides the equality.
    m = sx.Symbol('m')
    assert isinstance(sx.Eq(sx.MatrixSymbol('P', n, 1), sx.MatrixSymbol('Q', m, 1)), sx.Equality)
    with pytest.raises(IndexError):B[4]
    assert B[3] == B[1, 1]
    with pytest.raises(IndexError):sx.MatrixSymbol('C', 2, n)[1]


def test_matrix_element_and_trace_attributes_are_native():
    X = sx.MatrixSymbol('X', 3, 3)
    assert X[0, 1].indices == (0, 1) and X[0, 1].symbol == X
    t = sx.Trace(X)
    assert t.as_explicit() == X[0, 0] + X[1, 1] + X[2, 2]
    with pytest.raises(ValueError):sx.Trace(sx.MatrixSymbol('A', sx.Symbol('n'), sx.Symbol('n'))).as_explicit()
    A, B, C = (sx.MatrixSymbol(name, 3, 3) for name in 'ABC')
    assert sx.Trace(B*A*C)._normalize() == sx.Trace(A*C*B)
    assert sx.Trace(X).doit(deep=False) == sx.Trace(X)
    # Plain expressions keep the attribute protocol: no Trace methods.
    x = sx.Symbol('x')
    assert not hasattr(x + 1, 'as_explicit') and not hasattr(x, '_normalize')
    assert not hasattr(x, 'indices')


def test_array_comprehension_is_native(monkeypatch):
    calls = _spy(monkeypatch, '_array_comprehension_entries')
    i, j, k = sx.symbols('i j k')
    a = sx.ArrayComprehension(i*j, (i, 1, 3), (j, 2, 4))
    assert a.doit().tolist() == [[2, 3, 4], [4, 6, 8], [6, 9, 12]]
    assert calls
    assert len(sx.ArrayComprehension(i, (j, k + 1, k + 5))) == 5
    with pytest.raises(ValueError, match='Symbolic length'):len(sx.ArrayComprehension(i*j, (i, 1, 3), (j, 2, k + 4)))
    with pytest.raises(TypeError, match='Bounds should be'):sx.ArrayComprehension(i, (i, 1, k/2))
    with pytest.raises(ValueError, match='Lower bound'):sx.ArrayComprehension(i, (i, 2, 1))
    with pytest.raises(ValueError, match='Variable should not'):sx.ArrayComprehension(i, (i, 1, i + 1))
    assert str(a) == 'ArrayComprehension(i*j, (i, 1, 3), (j, 2, 4))'
    assert a.kind is sx.UndefinedKind


def test_array_comprehension_map_matches_sympy():
    i, k = sx.symbols('i k')
    expr = sx.ArrayComprehensionMap(lambda i: i + 1, (i, 1, k))
    assert expr.doit() == expr
    assert expr.subs(k, 4) == sx.ArrayComprehensionMap(lambda i: i + 1, (i, 1, 4))
    assert expr.subs(k, 4).doit() == sx.ImmutableDenseNDimArray([2, 3, 4, 5])
    assert str(expr) == 'ArrayComprehensionMap((i, 1, k))'
    assert expr.free_symbols == {k}
    with pytest.raises(ValueError):sx.ArrayComprehensionMap(lambda i, j: i + j, (i, 1, 5)).doit()
    with pytest.raises(ValueError, match='Data type not supported'):sx.ArrayComprehensionMap(i, (i, 1, 2))


def test_explicit_matrix_symbolic_index_matches_sympy():
    x, y, z, w, i, j = sx.symbols('x y z w i j')
    M, Ms = sx.Matrix([[x, y], [z, w]]), sx.SparseMatrix([[2, 3], [4, 5]])
    assert isinstance(M[i, j], sx.MatrixElement) and isinstance(Ms[i, j], sx.MatrixElement)
    for a in (0, 1):
        for b in (0, 1):
            assert M[i, j].subs({i: a, j: b}) == M[a, b]
            assert Ms[i, j].subs({i: a, j: b}) == Ms[a, b]
    A = sx.MatrixSymbol('A', 2, 2)
    assert A[i, j].subs(A, M) == M[i, j] and M[i, j].subs(M, A) == A[i, j]
    assert M[i, j].diff(x) == sx.Matrix([[1, 0], [0, 0]])[i, j]
    # SymPy's `(j < 0) is True` test only fires for a host int index.
    for key in [(i, 2), (i, -1), (2, i), (-1, i)]:
        with pytest.raises(ValueError, match='index out of boundary'):M[key]
    held = M[sx.Integer(-1), i]
    assert isinstance(held, sx.MatrixElement) and held.args[1:] == (-1, i)
    assert str(M[sx.Symbol('k', negative=True), 0]) == str(
        sp.Matrix([[*sp.symbols('x y')], [*sp.symbols('z w')]])[sp.Symbol('k', negative=True), 0])
    with pytest.raises(IndexError, match='indices out of range'):M[sx.Integer(-3), i]
