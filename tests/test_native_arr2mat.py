"""SymPy 1.14 parity of the native convert_array_to_matrix, its passes and
the _EditArrayContraction editor (crates/symox-core/src/arrayexpr_a2m.rs).

The expected values are SymPy 1.14's outputs. CASES are taken from the
differential corpus of audits/performance/2026-10-30-arr2mat (inputs that
the plain constructors rebuild identically); CURATED are hand-picked calls.
Each case also checks that the stock fallback (symox._arrayexpr) is not
used, so the result comes from the Rust port."""
import pytest
import symox as S
from symox.tensor.array.expressions import _stock
from symox.tensor.array.expressions import array_expressions as ae
from symox.tensor.array.expressions import from_array_to_matrix as a2m
from symox.tensor.array.expressions.array_expressions import (ArrayAdd, ArrayContraction, ArrayDiagonal,
                                                              ArraySymbol, ArrayTensorProduct, OneArray)
from symox.tensor.array.expressions.from_matrix_to_array import convert_matrix_to_array as cm

ca = a2m.convert_array_to_matrix
tp, con, dia, pd = ae._array_tensor_product, ae._array_contraction, ae._array_diagonal, ae._permute_dims
k = S.Symbol('k')
x, y = S.symbols('x y')
M, N, P = [S.MatrixSymbol(s, k, k) for s in 'MNP']
a, b = S.MatrixSymbol('a', k, 1), S.MatrixSymbol('b', k, 1)
c = S.MatrixSymbol('c', 1, k)
X, Y = S.MatrixSymbol('X', 3, 3), S.MatrixSymbol('Y', 3, 3)
u = S.MatrixSymbol('u', 3, 1)
I, I3 = S.Identity(k), S.Identity(3)
A = ArraySymbol('A', (k, k))
NS = dict(ArrayContraction=ArrayContraction, ArrayDiagonal=ArrayDiagonal, ArrayTensorProduct=ArrayTensorProduct,
          ArrayAdd=ArrayAdd, OneArray=OneArray, k=k, x=x, y=y, M=M, N=N, P=P, a=a, b=b, c=c, X=X, Y=Y, u=u,
          I=I, I3=I3, A=A)


def editor(e):
    ed = ae._EditArrayContraction(e)
    out = [ed.number_of_contraction_indices, [str(t.element) + str(t.indices) for t in ed.args_with_ind],
           ed.get_contraction_indices(), str(ed.get_contraction_indices_to_ind_rel_pos()),
           ed.number_of_diagonal_indices]
    ed.merge_scalars()
    out.append(str(ed.to_array_contraction()))
    return out


PASSES = {
    'convert': ca,
    'a2m': a2m._array2matrix,
    'rtd': a2m._remove_trivial_dims,
    'hadamard': a2m.identify_hadamard_products,
    'removable': a2m.identify_removable_identity_matrices,
    'remove_id': a2m.remove_identity_matrices,
    'multi_id': a2m._array_contraction_to_diagonal_multiple_identity,
    'diag2contr': a2m._array_diag2contr_diagmatrix,
    'editor': editor,
}


def fmt(v):
    if isinstance(v, (list, tuple)):
        return '[' + ', '.join(fmt(i) for i in v) + ']'
    return str(v)


def outcome(fn):
    try:
        return fmt(fn())
    except Exception as e:
        return 'EXC:' + type(e).__name__


@pytest.fixture
def no_stock(monkeypatch):
    def refuse(*args, **kwargs):
        raise AssertionError('stock fallback used')
    monkeypatch.setattr(_stock, 'to_stock', refuse)


CASES = {
    'convert': [
        ('ArrayDiagonal(ArrayTensorProduct(Y, I, c, N), (0, 1), (3, 6), (2, 5))',
         'PermuteDims(ArrayDiagonal(ArrayTensorProduct(Y, I, c, N), (0, 1), (2, 5), (3, 6)), (3 4))'),
        ('ArrayContraction(N, (0, 1))',
         'Trace(N)'),
        ('ArrayContraction(ArrayTensorProduct(x, N, OneArray(k), c, c), (0, 2), (1, 6))',
         'ArrayContraction(x*c.T*c*N.T, (1,))'),
        ('ArrayDiagonal(ArrayTensorProduct(P, N), (2, 3), (0, 1))',
         'PermuteDims(ArrayDiagonal(ArrayTensorProduct(P, N), (0, 1), (2, 3)), (0 1))'),
        ('ArrayContraction(ArrayTensorProduct(a, b, u), (0, 2), (1, 3))',
         'Trace(a*b.T)*u'),
        ('ArrayAdd(ArrayContraction(ArrayTensorProduct(N, P, Y, b), (0, 3), (1, 2), (4, 5)), ArrayContraction(ArrayTensorProduct(N, P, Y, b), (0, 3), (1, 2), (4, 5)))',
         '(2*Trace(N*P)*Trace(Y))*b'),
        ('ArrayContraction(ArrayTensorProduct(P, X), (0, 1), (2, 3))',
         'Trace(P)*Trace(X)'),
    ],
    'a2m': [
        ('ArrayDiagonal(ArrayTensorProduct(Y, I, c, N), (0, 1), (3, 6), (2, 5))',
         'PermuteDims(ArrayDiagonal(ArrayTensorProduct(Y, I, c, N), (0, 1), (2, 5), (3, 6)), (3 4))'),
        ('ArrayContraction(N, (0, 1))',
         'Trace(N)'),
        ('ArrayContraction(ArrayTensorProduct(x, N, OneArray(k), c, c), (0, 2), (1, 6))',
         'ArrayContraction(ArrayTensorProduct(OneArray(k), c, x*N*c.T), (0, 3))'),
        ('ArrayDiagonal(ArrayTensorProduct(P, N), (2, 3), (0, 1))',
         'PermuteDims(ArrayDiagonal(ArrayTensorProduct(P, N), (0, 1), (2, 3)), (0 1))'),
        ('ArrayContraction(ArrayTensorProduct(a, b, u), (0, 2), (1, 3))',
         'Trace(a*b.T)*u'),
        ('ArrayAdd(ArrayContraction(ArrayTensorProduct(N, P, Y, b), (0, 3), (1, 2), (4, 5)), ArrayContraction(ArrayTensorProduct(N, P, Y, b), (0, 3), (1, 2), (4, 5)))',
         '(2*Trace(N*P)*Trace(Y))*b'),
        ('ArrayContraction(ArrayTensorProduct(P, X), (0, 1), (2, 3))',
         'Trace(P)*Trace(X)'),
    ],
    'hadamard': [
        ('ArrayDiagonal(ArrayTensorProduct(Y, I, c, N), (0, 1), (3, 6), (2, 5))',
         'PermuteDims(ArrayDiagonal(ArrayTensorProduct(Y, I, c, N), (0, 1), (2, 5), (3, 6)), (3 4))'),
        ('ArrayContraction(ArrayTensorProduct(x, N, OneArray(k), c, c), (0, 2), (1, 6))',
         'ArrayContraction(ArrayTensorProduct(x*N, OneArray(k), c, c), (0, 2), (1, 6))'),
        ('ArrayDiagonal(ArrayTensorProduct(P, N), (2, 3), (0, 1))',
         'PermuteDims(ArrayDiagonal(ArrayTensorProduct(P, N), (0, 1), (2, 3)), (0 1))'),
        ('ArrayContraction(ArrayTensorProduct(a, b, u), (0, 2), (1, 3))',
         'Trace(a*b.T)*u'),
        ('ArrayAdd(ArrayDiagonal(ArrayTensorProduct(OneArray(1, k), c), (0, 2), (1, 3)), ArrayDiagonal(ArrayTensorProduct(OneArray(1, k), c), (0, 2), (1, 3)))',
         'EXC:NotImplementedError'),
        ('ArrayAdd(ArrayContraction(ArrayTensorProduct(x, P, Y), (0, 1), (2, 3)), ArrayContraction(ArrayTensorProduct(x, P, Y), (0, 1), (2, 3)))',
         'EXC:NotImplementedError'),
        ('ArrayDiagonal(ArrayTensorProduct(2, c, M), (1, 3))',
         'ArrayDiagonal(ArrayTensorProduct(2*c, M), (1, 3))'),
    ],
    'removable': [
        ('ArrayDiagonal(ArrayTensorProduct(Y, I, c, N), (0, 1), (3, 6), (2, 5))',
         'PermuteDims(ArrayDiagonal(ArrayTensorProduct(Y, I, c, N), (0, 1), (2, 5), (3, 6)), (3 4))'),
        ('ArrayContraction(ArrayTensorProduct(x, N, OneArray(k), c, c), (0, 2), (1, 6))',
         'ArrayContraction(ArrayTensorProduct(x*N, OneArray(k), c, c), (0, 2), (1, 6))'),
        ('ArrayDiagonal(ArrayTensorProduct(P, N), (2, 3), (0, 1))',
         'PermuteDims(ArrayDiagonal(ArrayTensorProduct(P, N), (0, 1), (2, 3)), (0 1))'),
        ('ArrayAdd(ArrayDiagonal(ArrayTensorProduct(OneArray(1, k), c), (0, 2), (1, 3)), ArrayDiagonal(ArrayTensorProduct(OneArray(1, k), c), (0, 2), (1, 3)))',
         'EXC:NotImplementedError'),
        ('ArrayAdd(ArrayContraction(ArrayTensorProduct(x, P, Y), (0, 1), (2, 3)), ArrayContraction(ArrayTensorProduct(x, P, Y), (0, 1), (2, 3)))',
         'EXC:NotImplementedError'),
        ('ArrayDiagonal(ArrayTensorProduct(2, c, M), (1, 3))',
         'ArrayDiagonal(ArrayTensorProduct(2*c, M), (1, 3))'),
        ('ArrayContraction(ArrayTensorProduct(2, c, I), (1, 2))',
         'ArrayContraction(ArrayTensorProduct(2*c, I), (1, 2))'),
    ],
    'diag2contr': [
        ('ArrayDiagonal(ArrayTensorProduct(Y, I, c, N), (0, 1), (3, 6), (2, 5))',
         'EXC:ValueError'),
        ('ArrayContraction(ArrayTensorProduct(x, N, OneArray(k), c, c), (0, 2), (1, 6))',
         'EXC:AttributeError'),
        ('ArrayDiagonal(ArrayTensorProduct(2, c, M), (1, 3))',
         'ArrayContraction(ArrayTensorProduct(2, OneArray(1), M, DiagMatrix(c)), (2, 4))'),
        ('ArrayDiagonal(ArrayTensorProduct(M, b), (0, 2))',
         'ArrayContraction(ArrayTensorProduct(M, OneArray(1), DiagMatrix(b)), (0, 3))'),
        ('ArrayDiagonal(ArrayTensorProduct(a, P), (0, 2))',
         'ArrayContraction(ArrayTensorProduct(OneArray(1), P, DiagMatrix(a)), (1, 3))'),
        ('ArrayDiagonal(ArrayTensorProduct(b, N, OneArray(k), A), (3, 5), (2, 4), (0, 6))',
         'ArrayDiagonal(ArrayContraction(ArrayTensorProduct(OneArray(1), N, OneArray(k), A, DiagMatrix(b)), (5, 6)), (3, 5), (2, 4))'),
        ('ArrayDiagonal(ArrayTensorProduct(x, c, P, c), (1, 2))',
         'ArrayContraction(ArrayTensorProduct(x, OneArray(1), P, c, DiagMatrix(c)), (1, 6))'),
    ],
    'rtd': [
        ('ArrayContraction(ArrayTensorProduct(x, N, OneArray(k), c, c), (0, 2), (1, 6))',
         '[ArrayContraction(x*c.T*c*N.T, (1,)), [0, 2]]'),
        ('ArrayTensorProduct(2, OneArray(1, k))',
         '[2, [0, 1]]'),
        ('ArrayDiagonal(ArrayTensorProduct(OneArray(1, k), N), (1, 3))',
         '[N, [0]]'),
        ('ArrayContraction(ArrayTensorProduct(P, P, OneArray(k), I), (0, 2), (1, 3), (4, 6))',
         '[ArrayContraction(ArrayTensorProduct(P, P), (0, 2), (1, 3)), [0]]'),
        ('ArrayDiagonal(ArrayTensorProduct(y, OneArray(k), M), (1, 2))',
         '[ArrayDiagonal(y*M, (0, 1)), [0]]'),
        ('ArrayContraction(ArrayTensorProduct(P, OneArray(1, k), P), (0, 4), (1, 3))',
         '[ArrayContraction(ArrayTensorProduct(P, P), (0, 2), (1,)), [0]]'),
        ('ArrayDiagonal(ArrayTensorProduct(2, OneArray(k), u, OneArray(k)), (0, 3))',
         '[2*u, [2]]'),
    ],
    'remove_id': [
        ('ArrayAdd(ArrayDiagonal(ArrayTensorProduct(OneArray(1, k), c), (0, 2), (1, 3)), ArrayDiagonal(ArrayTensorProduct(OneArray(1, k), c), (0, 2), (1, 3)))',
         'EXC:NotImplementedError'),
        ('ArrayAdd(ArrayContraction(ArrayTensorProduct(x, P, Y), (0, 1), (2, 3)), ArrayContraction(ArrayTensorProduct(x, P, Y), (0, 1), (2, 3)))',
         'EXC:NotImplementedError'),
    ],
    'multi_id': [
        ('ArrayDiagonal(ArrayTensorProduct(Y, I, c, N), (0, 1), (3, 6), (2, 5))',
         'EXC:ValueError'),
        ('ArrayDiagonal(N, (0, 1))',
         'EXC:ValueError'),
    ],
    'editor': [
        ('ArrayDiagonal(ArrayTensorProduct(Y, I, c, N), (0, 1), (3, 6), (2, 5))',
         'EXC:IndexError'),
        ('ArrayContraction(N, (0, 1))',
         '[1, [N[0, 0]], [[0, 1]], [[_IndPos(0, 0), _IndPos(0, 1)]], 0, ArrayContraction(N, (0, 1))]'),
        ('ArrayDiagonal(N, (0, 1))',
         'EXC:IndexError'),
        ('ArrayContraction(ArrayTensorProduct(x, N, OneArray(k), c, c), (0, 2), (1, 6))',
         '[2, [x[], N[0, 1], OneArray(k)[0], c[None, None], c[None, 1]], [[0, 2], [1, 6]], [[_IndPos(1, 0), _IndPos(2, 0)], [_IndPos(1, 1), _IndPos(4, 1)]], 0, ArrayContraction(ArrayTensorProduct(x*N, OneArray(k), c, c), (0, 2), (1, 6))]'),
        ('ArrayContraction(ArrayTensorProduct(A, M), (0, 3), (1, 2))',
         '[2, [A[0, 1], M[1, 0]], [[0, 3], [1, 2]], [[_IndPos(0, 0), _IndPos(1, 1)], [_IndPos(0, 1), _IndPos(1, 0)]], 0, ArrayContraction(ArrayTensorProduct(A, M), (0, 3), (1, 2))]'),
        ('ArrayContraction(ArrayTensorProduct(a, b, u), (0, 2), (1, 3))',
         '[2, [a[0, 1], b[0, 1], u[None, None]], [[0, 2], [1, 3]], [[_IndPos(0, 0), _IndPos(1, 0)], [_IndPos(0, 1), _IndPos(1, 1)]], 0, ArrayContraction(ArrayTensorProduct(a, b, u), (0, 2), (1, 3))]'),
        ('ArrayContraction(ArrayTensorProduct(P, X), (0, 1), (2, 3))',
         '[2, [P[0, 0], X[1, 1]], [[0, 1], [2, 3]], [[_IndPos(0, 0), _IndPos(0, 1)], [_IndPos(1, 0), _IndPos(1, 1)]], 0, ArrayContraction(ArrayTensorProduct(P, X), (0, 1), (2, 3))]'),
    ],
}

PARAMS = [(p, i, o) for p, lst in CASES.items() for i, o in lst]


@pytest.mark.parametrize('name,inp,expected', PARAMS, ids=['%s-%d' % (p[0], n) for n, p in enumerate(PARAMS)])
def test_corpus_cases_match_sympy(name, inp, expected, no_stock):
    e = eval(inp.replace('Identity(3)', 'I3'), dict(NS))
    assert str(e) == inp
    assert outcome(lambda: PASSES[name](e)) == expected


CURATED = [
    ('ca(cm(M*N))', 'M*N'),
    ('ca(cm((M*N).T))', 'N.T*M.T'),
    ('ca(cm(M*N*P))', 'M*N*P'),
    ('ca(cm(a.T*M*b))', 'a.T*M*b'),
    ('ca(cm(S.Trace(M*N)))', 'Trace(M*N)'),
    ('ca(cm(S.HadamardProduct(M, N)))', 'HadamardProduct(M, N)'),
    ('ca(cm(M*N + N.T*P))', 'M*N + N.T*P'),
    ('ca(tp(a, b))', 'a*b.T'),
    ('ca(dia(tp(M, a), (1, 2)))', 'M*DiagMatrix(a)'),
    ('ca(con(tp(M, I, N), (1, 2), (3, 4)))', 'M*N'),
    ('ca(con(tp(M, N, P), (0, 2), (1, 5)))', 'N.T*M*P.T'),
    ('ca(cm(2*S.Trace(N*P*P)*N*M))', '(2*Trace(P.T*N.T*P.T))*N*M'),
    ('ca(cm(S.Trace(N*M)*P))', 'Trace(M*N)*P'),
    ('ca(cm(x*M*N.T))', 'x*M*N.T'),
    ('ca(cm(S.HadamardPower(M, 2)))', 'HadamardPower(M, 2)'),
    ('ca(cm(S.KroneckerProduct(M, N)))', 'Reshape(PermuteDims(ArrayTensorProduct(M, N), (3)(1 2)), (k**2, k**2))'),
    ('ca(ae.ArrayElementwiseApplyFunc(S.sin, a.T*b))', 'sin(a.T*b)'),
    ('ca(ae.ArrayElementwiseApplyFunc(S.exp, M*N))', 'Lambda(_d, exp(_d)).(M*N)'),
    ('ca(pd(tp(M, N), [1, 0, 3, 2]))', 'ArrayTensorProduct(M.T, N.T)'),
    ('ca(ae.Reshape(M, (k*k,)))', 'Reshape(M, (k**2,))'),
    ('a2m._array_diag2contr_diagmatrix(dia(tp(M, a), (1, 2)))', 'ArrayContraction(ArrayTensorProduct(M, OneArray(1), DiagMatrix(a)), (1, 3))'),
    ('a2m._array_diag2contr_diagmatrix(con(tp(M, N), (1, 2)))', 'EXC:AttributeError'),
    ('a2m._array_diag2contr_diagmatrix(tp(M, N))', 'EXC:AttributeError'),
    ('a2m.identify_hadamard_products(dia(tp(M, N), (0, 2), (1, 3)))', 'HadamardProduct(M, N)'),
    ('a2m.remove_identity_matrices(con(tp(M, I, N), (1, 2), (3, 4)))', '[ArrayContraction(ArrayTensorProduct(I, M, N), (0, 3), (1, 4)), []]'),
    ('a2m._remove_trivial_dims(tp(a.T, b))', '[a*b.T, [0, 3]]'),
    ('a2m._support_function_tp1_recognize([(1, 2)], [M, N])', 'M*N'),
    ('a2m._support_function_tp1_recognize([(1, 2), (3, 4)], [M, I, N])', 'M*N'),
]


@pytest.mark.parametrize('call,expected', CURATED, ids=[c[0] for c in CURATED])
def test_curated_calls_match_sympy(call, expected, no_stock):
    def run():
        r = eval(call)
        return '[' + ', '.join(str(i) for i in r) + ']' if isinstance(r, (list, tuple)) else str(r)
    assert outcome(run) == expected


@pytest.mark.parametrize('expr,message', [
    (lambda: con(tp(M, N), (1, 2)), "'ArrayContraction' object has no attribute 'diagonal_indices'"),
    (lambda: tp(M, N), "'ArrayTensorProduct' object has no attribute 'expr'"),
    (lambda: M, "'MatrixSymbol' object has no attribute 'expr'"),
    (lambda: ArrayAdd(M, N), "'ArrayAdd' object has no attribute 'expr'"),
])
def test_diag2contr_attribute_errors_match_sympy(expr, message, no_stock):
    with pytest.raises(AttributeError) as info:
        a2m._array_diag2contr_diagmatrix(expr())
    assert str(info.value) == message


def test_element_wise_function_of_a_1x1_matrix_is_native(no_stock):
    # SymPy: expr.function(subexpr) for a 1x1 argument, removed axes [0, 1].
    f = (a.T * b).applyfunc(S.sin)
    assert fmt(a2m._remove_trivial_dims(f)) == '[sin(a.T*b), [0, 1]]'
    assert str(ca(f)) == 'sin(a.T*b)'
    # SymPy multiplies the scalar into the matrix: noncommutative scalar error.
    with pytest.raises(NotImplementedError, match='noncommutative scalars in MatMul are not supported.'):
        ca(tp(f, a))


def test_trace_of_matmul_keeps_sympy_order():
    # MatMul._eval_trace: no cyclic rotation of the factors.
    assert str(S.Trace(N * M) * P) == 'Trace(N*M)*P'
    assert str(ca(cm(S.Trace(N * M) * P))) == 'Trace(M*N)*P'
