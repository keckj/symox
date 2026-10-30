"""SymPy 1.14 parity of the native array-expression canonical forms
(crates/symox-core/src/arrayexpr.rs) and of convert_matrix_to_array
(crates/symox-core/src/arrayexpr_matrix.rs). The expected values are the
outputs of SymPy 1.14 on the same cases."""
import pytest
import symox as m
from symox.tensor.array.expressions import array_expressions as ae
from symox.tensor.array.expressions.from_matrix_to_array import convert_matrix_to_array as cm

MatrixSymbol, symbols = m.MatrixSymbol, m.symbols
ArraySymbol, ZeroArray = ae.ArraySymbol, ae.ZeroArray
tp, con, dia, pd, add, nest = ae._array_tensor_product, ae._array_contraction, ae._array_diagonal, ae._permute_dims, ae._array_add, ae.nest_permutation
PermuteDims, ArrayContraction, ArrayTensorProduct, ArrayDiagonal = ae.PermuteDims, ae.ArrayContraction, ae.ArrayTensorProduct, ae.ArrayDiagonal
k, x = symbols('k x')
M, N, P, Q = [MatrixSymbol(s, k, k) for s in 'MNPQ']
A = ArraySymbol('A', (k, k, k)); B = ArraySymbol('B', (k, k, k, k)); C = ArraySymbol('C', (k, k))
a = MatrixSymbol('a', k, 1); b = MatrixSymbol('b', k, 1)

CANON = {
 'tp_perm': lambda: tp(pd(M, [1, 0]), N),
 'tp_perm2': lambda: tp(M, pd(N, [1, 0]), pd(P, [1, 0])),
 'tp_contr': lambda: tp(Q, con(tp(M, N), (1, 2))),
 'tp_diag': lambda: tp(M, dia(tp(N, P), (1, 2))),
 'tp_zero': lambda: tp(M, ZeroArray(k, 2)),
 'con_nested': lambda: con(con(tp(M, N, P), (1, 2)), (0, 1)),
 'con_perm': lambda: con(pd(tp(M, N), [1, 0, 3, 2]), (1, 2)),
 'con_full_sort': lambda: con(tp(N, M, P), (0, 1), (2, 3)),
 'con_add': lambda: con(tp(add(M, N), P), (0, 1)),
 'con_diag': lambda: con(dia(tp(M, N, P), (1, 2)), (0, 3)),
 'con_trivial': lambda: con(tp(a, b), (1, 3)),
 'pd_nested': lambda: pd(pd(tp(M, N), [1, 0, 3, 2]), [2, 3, 0, 1]),
 'pd_tp_sort': lambda: pd(tp(M, N), [2, 3, 0, 1]),
 'pd_contr': lambda: pd(con(tp(M, N, P), (1, 2)), [1, 0, 3, 2]),
 'pd_contr2': lambda: pd(con(tp(M, N, P), (1, 2)), [3, 2, 1, 0]),
 'pd_zero': lambda: pd(ZeroArray(1, 2, 3), [2, 0, 1]),
 'dia_nested': lambda: dia(dia(tp(M, N, P), (0, 1)), (0, 1)),
 'dia_perm': lambda: dia(pd(tp(M, N), [1, 2, 3, 0]), (0, 1)),
 'dia_add': lambda: dia(add(tp(M, N), tp(N, P)), (1, 2)),
 'dia_trivial': lambda: dia(tp(M, N), (1,), (2, 3), allow_trivial_diags=True),
 'add_flat': lambda: add(add(M, N), ZeroArray(k, k), P),
 'add_zero': lambda: add(ZeroArray(k, k), ZeroArray(k, k)),
 'arr_mix': lambda: con(tp(A, B, C), (0, 3), (1, 7)),
 'arr_perm': lambda: pd(con(tp(A, C, B), (2, 3)), [4, 3, 2, 1, 0]),
 'nest_tp': lambda: nest(PermuteDims(tp(M, N, P), [1, 0, 2, 3, 5, 4])),
 'nest_tp_cross': lambda: nest(PermuteDims(tp(M, N), [2, 1, 0, 3])),
 'nest_con': lambda: nest(PermuteDims(ArrayContraction(ArrayTensorProduct(M, N, P), (1, 2)), [1, 0, 3, 2])),
 'doit1': lambda: ArrayContraction(ArrayTensorProduct(PermuteDims(M, [1, 0]), ArrayContraction(ArrayTensorProduct(N, P), (1, 2))), (0, 1)).doit(),
 'doit2': lambda: PermuteDims(PermuteDims(ArrayTensorProduct(M, N), [1, 0, 3, 2]), [1, 0, 3, 2]).doit(),
 'doit3': lambda: ArrayDiagonal(ArrayDiagonal(ArrayTensorProduct(M, N, P), (0, 1)), (0, 1)).doit(),
 'err_con': lambda: con(tp(M, a), (0, 3)),
 'err_dia': lambda: dia(tp(M, a), (0, 3)),
 'err_dia1': lambda: dia(M, (0,)),
 'err_perm': lambda: pd(M, [0, 1, 2]),
}

CANON_EXPECTED = {'tp_perm': 'PermuteDims(ArrayTensorProduct(M, N), [1, 0, 2, 3])',
    'tp_perm2': 'PermuteDims(ArrayTensorProduct(M, N, P), [0, 1, 3, 2, 5, 4])',
    'tp_contr': 'ArrayContraction(ArrayTensorProduct(Q, M, N), (3, 4))',
    'tp_diag': 'ArrayDiagonal(ArrayTensorProduct(M, N, P), (3, 4))',
    'tp_zero': 'ZeroArray(k, k, k, 2)',
    'con_nested': 'ArrayContraction(ArrayTensorProduct(M, N, P), (0, 3), (1, 2))',
    'con_perm': 'ArrayContraction(ArrayTensorProduct(M, N), (0, 3))',
    'con_full_sort': 'ArrayContraction(ArrayTensorProduct(M, N, P), (0, 1), (2, 3))',
    'con_add': 'ArrayTensorProduct(ArrayAdd(ArrayContraction(M, (0, 1)), ArrayContraction(N, (0, 1))), P)',
    'con_diag': 'ArrayDiagonal(ArrayContraction(ArrayTensorProduct(M, N, P), (0, 5)), (0, 1))',
    'con_trivial': 'ArrayContraction(ArrayTensorProduct(a, b), (1, 3))',
    'pd_nested': 'PermuteDims(ArrayTensorProduct(N, M), [1, 0, 3, 2])',
    'pd_tp_sort': 'ArrayTensorProduct(N, M)',
    'pd_contr': 'PermuteDims(ArrayContraction(ArrayTensorProduct(N, M, P), (0, 3)), [0, 1, 3, 2])',
    'pd_contr2': 'PermuteDims(ArrayContraction(ArrayTensorProduct(P, N, M), (2, 5)), [1, 0, 2, 3])',
    'pd_zero': 'ZeroArray(3, 1, 2)',
    'dia_nested': 'ArrayDiagonal(ArrayTensorProduct(M, N, P), (0, 1), (2, 3))',
    'dia_perm': 'PermuteDims(ArrayDiagonal(ArrayTensorProduct(M, N), (1, 2)), [1, 0, 2])',
    'dia_add': 'ArrayAdd(ArrayDiagonal(ArrayTensorProduct(M, N), (1, 2)), ArrayDiagonal(ArrayTensorProduct(N, P), (1, 2)))',
    'dia_trivial': 'ArrayDiagonal(ArrayTensorProduct(M, N), (2, 3))',
    'add_flat': 'ArrayAdd(M, N, P)',
    'add_zero': 'ZeroArray(k, k)',
    'arr_mix': 'ArrayContraction(ArrayTensorProduct(A, B, C), (0, 3), (1, 7))',
    'arr_perm': 'ERR ValueError',
    'nest_tp': 'PermuteDims(ArrayTensorProduct(M, N, P), [1, 0, 2, 3, 5, 4])',
    'nest_tp_cross': 'PermuteDims(ArrayTensorProduct(N, M), [0, 3, 2, 1])',
    'nest_con': 'PermuteDims(ArrayContraction(ArrayTensorProduct(N, M, P), (0, 3)), [0, 1, 3, 2])',
    'doit1': 'ArrayContraction(ArrayTensorProduct(M, N, P), (0, 1), (3, 4))',
    'doit2': 'ArrayTensorProduct(M, N)',
    'doit3': 'ArrayDiagonal(ArrayTensorProduct(M, N, P), (0, 1), (2, 3))',
    'err_con': 'ERR ValueError',
    'err_dia': 'ERR ValueError',
    'err_dia1': 'ERR ValueError',
    'err_perm': 'ERR ValueError'}

MATRIX = {
 'mul': lambda: M*N, 'mul3': lambda: M*N*P, 'smul': lambda: 2*M*N, 'xmul': lambda: 2*x*M*N, 'neg': lambda: -M,
 'add': lambda: M + N, 'add2': lambda: M*N + N.T, 'tr': lambda: M.T, 'trmul': lambda: (M*N).T,
 'trace': lambda: m.Trace(M), 'trace2': lambda: m.Trace(M*N), 'trpow': lambda: m.Trace(M)**2, 'trmulx': lambda: x*m.Trace(M)*m.Trace(N),
 'pow2': lambda: M**2, 'pow3': lambda: M**3, 'powk': lambda: M**k, 'inv': lambda: M**-1, 'powmn': lambda: (M*N)**2, 'powT': lambda: (M.T)**2,
 'had': lambda: m.hadamard_product(M, N), 'had3': lambda: m.hadamard_product(M, N, P), 'hp2': lambda: m.HadamardPower(M, 2), 'hpk': lambda: m.HadamardPower(M, k),
 'kron': lambda: m.KroneckerProduct(M, N), 'vec': lambda: a.T*M*a, 'zero': lambda: m.ZeroMatrix(k, k) + M, 'ident': lambda: m.Identity(k)*M,
 'scalar': lambda: x**2, 'scalar2': lambda: 2*x, 'atrace': lambda: m.Trace(M)*M,
}

MATRIX_EXPECTED = {'mul': 'ArrayContraction(ArrayTensorProduct(MatrixSymbol:M, MatrixSymbol:N), (1, 2))',
    'mul3': 'ArrayContraction(ArrayTensorProduct(MatrixSymbol:N, MatrixSymbol:M, MatrixSymbol:P), (0, 3), (1, 4))',
    'smul': 'ArrayContraction(ArrayTensorProduct(Integer:2, MatrixSymbol:M, MatrixSymbol:N), (1, 2))',
    'xmul': 'ArrayContraction(ArrayTensorProduct(Integer:2, Symbol:x, MatrixSymbol:M, MatrixSymbol:N), (1, 2))',
    'neg': 'ArrayTensorProduct(NegativeOne:-1, MatrixSymbol:M)',
    'add': 'ArrayAdd(MatrixSymbol:M, MatrixSymbol:N)',
    'add2': 'ArrayAdd(ArrayContraction(ArrayTensorProduct(MatrixSymbol:M, MatrixSymbol:N), (1, 2)), PermuteDims(MatrixSymbol:N, [1, 0]))',
    'tr': 'PermuteDims(MatrixSymbol:M, [1, 0])',
    'trmul': 'ArrayContraction(ArrayTensorProduct(MatrixSymbol:N, MatrixSymbol:M), (0, 3))',
    'trace': 'ArrayContraction(MatrixSymbol:M, (0, 1))',
    'trace2': 'ArrayContraction(ArrayTensorProduct(MatrixSymbol:M, MatrixSymbol:N), (0, 3), (1, 2))',
    'trpow': 'ArrayContraction(ArrayTensorProduct(MatrixSymbol:M, MatrixSymbol:M), (0, 1), (2, 3))',
    'trmulx': 'ArrayContraction(ArrayTensorProduct(Symbol:x, MatrixSymbol:M, MatrixSymbol:N), (0, 1), (2, 3))',
    'pow2': 'ArrayContraction(ArrayTensorProduct(MatrixSymbol:M, MatrixSymbol:M), (1, 2))',
    'pow3': 'ArrayContraction(ArrayTensorProduct(MatrixSymbol:M, MatrixSymbol:M, MatrixSymbol:M), (0, 3), (1, 4))',
    'powk': 'Apply(MatrixSymbol:M, d**k)',
    'inv': 'Inverse:M**(-1)',
    'powmn': "ERR AttributeError 'ArrayContraction' object has no attribu", 'powT': "ERR AttributeError 'PermuteDims' object has no attribute 'a", 'had': 'ArrayDiagonal(ArrayTensorProduct(MatrixSymbol:M, MatrixSymbol:N), (0, 2), (1, 3))',
    'had3': 'ArrayDiagonal(ArrayTensorProduct(MatrixSymbol:M, MatrixSymbol:N, MatrixSymbol:P), (0, 2, 4), (1, 3, 5))',
    'hp2': 'ArrayDiagonal(ArrayTensorProduct(MatrixSymbol:M, MatrixSymbol:M), (0, 2), (1, 3))',
    'hpk': 'Apply(MatrixSymbol:M, d**k)',
    'kron': 'Reshape(PermuteDims(ArrayTensorProduct(MatrixSymbol:M, MatrixSymbol:N), [0, 2, 1, 3]), (k**2, k**2))',
    'vec': 'ArrayContraction(ArrayTensorProduct(MatrixSymbol:M, MatrixSymbol:a, MatrixSymbol:a), (0, 2), (1, 4))',
    'zero': 'MatrixSymbol:M',
    'ident': 'MatrixSymbol:M',
    'scalar': 'ArrayTensorProduct(Symbol:x, Symbol:x)',
    'scalar2': 'ArrayTensorProduct(Integer:2, Symbol:x)',
    'atrace': 'ArrayContraction(ArrayTensorProduct(MatrixSymbol:M, MatrixSymbol:M), (0, 1))'}


def ser_canon(e):
    n = type(e).__name__
    if n == 'PermuteDims': return 'PermuteDims(%s, %s)' % (ser_canon(e.expr), list(e.permutation.array_form))
    if n in ('ArrayTensorProduct', 'ArrayAdd'): return '%s(%s)' % (n, ', '.join(ser_canon(a) for a in e.args))
    if n in ('ArrayContraction', 'ArrayDiagonal'): return '%s(%s, %s)' % (n, ser_canon(e.expr), ', '.join(str(tuple(g)) for g in e.args[1:]))
    return repr(e)

def ser_matrix(e):
    n = type(e).__name__
    if n == 'PermuteDims': return 'PermuteDims(%s, %s)' % (ser_matrix(e.expr), list(e.permutation.array_form))
    if n in ('ArrayTensorProduct', 'ArrayAdd'): return '%s(%s)' % (n, ', '.join(ser_matrix(a) for a in e.args))
    if n in ('ArrayContraction', 'ArrayDiagonal'): return '%s(%s, %s)' % (n, ser_matrix(e.expr), ', '.join(str(tuple(g)) for g in e.args[1:]))
    if n == 'Reshape': return 'Reshape(%s, %s)' % (ser_matrix(e.expr), tuple(e.shape))
    if n == 'ArrayElementwiseApplyFunc': return 'Apply(%s, %s)' % (ser_matrix(e.expr), e.function(m.Symbol('d')))
    return '%s:%s' % (n, e)


@pytest.mark.parametrize('name', sorted(CANON))
def test_canonical_forms_match_sympy(name):
    try:
        got = ser_canon(CANON[name]())
    except Exception as e:
        got = 'ERR ' + type(e).__name__
    assert got == CANON_EXPECTED[name]


@pytest.mark.parametrize('name', sorted(MATRIX))
def test_convert_matrix_to_array_matches_sympy(name):
    try:
        got = ser_matrix(cm(MATRIX[name]()))
    except Exception as e:
        got = 'ERR ' + type(e).__name__ + ' ' + str(e)[:40]
    assert got == MATRIX_EXPECTED[name]


def test_hadamard_power_symbolic_exponent_dummy_eq():
    d0 = m.Dummy('d0')
    result = ae.ArrayElementwiseApplyFunc(m.Lambda(d0, d0**k), M)
    assert cm(m.HadamardPower(M, k)).dummy_eq(result)


def _groups_digest(e):
    import hashlib
    return hashlib.md5(str([tuple(int(i) for i in g) for g in e.args[1:]]).encode()).hexdigest()[:12]


# SymPy 1.14 expands every positive integer power (no size limit); the
# digests are those of SymPy's contraction/diagonal groups for n = 5000.
@pytest.mark.parametrize('name,build,kind,groups,digest', [
    ('matpow', lambda n: M**n, 'ArrayContraction', 4999, '0bdba11c768b'),
    ('tracepow', lambda n: m.Trace(M)**n, 'ArrayContraction', 5000, '3972054d6935'),
    ('hadpow', lambda n: m.HadamardPower(M, n), 'ArrayDiagonal', 2, 'aa042e083334'),
])
def test_convert_large_integer_powers_match_sympy(name, build, kind, groups, digest):
    out = cm(build(5000))
    assert type(out).__name__ == kind
    assert type(out.expr).__name__ == 'ArrayTensorProduct'
    assert len(out.expr.args) == 5000 and set(out.expr.args) == {M}
    assert len(out.args) - 1 == groups
    assert _groups_digest(out) == digest


def test_convert_large_scalar_power_matches_sympy():
    out = cm(x**5000)
    assert type(out).__name__ == 'ArrayTensorProduct'
    assert len(out.args) == 5000 and set(out.args) == {x}


def test_convert_unbuildable_power_is_a_memory_error():
    # SymPy's range loop cannot build 10**30 factors either
    with pytest.raises(MemoryError):
        cm(M**(10**30))
