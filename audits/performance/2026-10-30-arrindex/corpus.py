"""Differential corpus of array_derive, convert_indexed_to_array and
convert_array_to_indexed: the same cases evaluated with SymPy 1.14 and with
symox, serialized structurally (Dummy-free)."""
import importlib


def namespace(mod):
    m = importlib.import_module(mod)
    ae = importlib.import_module(mod + '.tensor.array.expressions')
    ns = dict(m=m)
    for n in ['ArraySymbol', 'ArrayTensorProduct', 'ArrayContraction', 'ArrayDiagonal', 'PermuteDims', 'ArrayAdd',
              'ArrayElementwiseApplyFunc', 'Reshape', 'ZeroArray', 'OneArray', 'ArrayElement']:
        ns[n] = getattr(ae, n)
    try:
        ns['array_derive'] = importlib.import_module(mod + '.tensor.array.expressions.arrayexpr_derivatives').array_derive
    except ImportError:  # base tree: array_derive is only exported by the package
        ns['array_derive'] = ae.array_derive
    fi = importlib.import_module(mod + '.tensor.array.expressions.from_indexed_to_array')
    ns['c'], ns['c_'] = fi.convert_indexed_to_array, getattr(fi, '_convert_indexed_to_array', None)
    ns['a2i'] = importlib.import_module(mod + '.tensor.array.expressions.from_array_to_indexed').convert_array_to_indexed
    for n in ['MatrixSymbol', 'Identity', 'OneMatrix', 'ZeroMatrix', 'Trace', 'HadamardPower', 'hadamard_product',
              'KroneckerProduct', 'KroneckerDelta', 'IndexedBase', 'Sum', 'Lambda', 'sin', 'cos', 'exp', 'symbols',
              'Symbol', 'Inverse', 'Function', 'sqrt']:
        ns[n] = getattr(m, n)
    i, j, k, l, d, x, n_ = m.symbols('i j k l d x n')
    ns.update(i=i, j=j, k=k, l=l, d=d, x=x, n=n_)
    MS, AS = m.MatrixSymbol, ae.ArraySymbol
    ns.update(X=MS('X', k, k), Y=MS('Y', k, k), a=MS('a', k, 1), b=MS('b', k, 1),
              M=MS('M', d, d), N=MS('N', d, d), P=MS('P', d, d), Q=MS('Q', d, d), v=MS('v', d, 1), w=MS('w', d, 1),
              A=AS('A', (k, k, k)), B=AS('B', (k, k)), C=AS('C', (k,)), E=AS('E', (3, 2)), F=AS('F', (2, 3)),
              IA=m.IndexedBase('IA'), IB=m.IndexedBase('IB'), IC=m.IndexedBase('IC'), IM=m.IndexedBase('IM'),
              SA=m.IndexedBase('SA', shape=(d, d)), SB=m.IndexedBase('SB', shape=(d, d)), SV=m.IndexedBase('SV', shape=(d,)))
    ns['dd'] = m.Dummy('dd')
    return ns


def ser(e, m):
    n = type(e).__name__
    if n == 'PermuteDims':
        return 'PermuteDims(%s, %s)' % (ser(e.expr, m), list(e.permutation.array_form))
    if n in ('ArrayTensorProduct', 'ArrayAdd'):
        return '%s(%s)' % (n, ', '.join(ser(a, m) for a in e.args))
    if n in ('ArrayContraction', 'ArrayDiagonal'):
        return '%s(%s, %s)' % (n, ser(e.expr, m), ', '.join(str(tuple(int(j) for j in g)) for g in e.args[1:]))
    if n == 'Reshape':
        return 'Reshape(%s, %s)' % (ser(e.expr, m), tuple(e.shape))
    if n == 'ArrayElementwiseApplyFunc':
        return 'Apply(%s, %s)' % (ser(e.expr, m), e.function(m.Symbol('t')))
    if n == 'ArrayElement':
        return 'ArrayElement(%s, %s)' % (ser(e.name, m), list(e.indices))
    if n == 'ZeroArray':
        return 'ZeroArray%s' % (tuple(e.shape),)
    if n == 'OneArray':
        return 'OneArray%s' % (tuple(e.shape),)
    if n == 'ElementwiseApplyFunction':
        return 'EWF(%s, %s)' % (e.function(m.Symbol('t')), ser(e.expr, m))
    if isinstance(e, tuple):
        return '(%s)' % ', '.join(ser(a, m) for a in e)
    if isinstance(e, frozenset):
        return 'frozenset(%s)' % sorted(map(str, e))
    return '%s:%s' % (n, e)


DERIVE = {
 'msym_self': 'array_derive(X, X)', 'msym_other': 'array_derive(Y, X)', 'transpose': 'array_derive(X.T, X)',
 'mul': 'array_derive(X*Y, X)', 'mul_self': 'array_derive(X*X, X)', 'pow2': 'array_derive(X**2, X)',
 'pow3': 'array_derive(X**3, X)', 'inv': 'array_derive(X.I, X)', 'inv_mul': 'array_derive(Inverse(X*Y), X)',
 'bilinear': 'array_derive(a.T*X*b, X)', 'trace': 'array_derive(Trace(X), X)', 'trace_mul': 'array_derive(Trace(X)*X, X)',
 'ident': 'array_derive(Identity(k), X)', 'ones': 'array_derive(OneMatrix(k, k), X)', 'zero': 'array_derive(ZeroMatrix(k, k), X)',
 'ewf': 'array_derive(X.applyfunc(sin), X)', 'ewf_mul': 'array_derive(a.T*X.applyfunc(sin)*b, X)',
 'hadamard': 'array_derive(hadamard_product(X, Y), X)', 'hadpow': 'array_derive(HadamardPower(X, 2), X)',
 'kron': 'array_derive(KroneckerProduct(X, Y), X)', 'smul': 'array_derive(x*X, X)', 'axbt': 'array_derive(X*a*b.T, X)',
 'xtx': 'array_derive(X.T*X, X)', 'xax': 'array_derive(X*Y*X, X)', 'add': 'array_derive(X + Y*X, X)',
 'arr_self3': 'array_derive(A, A)', 'arr_self2': 'array_derive(B, B)', 'arr_self1': 'array_derive(C, C)',
 'arr_other': 'array_derive(A, B)', 'arr_wrt_matrix': 'array_derive(B, X)', 'matrix_wrt_arr': 'array_derive(X, B)',
 'tp': 'array_derive(ArrayTensorProduct(A, B), B)', 'tp_self': 'array_derive(ArrayTensorProduct(B, B), B)',
 'tp_scalar': 'array_derive(ArrayTensorProduct(x, B), B)', 'tp_int': 'array_derive(ArrayTensorProduct(2, B), B)',
 'tp_three': 'array_derive(ArrayTensorProduct(C, B, C), C)',
 'contr': 'array_derive(ArrayContraction(ArrayTensorProduct(B, B), (1, 2)), B)',
 'contr_trace': 'array_derive(ArrayContraction(B, (0, 1)), B)',
 'diag': 'array_derive(ArrayDiagonal(ArrayTensorProduct(B, B), (1, 2)), B)',
 'perm': 'array_derive(PermuteDims(A, [1, 2, 0]), A)', 'add_arr': 'array_derive(ArrayAdd(B, PermuteDims(B, [1, 0])), B)',
 'reshape': 'array_derive(Reshape(B, (k**2,)), B)', 'aeaf': 'array_derive(ArrayElementwiseApplyFunc(sin, B), B)',
 'aeaf_lambda': 'array_derive(ArrayElementwiseApplyFunc(Lambda(dd, dd**2), C), C)',
 'aeaf_exp': 'array_derive(ArrayElementwiseApplyFunc(exp, C), C)',
 'aeaf_other': 'array_derive(ArrayElementwiseApplyFunc(sin, C), B)',
 'zeroarray': 'array_derive(ZeroArray(k, k), B)', 'onearray': 'array_derive(OneArray(k), B)',
 'scalar_x': 'array_derive(x, B)', 'arr_wrt_symbol': 'array_derive(B, x)', 'tp_wrt_symbol': 'array_derive(ArrayTensorProduct(x, B), x)',
 'err_int': 'array_derive(2, B)', 'err_ewf_rank': 'array_derive(X.applyfunc(sin), C)',
 'contr_mat': 'array_derive(ArrayContraction(ArrayTensorProduct(X, Y), (1, 2)), X)',
 'diag_mat': 'array_derive(ArrayDiagonal(ArrayTensorProduct(X, Y), (0, 2)), Y)',
 'contr_vec': 'array_derive(ArrayContraction(ArrayTensorProduct(C, B, C), (0, 1), (2, 3)), C)',
 'perm_mat': 'array_derive(PermuteDims(ArrayTensorProduct(X, Y), [2, 3, 0, 1]), Y)',
}

INDEXED = {
 'mm': 'c(Sum(M[i, j]*N[j, k], (j, 0, d - 1)))', 'mm_first': 'c(Sum(M[i, j]*N[j, k], (j, 0, d - 1)), [k])',
 'trace': 'c(Sum(M[i, i], (i, 0, d - 1)))', 'diag': 'c(M[i, i])', 'elem': 'c(M[i, j])', 'elem_first': 'c(M[i, j], [j, i])',
 'tp': 'c(M[i, j]*N[k, l])', 'tp_first': 'c(M[i, j]*N[k, l], [l, k, j, i])',
 'mmm': 'c(Sum(M[i, j]*N[j, k]*P[k, l], (j, 0, d - 1), (k, 0, d - 1)))',
 'trace_mm': 'c(Sum(M[i, j]*N[j, i], (i, 0, d - 1), (j, 0, d - 1)))',
 'tr_mm2': 'c(Sum(M[i, j]*N[i, j], (i, 0, d - 1), (j, 0, d - 1)))',
 'hadamard': 'c(M[i, j]*N[i, j])', 'hadamard3': 'c(M[i, j]*N[i, j]*P[i, j])',
 'add': 'c(M[i, j] + N[i, j])', 'add_t': 'c(M[i, j] + N[j, i])', 'add_first': 'c(M[i, j] + N[j, i], [j, i])',
 'add_bcast': 'c(M[i, j] + v[i, 0])', 'add_scalar': 'c(M[i, j]*N[k, l] + M[i, l]*N[k, j])',
 'add_sum': 'c(Sum(M[i, j]*N[j, k], (j, 0, d - 1)) + P[i, k])',
 'kd': 'c(M[i, j]*KroneckerDelta(j, k))', 'kd_sum': 'c(Sum(M[i, j]*KroneckerDelta(j, k), (j, 0, d - 1)))',
 'kd2': 'c(M[i, j]*KroneckerDelta(i, j))', 'kd_first': 'c(M[i, j]*KroneckerDelta(j, k), [k])',
 'kd_chain': 'c(M[i, j]*KroneckerDelta(j, k)*KroneckerDelta(k, l))',
 'kd_alone': 'c(KroneckerDelta(i, j))', 'kd_two': 'c(M[i, j]*N[k, l]*KroneckerDelta(j, k))',
 'pow': 'c(M[i, j]**2)', 'pow3': 'c(M[i, j]**3)', 'pow_neg': 'c(M[i, j]**-1)', 'pow_half': 'c(M[i, j]**(1/2))',
 'pow_sym': 'c(M[i, j]**k)', 'func': 'c(sin(M[i, j]))', 'func_sum': 'c(Sum(sin(M[i, j])*N[j, k], (j, 0, d - 1)))',
 'scalar': 'c(2*M[i, j])', 'scalar_sym': 'c(x*M[i, j])', 'scalar_only': 'c(x)', 'int_index': 'c(M[i, 0])',
 'int_index_first': 'c(M[i, 0], [i])', 'int_mul': 'c(M[i, 0]*N[0, j])', 'int_mul_first': 'c(M[i, 0]*N[0, j], [j, i])',
 'vec_mm': 'c(Sum(v[i, 0]*M[i, j]*w[j, 0], (i, 0, d - 1), (j, 0, d - 1)))',
 'err_mismatch': 'c(Sum(M[i, j]*N[j, k], (j, 0, d - 2)))', 'err_start': 'c(Sum(M[i, j]*N[j, k], (j, 1, d - 1)))',
 'err_notfound': 'c(Sum(M[i, j], (k, 0, d - 1)))', 'err_base': 'c(IA)',
 'sum_partial': 'c(Sum(M[i, j], (j, 0, d - 1)))', 'sum_diag_free': 'c(Sum(M[i, i]*N[i, j], (j, 0, d - 1)))',
 'sum_diag_keep': 'c(Sum(M[i, j]*N[j, j], (i, 0, d - 1)))',
 'idx_shaped': 'c(Sum(SA[i, j]*SB[j, k], (j, 0, d - 1)))', 'idx_shaped_diag': 'c(SA[i, i])',
 'idx_shaped_vec': 'c(Sum(SA[i, j]*SV[j], (j, 0, d - 1)))', 'idx_shaped_add': 'c(SA[i, j] + SB[j, i])',
 'u_tp': 'c_(IA[i]*IB[j])', 'u_diag': 'c_(IA[i]*IB[i])', 'u_contr': 'c_(Sum(IA[i]*IB[i], (i, 0, 3)))',
 'u_mm': 'c_(IA[i, j]*IB[j, k])', 'u_mm_sum': 'c_(Sum(IA[i, j]*IB[j, k], (j, 0, 3)))', 'u_trace_d': 'c_(IA[i, i])',
 'u_trace': 'c_(Sum(IA[i, i], (i, 0, 3)))', 'u_mixed': 'c_(IA[i]*M[i, j])', 'u_add': 'c_(IA[i] + IB[i])',
 'u_add2': 'c_(IA[i] + IB[j])', 'u_add3': 'c_(IA[i, j] + IB[j, i])', 'u_pow': 'c_(IA[i]**2)', 'u_func': 'c_(sin(IA[i]))',
 'u_tp3': 'c_(IA[i]*IB[j]*IC[k])', 'u_sum3': 'c_(Sum(IA[i]*IB[i]*IC[j], (i, 0, 3)))', 'u_cross': 'c_(IA[i, j]*IB[j, i])',
 'u_cross_sum': 'c_(Sum(IA[i, j]*IB[j, i], (i, 0, 3), (j, 0, 3)))', 'u_scal': 'c_(2*IA[i])', 'u_diag3': 'c_(IA[i]*IB[j]*IC[i])',
 'u_sum1': 'c_(Sum(IA[i], (i, 0, 3)))', 'u_int': 'c_(IA[0])', 'u_int_mul': 'c_(IA[i, 0]*IB[i])', 'u_kd': 'c_(KroneckerDelta(i, j)*IA[i])',
 'u_vmv': 'c_(IA[i]*IB[j]*IM[i, j])', 'u_perm_err': 'c(IA[i]*IB[j], [j, i])', 'u_perm_err2': 'c(IA[i, j], [j, i])',
 'ae_elem': 'c(ArrayElement(B, (i, j)))', 'ae_mm': 'c(Sum(ArrayElement(B, (i, j))*ArrayElement(B, (j, l)), (j, 0, k - 1)))',
 'ae_diag': 'c(ArrayElement(A, (i, i, j)))', 'ae_first': 'c(ArrayElement(A, (i, j, l)), [l, i])',
 'raw_mm': 'c_(Sum(M[i, j]*N[j, k], (j, 0, d - 1)))', 'raw_kd': 'c_(M[i, j]*KroneckerDelta(j, k))',
 'raw_add': 'c_(M[i, j] + N[j, i])', 'raw_int': 'c_(M[i, 0])',
}

TO_INDEXED = {
 'sym': 'a2i(B, [i, j])', 'sym3': 'a2i(A, [i, j, l])', 'tp': 'a2i(ArrayTensorProduct(B, C), [i, j, l])',
 'contr': 'a2i(ArrayContraction(ArrayTensorProduct(X, Y), (1, 2)), [i, j])',
 'contr2': 'a2i(ArrayContraction(ArrayTensorProduct(X, Y), (0, 3), (1, 2)), [])',
 'diag': 'a2i(ArrayDiagonal(ArrayTensorProduct(X, Y), (1, 2)), [i, j, l])',
 'diag2': 'a2i(ArrayDiagonal(A, (0, 2)), [i, j])', 'perm': 'a2i(PermuteDims(A, [1, 2, 0]), [i, j, l])',
 'add': 'a2i(ArrayAdd(B, PermuteDims(B, [1, 0])), [i, j])', 'aeaf': 'a2i(ArrayElementwiseApplyFunc(sin, B), [i, j])',
 'reshape': 'a2i(Reshape(B, (k**2,)), [i])', 'reshape2': 'a2i(Reshape(E, (2, 3)), [i, j])',
 'reshape3': 'a2i(Reshape(A, (k, k**2)), [i, j])', 'reshape_flat': 'a2i(Reshape(E, (6,)), [i])',
 'zero': 'a2i(ZeroArray(k, k), [i, j])', 'one': 'a2i(OneArray(k), [i])', 'matrix': 'a2i(X, [i, j])',
 'numbers': 'a2i(ArrayTensorProduct(B, C), [0, 1, 2])', 'scalar_tp': 'a2i(ArrayTensorProduct(x, B), [i, j])',
 'nested': 'a2i(ArrayContraction(ArrayDiagonal(ArrayTensorProduct(A, B), (0, 3)), (0, 1)), [i, j])',
 'err_bounds': 'a2i(B, [k, 0])', 'err_neg': 'a2i(B, [-1, 0])', 'err_count': 'a2i(B, [i])',
 'contr_mat_vec': 'a2i(ArrayContraction(ArrayTensorProduct(X, a), (1, 2)), [i, j])',
 'diag_trace': 'a2i(ArrayContraction(ArrayDiagonal(ArrayTensorProduct(X, Y, X), (1, 2)), (0, 3)), [i, j])',
}

SECTIONS = {'derive': DERIVE, 'indexed': INDEXED, 'to_indexed': TO_INDEXED}


def evaluate(mod, section, name):
    ns = namespace(mod)
    try:
        return ser(eval(SECTIONS[section][name], ns), ns['m'])
    except Exception as e:
        return 'ERR %s: %s' % (type(e).__name__, str(e)[:80])


if __name__ == '__main__':
    import json
    import sys
    mod = sys.argv[1]
    out = {s: {n: evaluate(mod, s, n) for n in sorted(cases)} for s, cases in SECTIONS.items()}
    json.dump(out, sys.stdout, indent=1, sort_keys=True)
