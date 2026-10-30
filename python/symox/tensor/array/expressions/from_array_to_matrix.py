"""``sympy.tensor.array.expressions.from_array_to_matrix``.

SymPy 1.14's ``convert_array_to_matrix`` and its passes, computed in Rust
(``crates/symox-core/src/arrayexpr_a2m.rs``). The functions here only admit
the input, encode it (``_a2m``) and rebuild the result. An input outside the
native port (an object kind it does not model) runs the stock SymPy
algorithm (``symox._arrayexpr``) on a structural translation (``_stock``).
"""
from . import ArrayElement
from . import _a2m
from ._stock import stock_call
from .... import _core


def _native(op, exprs, extra=None, removed=False, stock_name=None, stock_args=None):
    try:
        return _a2m.run(op, exprs, extra, removed)
    except (_a2m._Decline, _a2m.Unsupported):
        return stock_call('from_array_to_matrix', stock_name, *(exprs if stock_args is None else stock_args))


def _array2matrix(expr):
    if isinstance(expr, ArrayElement):
        ret = _array2matrix(expr.name)
        if isinstance(ret, _core._MatrixExpr):
            from ....matrices.expressions.matexpr import MatrixElement
            return MatrixElement(ret, *expr.indices)
        return ArrayElement(ret, expr.indices)
    return _native('a2m', [expr], stock_name='_array2matrix')


def _remove_trivial_dims(expr):
    if isinstance(expr, ArrayElement):
        return expr, []
    return _native('rtd', [expr], removed=True, stock_name='_remove_trivial_dims')


def convert_array_to_matrix(expr):
    r"""
    Recognize matrix expressions in codegen objects.

    If more than one matrix multiplication line have been detected, return a
    list with the matrix expressions.

    Examples
    ========

    >>> from symox.tensor.array.expressions.from_indexed_to_array import convert_indexed_to_array
    >>> from symox.tensor.array import tensorcontraction, tensorproduct
    >>> from symox import MatrixSymbol, Sum
    >>> from symox.abc import i, j, k, l, N
    >>> from symox.tensor.array.expressions.from_matrix_to_array import convert_matrix_to_array
    >>> from symox.tensor.array.expressions.from_array_to_matrix import convert_array_to_matrix
    >>> M = MatrixSymbol("M", N, N)
    >>> P = MatrixSymbol("P", N, N)
    >>> expr = M*P
    >>> cg = convert_matrix_to_array(expr)
    >>> convert_array_to_matrix(cg)
    M*P
    >>> expr = M.T*P
    >>> cg = convert_matrix_to_array(expr)
    >>> convert_array_to_matrix(cg)
    M.T*P
    """
    if isinstance(expr, ArrayElement):
        return _remove_trivial_dims(_array2matrix(expr))[0]
    return _native('convert', [expr], stock_name='convert_array_to_matrix')


def _support_function_tp1_recognize(contraction_indices, args):
    ci = [[int(j) for j in g] for g in contraction_indices]
    return _native('tp1', list(args), extra=ci, stock_name='_support_function_tp1_recognize',
                   stock_args=[contraction_indices, args])


def identify_hadamard_products(expr):
    return _native('hadamard', [expr], stock_name='identify_hadamard_products')


def identify_removable_identity_matrices(expr):
    return _native('removable', [expr], stock_name='identify_removable_identity_matrices')


def remove_identity_matrices(expr):
    return _native('remove_identity', [expr], removed=True, stock_name='remove_identity_matrices')


def _array_contraction_to_diagonal_multiple_identity(expr):
    return _native('multi_identity', [expr], removed=True, stock_name='_array_contraction_to_diagonal_multiple_identity')


def _array_diag2contr_diagmatrix(expr):
    return _native('diag2contr', [expr], stock_name='_array_diag2contr_diagmatrix')


def _combine_removed(dim, removed1, removed2):
    return _a2m.combine_removed(removed1, removed2)


def _split_multiple_contractions(expr):
    try:
        return _a2m.run('split', [expr])
    except (_a2m._Decline, _a2m.Unsupported):
        from ._stock import to_stock, from_stock
        return from_stock(to_stock(expr).split_multiple_contractions())


def _flatten_contraction_of_diagonal(expr):
    try:
        return _a2m.run('flatten', [expr])
    except (_a2m._Decline, _a2m.Unsupported):
        from ._stock import to_stock, from_stock
        return from_stock(to_stock(expr).flatten_contraction_of_diagonal())


def __getattr__(name):
    # The other helpers of SymPy's module run as the stock implementation.
    import importlib
    mod = importlib.import_module('symox._arrayexpr.from_array_to_matrix')
    if name.startswith('__') or not callable(getattr(mod, name, None)):
        raise AttributeError(name)
    return lambda *a: stock_call('from_array_to_matrix', name, *a)
