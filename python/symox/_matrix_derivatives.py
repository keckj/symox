"""Derivatives of matrix expressions by a MatrixSymbol.

A port of SymPy's matrix-line algorithm (``_LeftRightArgs`` and the
``_eval_derivative_matrix_lines`` methods of the matrix-expression classes,
used by ``_matrix_derivative``). When an expression is derived by a matrix
``X[m, n]``, two lines of matrix products are built: the one contracted to
``m`` (first line) and the one contracted to ``n`` (second line).
Transposition swaps the lines, a trace joins their ends.

Results of rank at most two are converted back to matrix expressions, as
in SymPy. A derivative of a matrix by a matrix (rank four) is returned as
SymPy's array expression ``PermuteDims(ArrayTensorProduct(first, second))``
when each line is a plain matrix. Every other case raises ``_Unsupported``
and the caller keeps its previous (explicit-entry) derivative.
"""
from . import _core
from ._kind import matrix_tag


class _Unsupported(Exception):
    pass


def _S():
    import symox
    return symox


class _ExprBuilder:
    """SymPy's ``ExprBuilder``: a deferred call whose arguments are mutable."""

    def __init__(self, op, args, validator=None):
        self.op = op
        self.args = list(args)
        self.validator = validator

    def build(self):
        args = [a.build() if isinstance(a, _ExprBuilder) else a for a in self.args]
        if self.validator is not None:
            self.validator(*args)
        return self.op(*args)


def _build(expr):
    if isinstance(expr, _ExprBuilder):
        return expr.build()
    return expr


def _is_matrix(value):
    return isinstance(value, _core._MatrixExpr) or bool(getattr(value, 'is_Matrix', False))


def _times(a, b):
    """``a *= b`` on a line slot, with SymPy's evaluated products."""
    S = _S()
    if a == S.Integer(1):
        return b
    if b == S.Integer(1):
        return a
    if _is_matrix(a) or _is_matrix(b):
        # MatrixExpr.__mul__: MatMul(a, b).doit()
        return S.MatMul(a, b).doit()
    return a * b


class _LeftRightArgs:
    def __init__(self, lines, higher=None):
        S = _S()
        self._lines = list(lines)
        self._first_pointer_parent = self._lines
        self._first_pointer_index = 0
        self._first_line_index = 0
        self._second_pointer_parent = self._lines
        self._second_pointer_index = 1
        self._second_line_index = 1
        self.higher = S.Integer(1) if higher is None else higher

    @property
    def first_pointer(self):
        return self._first_pointer_parent[self._first_pointer_index]

    @first_pointer.setter
    def first_pointer(self, value):
        self._first_pointer_parent[self._first_pointer_index] = value

    @property
    def second_pointer(self):
        return self._second_pointer_parent[self._second_pointer_index]

    @second_pointer.setter
    def second_pointer(self, value):
        self._second_pointer_parent[self._second_pointer_index] = value

    def transpose(self):
        self._first_pointer_parent, self._second_pointer_parent = self._second_pointer_parent, self._first_pointer_parent
        self._first_pointer_index, self._second_pointer_index = self._second_pointer_index, self._first_pointer_index
        self._first_line_index, self._second_line_index = self._second_line_index, self._first_line_index
        return self

    def build(self):
        data = [_build(i) for i in self._lines]
        if not _is_one(self.higher):
            data.append(_build(self.higher))
        return data

    def append_first(self, other):
        self.first_pointer = _times(self.first_pointer, other)

    def append_second(self, other):
        self.second_pointer = _times(self.second_pointer, other)


def _is_one(value):
    return not isinstance(value, _ExprBuilder) and not _is_matrix(value) and value == 1


def _identity_or_one(n):
    S = _S()
    return S.Identity(n) if n != 1 else S.Integer(1)


def _matrix_symbol_lines(expr, x):
    S = _S()
    if expr != x:
        first = S.ZeroMatrix(x.shape[0], expr.shape[0]) if expr.shape[0] != 1 else S.Integer(0)
        second = S.ZeroMatrix(x.shape[1], expr.shape[1]) if expr.shape[1] != 1 else S.Integer(0)
        return [_LeftRightArgs([first, second])]
    return [_LeftRightArgs([_identity_or_one(expr.shape[0]), _identity_or_one(expr.shape[1])])]


def _matmul_lines(expr, x):
    S = _S()
    args = expr.args
    lines = []
    for ind, arg in enumerate(args):
        if not arg.has(x):
            continue
        left_args = args[:ind]
        right_args = args[ind + 1:]
        right_mat = _product(right_args) if right_args else S.Identity(expr.shape[1])
        if left_args:
            left_rev = _product([S.Transpose(i).doit() if _is_matrix(i) else i for i in reversed(left_args)])
        else:
            left_rev = S.Identity(expr.shape[0])
        for line in _lines(arg, x):
            line.append_first(left_rev)
            line.append_second(right_mat)
            lines.append(line)
    return lines


def _product(factors):
    S = _S()
    factors = list(factors)
    if len(factors) == 1:
        return factors[0]
    return S.MatMul(*factors).doit(deep=False) if any(_is_matrix(f) for f in factors) else S.Mul(*factors)


def _trace_lines(expr, x):
    from .tensor.array.expressions import ArrayTensorProduct, _array_contraction
    S = _S()
    lines = _lines(expr.args[0], x)
    for lr in lines:
        if _is_one(lr.higher):
            lr.higher = _ExprBuilder(
                _array_contraction,
                [_ExprBuilder(ArrayTensorProduct, [lr._lines[0], lr._lines[1]]), (1, 3)])
        else:
            lr.higher = _ExprBuilder(
                _array_contraction,
                [_ExprBuilder(ArrayTensorProduct, [lr._lines[0], lr._lines[1], lr.higher]), (1, 3), (0, 2)])
        lr._lines = [S.Integer(1), S.Integer(1)]
        lr._first_pointer_parent = lr._lines
        lr._second_pointer_parent = lr._lines
        lr._first_pointer_index = 0
        lr._second_pointer_index = 1
    return lines


def _inverse_lines(expr, x):
    S = _S()
    lines = _lines(expr.args[0], x)
    for line in lines:
        line.first_pointer = _times(line.first_pointer, -expr.T)
        line.second_pointer = _times(line.second_pointer, expr)
    return lines


def _matpow_lines(expr, x):
    S = _S()
    base, exp = expr.args
    if base.shape == (1, 1) and not exp.has(x):
        return _power_1x1_lines(base, exp, x)
    if (exp > 0) == True:
        newexpr = S.MatMul(*[base for _ in range(int(exp))])
    elif (exp == -1) == True:
        return _inverse_lines(S.Inverse(base), x)
    elif (exp < 0) == True:
        newexpr = S.MatMul(*[S.Inverse(base) for _ in range(int(-exp))])
    elif (exp == 0) == True:
        return _lines(expr.doit(), x)
    else:
        raise NotImplementedError("cannot evaluate %s derived by %s" % (expr, x))
    return _matmul_lines(newexpr, x)


def _power_1x1_chain(first, line0, coefficient, line1, second):
    S = _S()
    if not all(_is_one(p) or (_is_matrix(p) and p == S.Identity(1)) for p in (first, second)):
        raise _Unsupported
    if _is_matrix(line1):
        line1 = line1.T
    result = _times(line0, coefficient)
    if not (_is_one(line1) or (_is_matrix(line1) and line1 == S.Identity(1))):
        result = _times(result, line1)
    return result


def _power_1x1_lines(base, exp, x):
    """MatPow._eval_derivative_matrix_lines for a 1x1 base: one contracted chain."""
    S = _S()
    lines = _lines(base, x)
    for i in lines:
        if not _is_one(i.higher) or len(i._lines) != 2:
            raise _Unsupported
        subexpr = _ExprBuilder(_power_1x1_chain, [
            S.Identity(1), i._lines[0], exp * base**(exp - 1), i._lines[1], S.Identity(1)])
        i._first_pointer_parent = subexpr.args
        i._first_pointer_index = 0
        i._second_pointer_parent = subexpr.args
        i._second_pointer_index = 4
        i._lines = [subexpr]
    return lines


def _lines(expr, x):
    S = _S()
    tag = matrix_tag(expr)
    if tag == 'MatrixSymbol':
        return _matrix_symbol_lines(expr, x)
    if tag == 'MatAdd':
        return [line for arg in expr.args for line in _lines(arg, x)]
    if tag == 'MatMul':
        return _matmul_lines(expr, x)
    if tag == 'Transpose':
        return [line.transpose() for line in _lines(expr.args[0], x)]
    if tag == 'Inverse':
        return _inverse_lines(expr, x)
    if tag == 'MatPow':
        return _matpow_lines(expr, x)
    if tag is not None:
        if not expr.has(x):
            return [_LeftRightArgs([S.Integer(1), S.Integer(1)], higher=S.Integer(0))]
        raise _Unsupported
    if getattr(expr.func, '__name__', None) == 'Trace':
        return _trace_lines(expr, x)
    if (getattr(expr.func, '__name__', None) == 'Pow' and matrix_tag(expr.args[0]) is not None
            and tuple(expr.args[0].shape) == (1, 1) and not expr.args[1].has(x)):
        # SymPy keeps a power of a 1x1 matrix expression as a MatPow.
        return _power_1x1_lines(expr.args[0], expr.args[1], x)
    # Expr._eval_derivative_matrix_lines
    return [_LeftRightArgs([S.Integer(1), S.Integer(1)], higher=_scalar_derivative(expr, x))]


def _derive(expr, x):
    """``expr.diff(x)`` for a scalar factor of a product rule."""
    S = _S()
    if not getattr(expr, 'has', None) or not expr.has(x):
        return S.Integer(0)
    return matrix_derivative(expr, x)


def _scalar_derivative(expr, x):
    """SymPy's ``Expr._eval_derivative`` for the scalar heads that hold matrices."""
    S = _S()
    if not expr.has(x):
        return S.Integer(0)
    name = getattr(expr.func, '__name__', None)
    if name == 'Add':
        terms = [d for d in (_derive(a, x) for a in expr.args) if not _is_zero(d)]
        return _sum(terms)
    if name == 'Mul':
        args = list(expr.args)
        terms = []
        for i in range(len(args)):
            d = _derive(args[i], x)
            if not _is_zero(d):
                factors = args[:i] + [d] + args[i + 1:]
                term = S.Integer(1)
                for f in factors:
                    term = _times(term, f)
                terms.append(term)
        return _sum(terms)
    if name == 'Pow':
        base, exp = expr.args
        if exp.has(x):
            raise _Unsupported
        dbase = _derive(base, x)
        if _is_zero(dbase):
            return S.Integer(0)
        return expr * (dbase * exp / base)
    raise _Unsupported


def _is_zero(value):
    if _is_matrix(value):
        return matrix_tag(value) == 'ZeroMatrix'
    return value == 0


def _sum(terms):
    S = _S()
    if not terms:
        return S.Integer(0)
    result = terms[0]
    for t in terms[1:]:
        result = result + t
    return result


def _shape_of(elem):
    if _is_matrix(elem):
        return tuple(elem.shape)
    return (1, 1)


def _contract_one_dims(parts):
    S = _S()
    if len(parts) == 1:
        return parts[0]
    p1, p2 = parts[:2]
    if len(parts) == 2 and _is_matrix(p1) and _is_matrix(p2) and p1.shape[1] != p2.shape[1]:
        # A vector derived by a vector: SymPy's array conversion gives
        # these products (the lines are not both contracted to the
        # columns).
        if tuple(p2.shape) == (1, 1):
            return p1 if p2 == S.Identity(1) else S.KroneckerProduct(p1, p2)
        if p1.shape[1] == p2.shape[0]:
            return _times(p1, p2)
        raise _Unsupported
    if _is_matrix(p2):
        p2 = p2.T
    if _is_matrix(p1) and p1 == S.Identity(1):
        pbase = p2
    elif _is_matrix(p2) and p2 == S.Identity(1):
        pbase = p1
    else:
        pbase = _times(p1, p2)
    if len(parts) == 2:
        return pbase
    if _is_matrix(pbase):
        raise _Unsupported
    rest = S.Integer(1)
    for p in parts[2:]:
        rest = _times(rest, p)
    return _times(pbase, rest)


def _to_matrix(value):
    from .tensor.array.expressions import _chain_array_to_matrix, _ArrayExpr
    result = _chain_array_to_matrix(value) if isinstance(value, _ArrayExpr) else value
    if isinstance(result, _ArrayExpr):
        raise _Unsupported
    return result


def _rank_four(lines, x):
    """SymPy's array form of a matrix derived by a matrix, for plain lines."""
    from .tensor.array.expressions import ArrayAdd, ArrayTensorProduct, _permute_dims
    terms = []
    for line in lines:
        if len(line._lines) != 2:
            raise _Unsupported
        if not _is_one(line.higher):
            # A scalar factor derived by x (Trace(x)*M): the derivative of
            # the scalar (axes of x) times the matrix the two lines contract
            # to (axes of the expression).
            # (Trace keeps the line indices: the pointers give the sides.)
            if not (line._first_pointer_parent is line._lines and line._second_pointer_parent is line._lines):
                raise _Unsupported
            higher = _to_matrix(_build(line.higher))
            first = line.first_pointer
            second = line.second_pointer
            if not (_is_matrix(higher) and tuple(higher.shape) == tuple(x.shape)
                    and _is_matrix(first) and _is_matrix(second)):
                raise _Unsupported
            if _is_zero(higher) or _is_zero(first) or _is_zero(second):
                continue
            first = first.T
            matrix = second if first == _S().Identity(first.shape[0]) else _times(first, second)
            terms.append(ArrayTensorProduct(higher, matrix))
            continue
        first = line._lines[line._first_line_index]
        second = line._lines[line._second_line_index]
        if not (_is_matrix(first) and _is_matrix(second)):
            raise _Unsupported
        if matrix_tag(first) == 'ZeroMatrix' or matrix_tag(second) == 'ZeroMatrix':
            continue
        # lines[0] carries the x row axis, lines[1] the x column axis;
        # _first_line_index names the line that carries the expression's
        # row axis.  The result orders the axes (x row, x column, expr
        # row, expr column).
        order = [0, 2, 1, 3] if line._first_line_index == 0 else [0, 2, 3, 1]
        terms.append(_permute_dims(ArrayTensorProduct(line._lines[0], line._lines[1]), order))
    if not terms:
        raise _Unsupported
    return terms[0] if len(terms) == 1 else ArrayAdd(*terms)


_ELEMENTWISE = frozenset(('HadamardProduct', 'HadamardPower', 'ElementwiseApplyFunction'))


def _has_elementwise(expr):
    """Whether an element-wise matrix operation occurs in ``expr``."""
    stack = [expr]
    while stack:
        node = stack.pop()
        if matrix_tag(node) in _ELEMENTWISE:
            return True
        stack.extend(a for a in getattr(node, 'args', ()) if hasattr(a, 'args'))
    return False


def _array_matrix_derivative(expr, x):
    """SymPy's ``_matrix_derivative`` (array-expression algorithm).

    ``convert_matrix_to_array``, ``array_derive`` and
    ``convert_array_to_matrix`` of SymPy 1.14 (``symox._arrayexpr``).
    """
    from ._arrayexpr.from_matrix_to_array import convert_matrix_to_array
    from ._arrayexpr.arrayexpr_derivatives import array_derive
    from ._arrayexpr.from_array_to_matrix import convert_array_to_matrix
    try:
        return convert_array_to_matrix(array_derive(convert_matrix_to_array(expr), x))
    except IndexError:
        # SymPy's ArrayDiagonal flattening fails on a nested
        # ElementwiseApplyFunction. SymPy's MatrixExpr.__mul__ composes them
        # (MatMul(...).doit()), so a.T*X.applyfunc(sin).applyfunc(log)*b is
        # derived from the composed form there.
        evaluated = expr.doit()
        if evaluated == expr:
            raise
        return convert_array_to_matrix(array_derive(convert_matrix_to_array(evaluated), x))


_ARRAY_PATH_ERRORS = (RecursionError, NotImplementedError, IndexError, ValueError, TypeError, AttributeError,
                      _core.ArrayExprUnsupported)


def _native_array_derivative(expr, x):
    from .tensor.array.expressions._canon import convert_matrix_to_array
    from .tensor.array.expressions._native import array_derive
    from .tensor.array.expressions.from_array_to_matrix import convert_array_to_matrix
    array = convert_matrix_to_array(expr)
    if _opaque_scalar_has(array, x):
        # Admission only: SymPy 1.14's array_derive registers Expr -> ZeroArray,
        # so a scalar leaf depending on x (exp(Trace(X)), X[0, 1],
        # DotProduct(a, a)) would yield a wrong zero. Those inputs keep the
        # previous derivation chain, which differentiates them.
        raise NotImplementedError
    return convert_array_to_matrix(array_derive(array, x))


def _opaque_scalar_has(node, x):
    from .tensor.array.expressions import _ArrayExpr
    if isinstance(node, _ArrayExpr):
        return any(_opaque_scalar_has(a, x) for a in node.args)
    if isinstance(node, _core._MatrixExpr):
        return False
    return isinstance(node, _core.Expr) and node.has(x)


def matrix_derivative(expr, x):
    """``expr.diff(x)`` for a matrix expression or scalar and a MatrixSymbol ``x``."""
    S = _S()
    if not _is_matrix(expr) and not expr.has(x):
        return S.ZeroMatrix(*x.shape)
    try:
        # SymPy 1.14's _matrix_derivative: convert_matrix_to_array,
        # array_derive and convert_array_to_matrix, each computed in Rust
        # (arrayexpr_matrix.rs, arrayexpr_derive.rs, arrayexpr_a2m.rs).
        return _native_array_derivative(expr, x)
    except _ARRAY_PATH_ERRORS:
        # Where SymPy's array algorithm fails (e.g. RecursionError for
        # (X**-2).diff(X)), the matrix-line algorithm below still answers.
        pass
    native = _core._matrix_lines_derivative(expr, x)
    if native is not None:
        # The native matrix-line kernel (crates/symox-core matrixexpr/lines.rs)
        # covers the result; it declines on every other case.
        return native
    if _has_elementwise(expr):
        # Element-wise operations need the ArrayDiagonal forms of SymPy's
        # array-expression path (DiagMatrix and HadamardProduct results).
        return _array_matrix_derivative(expr, x)
    native = _core._matrix_lines_rank_four(expr, x)
    if native is not None:
        # Rank-four result computed natively (lines.rs, SymPy's _rank_four
        # array form); only the array-expression objects are built here.
        from .tensor.array.expressions import ArrayAdd, ArrayTensorProduct, PermuteDims
        terms = [ArrayTensorProduct(a, b) if permutation is None
                 else PermuteDims(ArrayTensorProduct(a, b), permutation)
                 for a, b, permutation in native]
        return terms[0] if len(terms) == 1 else ArrayAdd(*terms)
    lines = _lines(expr, x)
    parts = [[_to_matrix(j) for j in line.build()] for line in lines]

    def rank(part):
        return sum(j not in (1, None) for i in part for j in _shape_of(i))

    if parts and rank(parts[0]) <= 2:
        return _sum([_contract_one_dims(p) for p in parts])
    return _rank_four(lines, x)
