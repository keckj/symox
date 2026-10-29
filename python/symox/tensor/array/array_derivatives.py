"""SymPy's held ``ArrayDerivative`` of a matrix expression by a matrix.

``ArrayDerivative._dispatch_eval_derivative_n_times`` derives a matrix
expression by a matrix symbol only when the expression contains it
(``_call_derive_default``); otherwise the derivative stays unevaluated.
Only that held form is represented here: evaluated matrix derivatives go
through the matrix-line algorithm (``symox._matrix_derivatives``).
"""


def _symox():
    import symox
    return symox


class ArrayDerivative:
    """``Derivative(expr, (v, count))`` of a matrix expression by a matrix."""
    is_commutative = False
    is_scalar = False
    is_Derivative = True
    _symox_held_matrix_derivative = True
    __slots__ = ('_expr', '_variable_count')

    def __init__(self, expr, *variables):
        S = _symox()
        counts = []
        for v in variables:
            if isinstance(v, tuple) or (hasattr(v, 'is_Tuple') and getattr(v, 'is_Tuple', False)):
                var, count = v
                count = S.sympify(count)
            else:
                var, count = v, S.Integer(1)
            if counts and counts[-1][0] == var:
                counts[-1] = (var, counts[-1][1] + count)
            else:
                counts.append((var, count))
        self._expr = expr
        self._variable_count = tuple(S.Tuple(v, c) for v, c in counts)

    @property
    def expr(self):
        return self._expr

    @property
    def variable_count(self):
        return self._variable_count

    @property
    def variables(self):
        return tuple(v for v, c in self._variable_count for _ in range(int(c)))

    @property
    def args(self):
        return (self._expr,) + self._variable_count

    @property
    def func(self):
        return type(self)

    @property
    def shape(self):
        # Scalar variables add no axes (symbolic orders only occur there).
        shape = ()
        for v, c in self._variable_count:
            if hasattr(v, 'shape'):
                shape += tuple(v.shape) * int(c)
        return shape + tuple(self._expr.shape)

    @property
    def free_symbols(self):
        # Derivative.free_symbols: the expression's symbols and those of
        # the counts (the variables are bound).
        out = set(self._expr.free_symbols)
        for _, c in self._variable_count:
            out |= c.free_symbols
        return out

    def has(self, *patterns):
        return any(arg.has(*patterns) if hasattr(arg, 'has') else arg in patterns for arg in self.args)

    def doit(self, **hints):
        expr = self._expr.doit(**hints) if hints.get('deep', True) else self._expr
        if _symbolic_orders(self._variable_count):
            # Derivative.doit re-evaluates; a symbolic order stays held.
            return _evaluate(expr, self._variable_count, type(self))
        result = expr
        for v, c in self._variable_count:
            result = result.diff((v, c))
        return result

    def diff(self, *variables, **options):
        # Derivative.__new__ of a held derivative: the zero exit uses this
        # object's free symbols (S.Zero), else the orders are merged.
        pairs = _parse(variables)
        if pairs is None:
            raise NotImplementedError('derivatives of a held array derivative')
        free = self.free_symbols
        if any(c.is_positive and v not in free for v, c in pairs):
            return _symox().S.Zero
        return _evaluate(self._expr, list(self._variable_count) + pairs, ArrayDerivative)

    def _substitute(self, method, *args, **kwargs):
        bound = set(v for v, _ in self._variable_count)
        expr = self._expr
        if hasattr(expr, method):
            changed = getattr(expr, method)(*args, **kwargs)
        else:
            changed = expr
        counts = [(v, getattr(c, method)(*args, **kwargs)) for v, c in self._variable_count]
        if bound & set(getattr(changed, 'free_symbols', ())) != bound & set(getattr(expr, 'free_symbols', ())):
            return self
        # SymPy rebuilds Derivative(*newargs) unevaluated: an ArrayDerivative.
        return ArrayDerivative(changed, *counts)

    def subs(self, *args, **kwargs):
        return self._substitute('subs', *args, **kwargs)

    def xreplace(self, rule):
        return self._substitute('xreplace', rule)

    def __eq__(self, other):
        return type(other) is type(self) and self.args == other.args

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash((type(self).__name__,) + self.args)

    def __str__(self):
        parts = [str(v) if c == 1 else '(%s, %s)' % (v, c) for v, c in self._variable_count]
        return 'Derivative(%s)' % ', '.join([str(self._expr)] + parts)

    __repr__ = __str__

    def _sympystr(self, printer, *args):
        return str(self)

    def _latex(self, printer, *args, **kwargs):
        # LatexPrinter._print_Derivative with requires_partial.
        partial = sum(not getattr(s, 'is_integer', False) for s in self._expr.free_symbols) > 1
        diff_symbol = r'\partial' if partial else 'd'
        total = 0
        tex = ''
        for v, c in reversed(self._variable_count):
            total += c
            if c == 1:
                tex += r'%s %s' % (diff_symbol, printer._print(v))
            else:
                tex += r'%s %s^{%s}' % (diff_symbol, printer._print(v), printer._print(c))
        top = diff_symbol if total == 1 else r'%s^{%s}' % (diff_symbol, printer._print(total))
        body = printer._print(self._expr)
        if getattr(self._expr, '_kind_tag', None) == 'MatAdd':
            body = r'\left(%s\right)' % body
        return r'\frac{%s}{%s} %s' % (top, tex, body)

    def srepr(self):
        S = _symox()
        return 'ArrayDerivative(%s)' % ', '.join(S.srepr(arg) for arg in self.args)


class Derivative(ArrayDerivative):
    """SymPy's held ``Derivative(X, v)`` of a matrix expression.

    Without the variables in ``X``, ``Derivative.doit`` reaches
    ``MatrixExpr._eval_derivative``, which gives ``ZeroMatrix(*X.shape)``.
    With them, SymPy's ``doit`` is wrong or raises (``Derivative(X, X)``
    gives 0, ``Derivative(X*A, X)`` a TypeError); here it is the derivative
    ``X.diff(v)``, and the held form stays when that is unsupported.
    """
    __slots__ = ()

    @property
    def shape(self):
        from ... import _core
        if not isinstance(self._expr, _core._MatrixExpr):
            # SymPy's Derivative of an explicit matrix or array has no shape.
            raise AttributeError("'Derivative' object has no attribute 'shape'")
        if _symbolic_orders(self._variable_count):
            raise AttributeError("'Derivative' object has no attribute 'shape'")
        # Scalar variables add no axes.
        shape = ()
        for v, c in self._variable_count:
            shape += tuple(getattr(v, 'shape', ())) * int(c)
        return shape + tuple(self._expr.shape)

    def doit(self, **hints):
        S = _symox()
        from ... import _core
        if _symbolic_orders(self._variable_count) or not isinstance(self._expr, _core._MatrixExpr):
            expr = self._expr.doit(**hints) if hints.get('deep', True) and hasattr(self._expr, 'doit') else self._expr
            return _evaluate(expr, self._variable_count, Derivative)
        if not any(self._expr.has(v) for v, _ in self._variable_count):
            return S.ZeroMatrix(*self._expr.shape)
        expr = self._expr.doit(**hints) if hints.get('deep', True) else self._expr
        try:
            for v in self.variables:
                expr = expr.diff(v)
        except (NotImplementedError, TypeError, AttributeError):
            return self
        return expr

    def srepr(self):
        S = _symox()
        return 'Derivative(%s)' % ', '.join(S.srepr(arg) for arg in self.args)


def held_matrix_derivative(expr, variables, cls=ArrayDerivative, contained=False):
    """The held SymPy result of ``expr.diff(*variables)``, or None.

    ``expr`` is a matrix expression and every variable a MatrixSymbol (or a
    ``(MatrixSymbol, count)`` pair) that ``expr`` does not contain.  With
    ``contained`` (the ``Derivative`` constructor), the variables may occur
    in ``expr`` and may also be scalar Symbols.
    """
    from ... import _core
    pairs = []
    for v in variables:
        count = 1
        if isinstance(v, tuple):
            if len(v) != 2:
                return None
            v, count = v
        if contained and getattr(v, 'is_Symbol', False) and not isinstance(v, _core._MatrixExpr):
            pass
        elif not (isinstance(v, _core._MatrixExpr) and v._kind_tag == 'MatrixSymbol'):
            return None
        if not (isinstance(count, int) or getattr(count, 'is_Integer', False)) or int(count) <= 0:
            return None
        pairs.append((v, count))
    if not pairs or not isinstance(expr, _core._MatrixExpr):
        return None
    if not contained and any(expr.has(v) for v, _ in pairs):
        return None
    return cls(expr, *pairs)


def _symbolic_orders(pairs):
    return any(not getattr(c, 'is_Integer', False) for _, c in pairs)


def _is_explicit(expr):
    from ... import _core
    from ._array import NDimArray
    return (isinstance(expr, (NDimArray, _core.Matrix))
            or bool(getattr(expr, '_symox_sparse', False)))


def _is_mutable_matrix(expr):
    from ... import _core
    return ((isinstance(expr, _core.Matrix) or getattr(expr, '_symox_sparse', False))
            and not getattr(expr, '_immutable', False))


def _parse(variables):
    """Merged ``(Symbol, count)`` pairs of scalar Symbol variables, or None."""
    S = _symox()
    from ... import _core
    pairs = []
    for v in variables:
        count = S.Integer(1)
        if isinstance(v, int) or getattr(v, 'is_Integer', False):
            # Derivative(expr, x, 2): an order following its variable.
            if not pairs or pairs[-1][1] != 1 or int(v) < 0:
                return None
            pairs[-1] = (pairs[-1][0], S.Integer(int(v)))
            if int(v) == 0:
                pairs.pop()
            continue
        if isinstance(v, tuple) or getattr(v, 'is_Tuple', False):
            if len(v) != 2:
                return None
            v, count = v
            try:
                count = S.sympify(count)
            except Exception:
                return None
        if (isinstance(v, _core._MatrixExpr) or not getattr(v, 'is_Symbol', False)
                or getattr(count, 'is_Number', False) and not getattr(count, 'is_Integer', False)
                or not hasattr(count, 'free_symbols')):
            return None
        if count.is_negative:
            raise ValueError('order of differentiation must be nonnegative')
        if pairs and pairs[-1][0] == v:
            count = pairs[-1][1] + count
            if count == 0:
                pairs.pop()
            else:
                pairs[-1] = (v, count)
        elif count != 0:
            pairs.append((v, count))
    return pairs


def _zero_like(expr, cls):
    S = _symox()
    from ... import _core
    if cls is not ArrayDerivative:
        return S.S.Zero
    if isinstance(expr, _core._MatrixExpr):
        return S.ZeroMatrix(*expr.shape)
    return type(expr).zeros(*expr.shape)


def _evaluate(expr, pairs, cls):
    """Derivative.__new__(cls, expr, *pairs, evaluate=True) with scalar
    Symbol variables: the zero exit, the canonical order, the integer
    orders evaluated, then the held rest from the first symbolic order."""
    S = _symox()
    pairs = [(S.sympify(v), S.sympify(c)) for v, c in pairs]
    if not pairs:
        return expr
    free = expr.free_symbols
    if any(c.is_positive and v not in free for v, c in pairs):
        return _zero_like(expr, cls)
    # Derivative._sort_variable_count: commuting symbols in canonical order,
    # equal neighbours merged.
    ordered = []
    for v, c in sorted(pairs, key=lambda pair: S.default_sort_key(pair[0])):
        if ordered and ordered[-1][0] == v:
            ordered[-1] = (v, ordered[-1][1] + c)
        else:
            ordered.append((v, c))
    for i, (v, c) in enumerate(ordered):
        if not getattr(c, 'is_Integer', False):
            return cls(expr, *ordered[i:])
        expr = expr.diff((v, c))
    return expr


def _order_hint(v):
    """Cheap filter: a ``(v, count)`` pair whose count is not a nonnegative
    integer (symbolic, negative or non-integer) needs ``_parse``."""
    if not (isinstance(v, tuple) or getattr(v, 'is_Tuple', False)) or len(v) != 2:
        return False
    count = v[1]
    if type(count) is int:
        return count < 0
    return not (getattr(count, 'is_Integer', False) and not count.is_negative)


def symbolic_order_derivative(expr, variables, cls=ArrayDerivative, evaluate=True):
    """SymPy's result of ``expr.diff(*variables)`` (``cls=ArrayDerivative``)
    or ``Derivative(expr, *variables, evaluate=...)`` (``cls=Derivative``)
    for a matrix expression, matrix or array derived by scalar Symbols when
    an order is symbolic (or, for ``Derivative``, ``expr`` is explicit).
    None when another path applies."""
    from ... import _core
    hint = any(map(_order_hint, variables))
    if not hint and cls is ArrayDerivative:
        return None
    explicit = _is_explicit(expr)
    if not (explicit or (hint and isinstance(expr, _core._MatrixExpr))):
        return None
    pairs = _parse(variables)
    if pairs is None or (not _symbolic_orders(pairs) and not (explicit and cls is Derivative)):
        return None
    mutable = _is_mutable_matrix(expr)
    if explicit:
        expr = expr.as_immutable()
    if not pairs:
        result = expr
    elif cls is Derivative and not evaluate:
        return Derivative(expr, *pairs)
    else:
        result = _evaluate(expr, pairs, cls)
    if mutable and cls is ArrayDerivative:
        # MatrixBase.diff: ArrayDerivative(...).as_mutable() for a mutable matrix.
        if isinstance(result, ArrayDerivative):
            raise AttributeError("'ArrayDerivative' object has no attribute 'as_mutable'")
        if hasattr(result, 'as_mutable'):
            result = result.as_mutable()
    return result


__all__ = ['ArrayDerivative']
