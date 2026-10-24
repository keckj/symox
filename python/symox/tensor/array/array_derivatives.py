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
                count = S.Integer(count)
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
        shape = ()
        for v, c in self._variable_count:
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
        result = expr
        for v, c in self._variable_count:
            result = result.diff((v, c))
        return result

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
            total += int(c)
            if c == 1:
                tex += r'%s %s' % (diff_symbol, printer._print(v))
            else:
                tex += r'%s %s^{%s}' % (diff_symbol, printer._print(v), printer._print(c))
        top = diff_symbol if total == 1 else r'%s^{%s}' % (diff_symbol, total)
        body = printer._print(self._expr)
        if getattr(self._expr, '_kind_tag', None) == 'MatAdd':
            body = r'\left(%s\right)' % body
        return r'\frac{%s}{%s} %s' % (top, tex, body)

    def srepr(self):
        S = _symox()
        return 'ArrayDerivative(%s)' % ', '.join(S.srepr(arg) for arg in self.args)


class Derivative(ArrayDerivative):
    """SymPy's ``Derivative(X, A)`` of a matrix expression by matrices it does not contain.

    ``Derivative.doit`` reaches ``MatrixExpr._eval_derivative``, which gives
    ``ZeroMatrix(*X.shape)`` for every such variable.
    """
    __slots__ = ()

    def doit(self, **hints):
        S = _symox()
        return S.ZeroMatrix(*self._expr.shape)

    def srepr(self):
        S = _symox()
        return 'Derivative(%s)' % ', '.join(S.srepr(arg) for arg in self.args)


def held_matrix_derivative(expr, variables, cls=ArrayDerivative):
    """The held SymPy result of ``expr.diff(*variables)``, or None.

    ``expr`` is a matrix expression and every variable a MatrixSymbol (or a
    ``(MatrixSymbol, count)`` pair) that ``expr`` does not contain.
    """
    from ... import _core
    pairs = []
    for v in variables:
        count = 1
        if isinstance(v, tuple):
            if len(v) != 2:
                return None
            v, count = v
        if not (isinstance(v, _core._MatrixExpr) and v._kind_tag == 'MatrixSymbol'):
            return None
        if not (isinstance(count, int) or getattr(count, 'is_Integer', False)) or int(count) <= 0:
            return None
        pairs.append((v, count))
    if not pairs or not isinstance(expr, _core._MatrixExpr):
        return None
    if any(expr.has(v) for v, _ in pairs):
        return None
    return cls(expr, *pairs)


__all__ = ['ArrayDerivative']
