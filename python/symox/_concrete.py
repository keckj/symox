"""SymPy ``ExprWithIntLimits`` protocol for held ``Sum``/``Product`` nodes.

Convergence tests and limit manipulations are native; this module only
resolves symbolic indices, attaches the methods to ``Expr`` and raises
SymPy's ``ReorderError``.
"""
from . import _core
from ._core import Expr


class ReorderError(NotImplementedError):
    """Exception raised when trying to reorder dependent limits."""
    def __init__(self, expr, msg):
        super().__init__("%s could not be reordered: %s." % (expr, msg))


def _require(self, name, *heads):
    if not any(getattr(self, 'is_' + h) for h in heads):
        raise AttributeError("'%s' object has no attribute '%s'" % (type(self).__name__, name))


def _is_convergent(self):
    """Convergence of a one-limit Sum or Product."""
    return _core._is_convergent(self)


def _is_absolutely_convergent(self):
    """``Sum(Abs(f)).is_convergent()`` for a one-limit Sum."""
    return _core._is_absolutely_convergent(self)


def _index(self, x):
    _require(self, 'index', 'Sum', 'Product')
    return _core._limits_index(self, x)


def _resolve(self, indices):
    return [i if isinstance(i, int) else self.index(i) for i in indices]


def _reverse_order(self, *indices):
    _require(self, 'reverse_order', 'Sum', 'Product')
    return _core._reverse_order(self, _resolve(self, indices))


def _change_index(self, var, trafo, newvar=None):
    _require(self, 'change_index', 'Sum', 'Product')
    return _core._change_index(self, var, trafo, newvar)


def _reorder_limit(self, x, y):
    _require(self, 'reorder_limit', 'Sum', 'Product')
    try:
        return _core._reorder_limit(self, x, y)
    except NotImplementedError as error:
        raise ReorderError(self, str(error)) from None


def _reorder(self, *arg):
    _require(self, 'reorder', 'Sum', 'Product')
    new_expr = self
    for r in arg:
        if len(r) != 2:
            raise ValueError(r, "Invalid number of arguments")
        index1, index2 = _resolve(self, r)
        new_expr = new_expr.reorder_limit(index1, index2)
    return new_expr


def _limits_property(name, *heads):
    def getter(self):
        _require(self, name, *heads)
        return _core._limits_property(self, name)
    getter.__name__ = name
    return property(getter)


def _term(self):
    _require(self, 'term', 'Product')
    return self.args[0]


Expr.is_convergent = _is_convergent
Expr.is_absolutely_convergent = _is_absolutely_convergent
Expr.index = _index
Expr.reverse_order = _reverse_order
Expr.change_index = _change_index
Expr.reorder_limit = _reorder_limit
Expr.reorder = _reorder
Expr.has_finite_limits = _limits_property('has_finite_limits', 'Sum', 'Product', 'Integral')
Expr.has_reversed_limits = _limits_property('has_reversed_limits', 'Sum', 'Product', 'Integral')
Expr.has_empty_sequence = _limits_property('has_empty_sequence', 'Sum', 'Product')
Expr.term = property(_term)


# Indexed objects enter native arithmetic as symbols named like ``f[n]``
# (codegen.ast.Indexed._symox_), which hides their dependence on the index.
# The symbols are recorded by _public_matrix, and the evaluation of held sums and products
# replaces those whose indices are not constant by applications of an
# undefined function of the indices, so the summation engine substitutes the
# summation variable into them as SymPy substitutes into Indexed.
_INDEXED_HEADS = {}
_REGISTRY = []


def _indexed_mapping(expr):
    if not _REGISTRY:
        from ._public_matrix import _indexed
        _REGISTRY.append(_indexed)
    _INDEXED_SYMBOLS = _REGISTRY[0]
    if not _INDEXED_SYMBOLS:
        return None
    if type(expr).__name__ == 'Indexed' and hasattr(expr, '_symox_'):
        expr = expr._symox_()
    if not isinstance(expr, Expr):
        return None
    try:
        free = expr.free_symbols
    except Exception:
        return None
    flattened = [s for s in free if s in _INDEXED_SYMBOLS]
    if not flattened:
        return None
    from . import Function
    forward, heads = {}, {}
    for symbol in flattened:
        indexed = _INDEXED_SYMBOLS[symbol]
        if not any(getattr(i, 'free_symbols', None) for i in indexed.indices):
            continue
        key = (indexed.base, len(indexed.indices))
        head = _INDEXED_HEADS.get(key)
        if head is None:
            head = _INDEXED_HEADS[key] = Function('_symox_indexed_%d' % len(_INDEXED_HEADS))
        heads[head] = key[0]
        forward[symbol] = head(*indexed.indices)
    return (forward, heads, expr) if forward else None


def _indexed_restore(value, heads):
    from .codegen.ast import Indexed
    if not isinstance(value, Expr):
        return value
    for head, base in heads.items():
        value = value.replace(head, lambda *indices, _base=base: Indexed(_base, *indices)._symox_())
    return value


def _through_indexed(native):
    def call(expr, *args, **kwargs):
        mapping = _indexed_mapping(expr)
        if mapping is None:
            return native(expr, *args, **kwargs)
        forward, heads, expr = mapping
        return _indexed_restore(native(expr.xreplace(forward), *args, **kwargs), heads)
    call.__name__ = getattr(native, '__name__', 'call')
    call.__doc__ = getattr(native, '__doc__', None)
    call._symox_native = native
    return call


Expr.doit = _through_indexed(Expr.doit)
summation = _through_indexed(_core.summation)
product = _through_indexed(_core.product)
