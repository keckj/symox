"""``Order.contains`` and ``Order.__contains__``, ported from SymPy 1.14
``sympy/series/order.py``.

The native ``_order_contains`` only answers True or False; the reference
answers None when the inclusion cannot be decided (different symbols) and
compares the ``(variable, point)`` pairs when both orders have the same
expression, so ``O(1) in O(1, x)`` holds.
"""

from symox import S, sympify


def _pairs(order):
    return tuple(zip(order.variables, order.point))


def _rebuild(order, expr):
    from symox import Order
    return Order(expr, *_pairs(order))


def _contains(self, expr):
    expr = sympify(expr)
    if expr.is_zero:
        return True
    if expr is S.NaN:
        return False
    point = self.point[0] if self.point else S.Zero
    if expr.is_Order:
        if (any(p != point for p in expr.point) or
                any(p != point for p in self.point)):
            return None
        if expr.expr == self.expr:
            # O(1) + O(1), O(1) + O(1, x), etc.
            mine = _pairs(self)
            return all(x in mine for x in _pairs(expr))
        if expr.expr.is_Add:
            return all(self.contains(x) for x in expr.expr.args)
        if self.expr.is_Add and point.is_zero:
            return any(_rebuild(self, x).contains(expr)
                       for x in self.expr.args)
        if self.variables and expr.variables:
            common_symbols = tuple(
                [s for s in self.variables if s in expr.variables])
        elif self.variables:
            common_symbols = self.variables
        else:
            common_symbols = expr.variables
        if not common_symbols:
            return None
        if (self.expr.is_Pow and len(self.variables) == 1
                and self.variables == expr.variables):
            symbol = self.variables[0]
            other = expr.expr.as_independent(symbol, as_Add=False)[1]
            if (other.is_Pow and other.base == symbol and
                    self.expr.base == symbol):
                rv = None
                if point.is_zero:
                    rv = (self.expr.exp - other.exp).is_nonpositive
                if point.is_infinite:
                    rv = (self.expr.exp - other.exp).is_nonnegative
                if rv is not None:
                    return rv
        from symox.simplify.powsimp import powsimp
        from symox import Limit
        r = None
        ratio = self.expr/expr.expr
        ratio = powsimp(ratio, deep=True, combine='exp')
        for s in common_symbols:
            # The native Limit.doit takes no heuristics hint; its answer is
            # the gruntz value whenever that one exists.
            try:
                l = Limit(ratio, s, point).doit()
            except NotImplementedError:
                l = Limit(ratio, s, point)
            if not isinstance(l, Limit):
                l = l != 0
            else:
                l = None
            if r is None:
                r = l
            else:
                if r != l:
                    return
        return r
    if self.expr.is_Pow and len(self.variables) == 1:
        symbol = self.variables[0]
        other = expr.as_independent(symbol, as_Add=False)[1]
        if (other.is_Pow and other.base == symbol and
                self.expr.base == symbol):
            rv = None
            if point.is_zero:
                rv = (self.expr.exp - other.exp).is_nonpositive
            if point.is_infinite:
                rv = (self.expr.exp - other.exp).is_nonnegative
            if rv is not None:
                return rv
    obj = _rebuild(self, expr)
    return self.contains(obj)


def _dunder_contains(self, other):
    result = self.contains(other)
    if result is None:
        raise TypeError('contains did not evaluate to a bool')
    return result


def install():
    from symox import _node_types
    cls = _node_types.node_class('Order')
    cls.contains = _contains
    cls.__contains__ = _dunder_contains
