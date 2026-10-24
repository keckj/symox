"""SymPy's derived relations, Boolean atom attributes and class flags on
native nodes.

Native relations and Boolean atoms are instances of the nominal classes of
``_node_types``. The properties below are installed on those classes and
follow sympy/core/relational.py and sympy/logic/boolalg.py: ``reversed``,
``reversedsign`` and ``negated`` build their relation without evaluating it
(``Relational.__new__`` of a subclass is ``Basic.__new__``), ``weak``,
``strict`` and ``canonical`` call the evaluating constructors.
"""
from . import _core, _node_types

_REVERSED = {'==': '==', '!=': '!=', '<': '>', '<=': '>=', '>': '<', '>=': '<='}
_NEGATED = {'==': '!=', '!=': '==', '<': '>=', '<=': '>', '>': '<=', '>=': '<'}
_CONSTRUCTORS = {'==': 'Eq', '!=': 'Ne', '<': 'Lt', '<=': 'Le', '>': 'Gt', '>=': 'Ge'}


def _unevaluated(op, lhs, rhs):
    return _core.Rel(op, lhs, rhs, evaluate=False)


def _constructor(op):
    import symox
    return getattr(symox, _CONSTRUCTORS[op])


def _is_boolean_atom(value):
    return isinstance(value, _node_types.node_class('BooleanAtom'))


def reversed(self):
    """Return the relationship with sides reversed.

    >>> from symox import Eq
    >>> from symox.abc import x
    >>> Eq(x, 1).reversed
    Eq(1, x)
    >>> (x < 1).reversed
    1 > x
    """
    a, b = self.args
    return _unevaluated(_REVERSED[self.rel_op], b, a)


def reversedsign(self):
    """Return the relationship with signs reversed.

    >>> from symox import Eq
    >>> from symox.abc import x
    >>> Eq(x, 1).reversedsign
    Eq(-x, -1)
    >>> (x < 1).reversedsign
    -x > -1
    """
    a, b = self.args
    if _is_boolean_atom(a) or _is_boolean_atom(b):
        return self
    return _unevaluated(_REVERSED[self.rel_op], -a, -b)


def negated(self):
    """Return the negated relationship, unaffected by ``evaluate=False``.

    >>> from symox import Eq
    >>> from symox.abc import x
    >>> Eq(x, 1).negated
    Ne(x, 1)
    >>> (x < 1).negated
    x >= 1
    """
    return _unevaluated(_NEGATED[self.rel_op], *self.args)


def _self(self):
    return self


def _strict(op):
    def strict(self):
        return _constructor(op)(*self.args)
    return strict


def canonical(self):
    """Return a canonical form of the relational: a number on the rhs, a
    sign removed from the lhs, or else canonically ordered args.

    >>> from symox.abc import x, y
    >>> (x < 2).reversed.canonical
    x < 2
    >>> (-y < x).canonical
    x > -y
    >>> (-y > x).canonical
    x < -y
    >>> (-y < -x).canonical
    x < y
    """
    from symox.core.sorting import ordered
    relational = _node_types.node_class('Relational')
    args = tuple([i.canonical if isinstance(i, relational) else i for i in self.args])
    if args != self.args:
        r = _constructor(self.rel_op)(*args)
        if not isinstance(r, relational):
            return r
    else:
        r = self
    if r.rhs.is_number:
        if r.rhs.is_Number and r.lhs.is_Number and r.lhs > r.rhs:
            r = r.reversed
    elif r.lhs.is_number:
        r = r.reversed
    elif tuple(ordered(args)) != args:
        r = r.reversed

    LHS_CEMS = getattr(r.lhs, 'could_extract_minus_sign', None)
    RHS_CEMS = getattr(r.rhs, 'could_extract_minus_sign', None)

    if _is_boolean_atom(r.lhs) or _is_boolean_atom(r.rhs):
        return r

    # Check if first value has negative sign
    if LHS_CEMS and LHS_CEMS():
        return r.reversedsign
    elif not r.rhs.is_number and RHS_CEMS and RHS_CEMS():
        # Right hand side has a minus, but not lhs: Eq(x, -y) and Eq(-x, y)
        # have the same canonical representation.
        expr1, _ = ordered([r.lhs, -r.rhs])
        if expr1 != r.lhs:
            return r.reversed.reversedsign
    return r


def _first(self):
    return self.args[0]


def _second(self):
    return self.args[1]


def _side(index):
    def side(self):
        value = self.args[index]
        # An explicit matrix side is stored as a native ImmutableMatrix
        # node; SymPy's lhs/rhs is the ImmutableDenseMatrix itself.
        if getattr(value, '_kind_tag', None) == 'ImmutableMatrix':
            return value.as_explicit()
        return value
    return side


def _set(name, attribute, value):
    type.__setattr__(_node_types.node_class(name), attribute, value)


# Class flags that SymPy defines on these classes (sympy/core/basic.py,
# expr.py, numbers.py, symbol.py, logic/boolalg.py, tensor/indexed.py) and
# that no native getter provides.
_FLAGS = (
    ('Basic', ('is_AlgebraicNumber', 'is_Indexed', 'is_MatAdd', 'is_MatMul', 'is_Not', 'is_NumberSymbol',
               'is_Point', 'is_Poly', 'is_Vector', 'is_symbol'), False),
    ('Expr', ('is_scalar',), True),
    ('AlgebraicNumber', ('is_AlgebraicNumber',), True),
    ('Indexed', ('is_Indexed', 'is_symbol'), True),
    ('IndexedBase', ('is_symbol',), True),
    ('Idx', ('is_symbol',), True),
    ('Not', ('is_Not',), True),
    ('NumberSymbol', ('is_NumberSymbol',), True),
    ('Symbol', ('is_symbol',), True),
)


def install():
    for name, attributes, value in _FLAGS:
        if _node_types.node_class(name) is not None:
            for attribute in attributes:
                _set(name, attribute, value)
    for attribute, function in (('reversed', reversed), ('reversedsign', reversedsign),
                                ('negated', negated), ('canonical', canonical),
                                ('weak', _self), ('strict', _self)):
        _set('Relational', attribute, property(function, doc=function.__doc__))
    _set('Relational', 'lhs', property(_side(0)))
    _set('Relational', 'rhs', property(_side(1)))
    _set('_Greater', 'gts', property(_first))
    _set('_Greater', 'lts', property(_second))
    _set('_Less', 'gts', property(_second))
    _set('_Less', 'lts', property(_first))
    _set('GreaterThan', 'strict', property(_strict('>')))
    _set('LessThan', 'strict', property(_strict('<')))
    _set('StrictGreaterThan', 'weak', property(_strict('>=')))
    _set('StrictLessThan', 'weak', property(_strict('<=')))

    from . import S
    _set('BooleanAtom', '_op_priority', 11)
    _set('BooleanAtom', 'canonical', property(_self))
    _set('BooleanTrue', 'negated', property(lambda self: S.false))
    _set('BooleanFalse', 'negated', property(lambda self: S.true))
