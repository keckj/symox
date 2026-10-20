"""Interned scalar kinds used by typed array metadata."""


class KindMeta(type):
    def __new__(mcls,name,bases,namespace,**kwargs):
        cls=super().__new__(mcls,name,bases,namespace,**kwargs);cls._kind_instances={};return cls
    def __call__(cls,*args,**kwargs):
        key=args,tuple(sorted(kwargs.items()))
        if key not in cls._kind_instances:cls._kind_instances[key]=super().__call__(*args,**kwargs)
        return cls._kind_instances[key]


class Kind(metaclass=KindMeta):
    __slots__ = ()


class _ScalarKind(Kind):
    __slots__ = ('name',)
    def __init__(self,name):self.name=name
    def __repr__(self):return self.name
    def __reduce__(self):return _scalar_kind,(self.name,)
    def to_sympy(self):
        from sympy.core import kind
        return getattr(kind,self.name)


NumberKind = _ScalarKind('NumberKind')
BooleanKind = _ScalarKind('BooleanKind')
UndefinedKind = _ScalarKind('UndefinedKind')


def _scalar_kind(name):
    return {'NumberKind':NumberKind,'BooleanKind':BooleanKind,
            'UndefinedKind':UndefinedKind}[name]


def scalar_kind(value,_depth=0):
    from . import true,false
    if _depth>128:return UndefinedKind
    if value in (true,false) or getattr(value,'is_Relational',False) or getattr(value,'is_Boolean',False):
        return BooleanKind
    if (not value.args or value.is_Number) and value.is_commutative is True:
        return NumberKind
    if getattr(value.func,'__name__','')=='exp' and value.is_commutative is True:
        return NumberKind
    if value.is_Add or value.is_Mul or value.is_Pow:
        kinds={scalar_kind(arg,_depth+1) for arg in value.args}
        return next(iter(kinds)) if len(kinds)==1 else UndefinedKind
    if value.is_Integral or value.is_Sum or value.is_Product:
        return scalar_kind(value.function,_depth+1)
    return UndefinedKind


__all__ = ['Kind','NumberKind','BooleanKind','UndefinedKind']


# The ``kind`` attribute of native nodes (sympy/core/*.py and
# sympy/matrices/repmatrix.py). Each rule is keyed by the SymPy class that
# defines it and applies to the instances of that class and of its
# subclasses, found in the nominal class hierarchy of the node.

def _number_kind(value, kind):
    return NumberKind


def _symbol_kind(value, kind):
    return NumberKind if value.is_commutative else UndefinedKind


def _add_kind(value, kind):
    kinds = frozenset(kind(arg) for arg in value.args)
    if len(kinds) != 1:
        # Since addition is group operator, kind must be same.
        return UndefinedKind
    result, = kinds
    return result


def mul_kind(*kinds):
    """SymPy's ``Mul._kind_dispatcher`` with the rules of the core and of
    the matrices: equal kinds give that kind, numbers times matrices give
    matrices of the product of the element kinds, other pairs are
    undefined."""
    from .matrices.kind import MatrixKind
    kinds = frozenset(kinds)
    if not kinds:
        return UndefinedKind
    result = None
    for kind in kinds:
        if result is None:
            result = kind
        elif isinstance(result, MatrixKind) and isinstance(kind, MatrixKind):
            result = MatrixKind(mul_kind(result.element_kind, kind.element_kind))
        elif isinstance(result, MatrixKind) and kind is NumberKind:
            result = MatrixKind(mul_kind(kind, result.element_kind))
        elif result is NumberKind and isinstance(kind, MatrixKind):
            result = MatrixKind(mul_kind(result, kind.element_kind))
        else:
            result = UndefinedKind
    return result


def _mul_kind(value, kind):
    return mul_kind(*(kind(arg) for arg in value.args))


def _pow_kind(value, kind):
    if kind(value.exp) is NumberKind:
        return kind(value.base)
    return UndefinedKind


def _first_argument_kind(value, kind):
    return kind(value.args[0])


def _exp_kind(value, kind):
    return kind(value.exp)


def _function_kind(value, kind):
    return kind(value.function)


def _tuple_kind(value, kind):
    from .core.containers import TupleKind
    return TupleKind(*(kind(arg) for arg in value.args))


def _boolean_kind(value, kind):
    return BooleanKind


_RULES = {
    'Number': _number_kind, 'AlgebraicNumber': _number_kind,
    'ComplexInfinity': _number_kind, 'NumberSymbol': _number_kind,
    'ImaginaryUnit': _number_kind, 'Mod': _number_kind, 'Rem': _number_kind,
    'Symbol': _symbol_kind, 'Add': _add_kind, 'Mul': _mul_kind, 'Pow': _pow_kind,
    'Derivative': _first_argument_kind, 'ExpBase': _exp_kind,
    'ExprWithLimits': _function_kind, 'Tuple': _tuple_kind, 'Boolean': _boolean_kind,
}


def _rule(cls):
    from ._node_types import is_symox_class
    for base in cls.__mro__:
        rule = _RULES.get(base.__name__)
        if rule is not None and is_symox_class(base):
            return rule
    return None


def expression_kind(value):
    """The kind of a native node or of an instance of a symox class.

    Native trees share repeated subexpressions, so the kinds of the
    subexpressions found while computing one kind are remembered for that
    computation: a tree with many occurrences of a shared node is not
    unfolded. Arguments whose class defines its own ``kind`` are asked.
    """
    memo = {}
    rules = {}

    def kind(node):
        cls = type(node)
        if getattr(cls, 'kind', None) is not _PROPERTY:
            return node.kind
        try:
            return memo[node]
        except KeyError:
            pass
        except TypeError:
            return compute(node)
        result = memo[node] = compute(node)
        return result

    def compute(node):
        cls = type(node)
        try:
            rule = rules[cls]
        except KeyError:
            rule = rules[cls] = _rule(cls)
        return UndefinedKind if rule is None else rule(node, kind)

    return compute(value)


_PROPERTY = property(expression_kind)


def matrix_tag(value):
    """The string ``kind`` tag of a native matrix expression, or None.

    Scalar nodes have SymPy kinds rather than tags, so their kind, which
    walks the whole tree, is not computed by tag tests.
    """
    from . import _core
    if isinstance(value, _core.Expr) and getattr(type(value), 'kind', None) is _PROPERTY:
        return None
    return getattr(value, 'kind', None)


def matrix_kind(matrix):
    """The kind of an explicit matrix, as SymPy's ``RepMatrix.kind``: the
    matrices of rationals have numbers, the others the common kind of
    their nonzero entries."""
    from .matrices.kind import MatrixKind
    entries = list(matrix)
    if all(entry.is_Rational for entry in entries):
        return MatrixKind(NumberKind)
    kinds = {entry.kind for entry in entries if entry != 0}
    return MatrixKind(kinds.pop() if len(kinds) == 1 else UndefinedKind)


def install():
    from . import _core
    _core.Expr.kind = _PROPERTY
    _core.Matrix.kind = property(matrix_kind)
