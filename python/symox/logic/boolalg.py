"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox', 'And Or Not Xor Nand Nor Implies Equivalent ITE true false to_cnf to_dnf to_nnf POSform SOPform simplify_logic bool_map')


def gateinputcount(expr):
    """Count wires entering Boolean gates, including nested gate occurrences."""
    from symox import sympify
    expr=sympify(expr)
    gates={'And','Or','Not','Xor','Nand','Nor','ITE','Implies','Equivalent'}
    if expr.class_name not in gates:return 0
    return len(expr.args)+sum(gateinputcount(a) for a in expr.args)
__all__ += ['gateinputcount']


from symox._basic import Basic, _BasicMeta


class _BooleanType(_BasicMeta):
    def __instancecheck__(cls, value):
        if type.__instancecheck__(cls, value):
            return True
        from symox import _core
        if not isinstance(value, _core.Expr):
            return False
        native_names = cls.__dict__.get('_native_names')
        if native_names is not None:
            return value.class_name in native_names
        if cls not in (Boolean, BooleanFunction):
            return False
        if cls is BooleanFunction:
            return value.is_Boolean and value.is_Function
        return value._is_boolean_category

    # The Boolean atoms are the nominal classes of the native singletons:
    # class dispatch on them accepts those node types, and attributes patched
    # onto them after import (SageMath's ``_sage_``) reach ``true``/``false``.
    def _native_counterpart(cls):
        if '_native_names' not in cls.__dict__ or globals().get(cls.__name__) is not cls:
            return None
        from symox._node_types import node_class
        return node_class(cls.__name__)

    def __subclasscheck__(cls, subclass):
        if super().__subclasscheck__(subclass):
            return True
        from symox._node_types import NodeType
        nominal = cls._native_counterpart() if isinstance(subclass, NodeType) else None
        return nominal is not None and issubclass(subclass, nominal)

    def __setattr__(cls, name, value):
        super().__setattr__(name, value)
        if cls._native_counterpart() is not None:
            from symox._node_types import forward_setattr
            forward_setattr(cls, cls.__name__, name, value)

    def __delattr__(cls, name):
        super().__delattr__(name)
        if cls._native_counterpart() is not None:
            from symox._node_types import forward_delattr
            forward_delattr(cls, cls.__name__, name)


class Boolean(Basic, metaclass=_BooleanType):
    """Boolean category and formal Basic container.

    Native gate, relation and Boolean-atom operations remain in the native
    kernel. Algebra involving a directly constructed formal base container
    is not yet implemented by that kernel.
    """
    __slots__ = ()
    is_Boolean = False
    is_Function = False
    is_Atom = False
    is_Symbol = False
    is_Number = False
    is_Relational = False



class BooleanFunction(Boolean):
    """Category for native gates and formal Boolean function applications."""
    is_Boolean = True
    is_Function = True

    def __new__(cls, *args):
        from symox import sympify
        return Boolean.__new__(cls, *(sympify(arg) for arg in args))


class BooleanAtom(Boolean):
    """Base class of :py:class:`~.BooleanTrue` and :py:class:`~.BooleanFalse`.

    The two atoms are the native singletons ``true`` and ``false``.
    """
    _native_names = ('BooleanTrue', 'BooleanFalse')
    is_Boolean = True
    is_Atom = True
    _op_priority = 11  # higher than Expr


class BooleanTrue(BooleanAtom):
    """SymPy version of ``True``, a singleton that can be accessed via ``S.true``."""
    _native_names = ('BooleanTrue',)

    def __new__(cls):
        from symox import true
        return true


class BooleanFalse(BooleanAtom):
    """SymPy version of ``False``, a singleton that can be accessed via ``S.false``."""
    _native_names = ('BooleanFalse',)

    def __new__(cls):
        from symox import false
        return false


__all__ += ['Boolean', 'BooleanFunction', 'BooleanAtom', 'BooleanTrue', 'BooleanFalse']
