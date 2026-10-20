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
        if cls not in (Boolean, BooleanFunction) or not isinstance(value, _core.Expr):
            return False
        if cls is BooleanFunction:
            return value.is_Boolean and value.is_Function
        return value._is_boolean_category


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


__all__ += ['Boolean', 'BooleanFunction']


from symox._compat import _ExprType


class BooleanAtom(metaclass=_ExprType):
    """Category of the native Boolean singletons ``true`` and ``false``."""
    _matches = staticmethod(lambda expr: expr.class_name in ('BooleanTrue', 'BooleanFalse'))


class BooleanTrue(BooleanAtom):
    """Class of ``S.true``; instantiation returns that singleton."""
    _matches = staticmethod(lambda expr: expr.class_name == 'BooleanTrue')

    def __new__(cls):
        from symox import true
        return true


class BooleanFalse(BooleanAtom):
    """Class of ``S.false``; instantiation returns that singleton."""
    _matches = staticmethod(lambda expr: expr.class_name == 'BooleanFalse')

    def __new__(cls):
        from symox import false
        return false


__all__ += ['BooleanAtom', 'BooleanTrue', 'BooleanFalse']
