"""Public ``Expr`` class: SymPy's expression category over native nodes.

Native expressions, Boolean values, relations, Tuples and the other Basic
nodes share one native storage type, so that type cannot stand for SymPy's
``Expr``: ``isinstance(false, Expr)`` and ``isinstance(Eq(x, 1), Expr)`` are
False in SymPy. The public class below is a subclass of the storage type
whose instance and subclass checks follow the nominal class of each value,
as SymPy's ``isinstance`` follows its class hierarchy.

User classes deriving from the public class derive from the storage type
instead (see ``_compat._ExprType.__new__``): their instances keep the native
layout and the metaclass the user class would have had in symox before.
"""
import weakref

from . import _core, _node_types
from ._compat import _ExprType


class Expr(_core.Expr, metaclass=_ExprType):
    """Base class for algebraic expressions."""
    __slots__ = ()
    __module__ = 'symox.core.expr'
    _symox_expr_facade = True


_nominal = []
_type_results = weakref.WeakKeyDictionary()


def _nominal_expr():
    if not _nominal:
        _nominal.append(_node_types.node_class('Expr'))
    return _nominal[0]


def _symox_class_is_expr(cls):
    # A symox class stands for the SymPy class of its name. The storage
    # types (the native Expr type, HostTree) stand for no SymPy class: they
    # also hold Boolean and Basic nodes.
    if cls is _core.Expr:
        return False
    if type.__subclasscheck__(_core._MatrixExpr, cls):
        return True
    name = type.__dict__['__dict__'].__get__(cls).get('_symox_nominal', cls.__name__)
    nominal = _node_types.node_class(name) if name in _node_types.CLASSES else None
    return nominal is not None and issubclass(nominal, _nominal_expr())


def _user_class_is_expr(cls):
    # A class defined outside symox is an expression class when one of its
    # symox ancestors is, or when it derives from the public Expr class,
    # whose place the native storage type takes in its bases.
    for base in cls.__mro__:
        if base is object:
            continue
        if _node_types.is_symox_class(base):
            if _symox_class_is_expr(base):
                return True
        elif _core.Expr in base.__bases__:
            return True
    return False


def _compute(cls):
    if _node_types.is_symox_class(cls):
        # A Function subclass defined in a symox module (such as the ODE
        # solvers' diffx) is an expression class, as a user one is.
        from ._host import _HostFunction
        if type.__subclasscheck__(_HostFunction, cls):
            return True
        return _symox_class_is_expr(cls)
    return _user_class_is_expr(cls)


def is_expr_class(cls):
    """Whether ``cls`` is a subclass of SymPy's ``Expr``."""
    if isinstance(cls, _node_types.NodeType):
        return issubclass(cls, _nominal_expr())
    try:
        return _type_results[cls]
    except KeyError:
        pass
    except TypeError:
        return _compute(cls)
    result = _type_results[cls] = _compute(cls)
    return result


def is_expr(value):
    """Whether ``value`` is an instance of SymPy's ``Expr``."""
    cls = type(value)
    if isinstance(cls, _node_types.NodeType):
        return issubclass(cls, _nominal_expr())
    if cls is _core.Expr:
        # A native node without a nominal class (Trace, Determinant,
        # DotProduct: scalar matrix expressions) is an expression unless it
        # is a Boolean value.
        name = value.class_name
        if name in _node_types.CLASSES:
            return issubclass(_node_types.node_class(name), _nominal_expr())
        return not getattr(value, 'is_Boolean', False)
    if cls is _core.Matrix or cls in _native_dense_classes():
        # The immutable matrices are matrix expressions in SymPy.
        return bool(value._immutable)
    if cls is _Function[0]:
        # Undefined function heads are classes in SymPy.
        return False
    if cls is _CallableExpr[0]:
        # Lambda nodes are expressions, applied predicates are not.
        nominal = _node_types.node_class(value.class_name)
        return nominal is not None and issubclass(nominal, _nominal_expr())
    return is_expr_class(cls)


_CallableExpr = [None]
_Function = [None]


def _install():
    from . import _CallableExpr as callable_class, Function
    _CallableExpr[0] = callable_class
    _Function[0] = Function


def _native_dense_classes():
    """The native dense matrix class and its two public SymPy-named subclasses."""
    import symox
    return (symox.MutableDenseMatrix, symox.ImmutableDenseMatrix)
