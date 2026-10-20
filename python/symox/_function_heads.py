"""Builtin function heads as classes.

SymPy's function heads (``sin``, ``Abs``, ``Mod``...) are classes: code
registers them with ``functools.singledispatch``, tests ``issubclass`` and
derives its own classes from them (Devito's ``class Abs(DifferentiableFunction,
sympy.Abs)``). Each native head below is such a class:

* calling the head itself runs its native constructor, which evaluates in
  the kernel and returns native nodes;
* native nodes are instances of the head through their nominal class, so
  ``isinstance(asin(x), InverseTrigonometricFunction)`` and
  ``isinstance(sin(x), sin)`` hold as in SymPy;
* a class derived from a head outside symox is a user expression class (see
  ``_host.expression_class``). Constructing it runs the native evaluation as
  SymPy's ``Function.__new__`` runs ``eval``: an evaluated result is
  returned, while an unevaluated application becomes an instance of the user
  class, which differentiates as the native head does.
"""
import inspect

from . import _core, _node_types
from ._compat import _ExprType
from ._expr_facade import Expr


def _own(cls):
    return type.__dict__['__dict__'].__get__(cls)


def head_of(cls):
    """The builtin head class that ``cls`` is or derives from, or None."""
    for base in cls.__mro__:
        if isinstance(base, BuiltinHeadType) and '_symox_own_head' in _own(base):
            return base
    return None


def nominal_of(cls):
    """The nominal node class standing for a builtin head class or a class
    derived from one, or None."""
    head = head_of(cls) if isinstance(cls, type) else None
    if head is None:
        return None
    return _node_types.node_class(_own(head)['_symox_own_head'])


class _HeadName:
    """``_symox_head_name`` of a head class.

    Only the head itself is the native head: a user class derived from it
    builds its own applications, so structural queries keyed by the user
    class (``has``, ``atoms``, ``replace``) must not match native nodes.
    """

    def __get__(self, cls, owner=None):
        if cls is None:
            return self
        name = _own(cls).get('_symox_own_head')
        if name is None:
            raise AttributeError('_symox_head_name')
        return name

    def __set__(self, cls, value):
        type.__setattr__(cls, '_symox_own_head', value)


def _as_nargs(value):
    # Reference: FunctionClass.__init__ and FunctionClass.nargs.
    from . import FiniteSet, Naturals0
    if value is None:
        return Naturals0
    if isinstance(value, int):
        value = (value,)
    return FiniteSet(*value)


class _Nargs:
    """``nargs`` of a head class or of a class derived from one: the
    argument counts SymPy's class accepts, or a user class's own ``nargs``."""

    def __get__(self, cls, owner=None):
        if cls is None:
            return self
        from ._function_nargs import NARGS
        for base in cls.__mro__:
            namespace = _own(base)
            if '_symox_own_head' in namespace:
                return _as_nargs(NARGS.get(namespace['_symox_own_head']))
            if 'nargs' in namespace and not _node_types.is_symox_class(base):
                value = namespace['nargs']
                return value if hasattr(value, 'is_Set') else _as_nargs(value)
        raise AttributeError('nargs')

    def __set__(self, cls, value):
        type.__setattr__(cls, 'nargs', value)


class BuiltinHeadType(_ExprType):
    """Metaclass of the builtin function heads (SymPy's FunctionClass)."""

    _symox_head_name = _HeadName()
    nargs = _Nargs()

    def __new__(metaclass, name, *rest, **kwargs):
        if len(rest) == 1 and not kwargs:
            return _head_class(metaclass, name, rest[0])
        return super().__new__(metaclass, name, *rest, **kwargs)

    def __init__(cls, *args, **kwargs):
        pass

    def __call__(cls, *args, **kwargs):
        factory = _own(cls).get('_symox_factory')
        if factory is not None:
            return factory(*args, **kwargs)
        return super().__call__(*args, **kwargs)

    def __getattr__(cls, name):
        try:
            return super().__getattr__(name)
        except AttributeError:
            pass
        factory = _own(cls).get('_symox_factory')
        if factory is None or name.startswith('__'):
            raise AttributeError(name)
        return getattr(factory, name)

    def __instancecheck__(cls, value):
        if type.__instancecheck__(cls, value):
            return True
        name = _own(cls).get('_symox_own_head')
        if name is None or not isinstance(value, _core.Expr):
            return False
        kind = type(value)
        if isinstance(kind, _node_types.NodeType):
            nominal = _node_types.node_class(name)
            if nominal is not None:
                return issubclass(kind, nominal)
        return kind is _core.Expr and value.func is cls

    def __repr__(cls):
        if '_symox_own_head' in _own(cls):
            return cls.__name__
        return super().__repr__()

    __str__ = __repr__


def _head_new(cls, *args, **options):
    head = head_of(cls)
    if cls is head:
        return _own(head)['_symox_factory'](*args, **options)
    return _subclass_instance(cls, head, args, options)


def _user_eval(cls):
    from ._host import _user_classes
    for base in _user_classes(cls):
        if 'eval' in base.__dict__:
            return getattr(cls, 'eval')
    return None


def _subclass_instance(cls, head, args, options):
    # Reference: Application.__new__, with the native constructor of the
    # head in place of the inherited ``eval``.
    from . import sympify
    from ._public_core import _evaluation
    evaluate = options.pop('evaluate', None)
    if evaluate is None:
        evaluate = _evaluation.get()
    options.pop('nargs', None)
    if options:
        raise ValueError("Unknown options: %s" % options)
    args = tuple(sympify(arg) for arg in args)
    factory = _own(head)['_symox_factory']
    if evaluate:
        evaluator = _user_eval(cls)
        if evaluator is not None:
            evaluated = evaluator(*args)
            if evaluated is not None:
                return evaluated
        result = factory(*args)
        if not _unevaluated(result, head, args):
            return result
    return _core.Expr.__new__(cls, *args)


def _unevaluated(result, head, args):
    return (isinstance(result, _core.Expr) and result.class_name == _own(head)['_symox_own_head']
            and not _core._host_contains(result)[1] and tuple(result.args) == args)


def _twin_derivative(self, symbol):
    # SymPy's Function._eval_derivative over the inherited fdiff: the
    # derivative of the native application with the same arguments.
    head = head_of(type(self))
    try:
        twin = _own(head)['_symox_factory'](*self.args, evaluate=False)
    except TypeError:
        twin = _core._apply(_own(head)['_symox_own_head'], *self.args, evaluate=False)
    # The head's fdiff returns the application itself where SymPy's does
    # (exp, unevaluated derivatives): that application is ``self``.
    return twin.diff(symbol).xreplace({twin: self})


def _defines(cls, name):
    from ._host import _defines as defines
    return defines(cls, name)


def _subclass_hook(head):
    def __init_subclass__(cls, **kwargs):
        # User applications are SymPy Function instances that differentiate
        # by the chain rule of the head they derive from.
        if not _defines(cls, 'is_Function'):
            type.__setattr__(cls, 'is_Function', True)
        if not _defines(cls, '_eval_derivative'):
            type.__setattr__(cls, '_eval_derivative', _twin_derivative)
        # Rebuilt applications (subs, xreplace, func(*args)) evaluate again,
        # as SymPy's Function.__new__ runs eval on every construction.
        if not _defines(cls, '__new__'):
            type.__setattr__(cls, '__new__', _head_new)
        super(head, cls).__init_subclass__(**kwargs)
    return classmethod(__init_subclass__)


def _head_class(metaclass, name, factory):
    namespace = {
        '__module__': 'symox',
        '__qualname__': name,
        '__doc__': factory.__doc__,
        '__slots__': (),
        '__new__': _head_new,
        '_symox_own_head': name,
        '_symox_factory': factory,
        '_symox_expression_kind': 'expr',
    }
    try:
        namespace['__signature__'] = inspect.signature(factory)
    except (TypeError, ValueError):
        pass
    if isinstance(factory, type):
        # A head over a Python class: the class's own documentation.
        namespace['__doc__'] = type.__dict__['__doc__'].__get__(factory)
    cls = _ExprType.__new__(metaclass, name, _head_bases(name), namespace)
    type.__setattr__(cls, '__init_subclass__', _subclass_hook(cls))
    return cls


class FunctionHead:
    """Stable public callable identity for native heads that SymPy defines
    as plain functions (``sqrt``, ``cbrt``, ``E1``)."""

    def __init__(self, name, factory):
        self.__name__ = self.__qualname__ = self._symox_head_name = name
        self.__doc__ = factory.__doc__
        self.__module__ = factory.__module__
        self.__signature__ = inspect.signature(factory)
        self._factory = factory

    def __call__(self, *args, **kwargs):
        return self._factory(*args, **kwargs)

    def __getattr__(self, name):
        if name == '_factory':
            raise AttributeError(name)
        try:
            return _node_types.class_attribute(self.__name__, name)
        except AttributeError:
            return getattr(self._factory, name)

    def __setattr__(self, name, value):
        object.__setattr__(self, name, value)
        if '_factory' in self.__dict__:
            _node_types.forward_setattr(self, self.__name__, name, value)

    def __delattr__(self, name):
        object.__delattr__(self, name)
        _node_types.forward_delattr(self, self.__name__, name)

    def __repr__(self):
        return self.__name__

    __str__ = __repr__

    def __instancecheck__(self, value):
        return isinstance(value, _core.Expr) and value.func is self

    def __reduce__(self):
        from . import _public_function_head
        return _public_function_head, (self.__name__,)


def builtin_head(name, factory):
    """The public head of the native function ``name``: a class when SymPy
    defines one of that name."""
    from ._sympy_classes import CLASSES
    if name in CLASSES:
        return BuiltinHeadType(name, factory)
    return FunctionHead(name, factory)


def is_builtin_head(value):
    return isinstance(value, FunctionHead) or (
        isinstance(value, BuiltinHeadType) and '_symox_own_head' in _own(value))



# --- SymPy's function class hierarchy -----------------------------------------

class Application(metaclass=_ExprType):
    """Base class of applied functions (sympy.core.function.Application)."""
    __slots__ = ()
    __module__ = 'symox.core.function'
    _matches = staticmethod(lambda e: e.is_Function)

    def __new__(cls, *args):
        if cls is Application:
            raise TypeError('instantiate a concrete function application')
        from . import Function
        return Function(cls.__name__)(*args)


class DefinedFunction(Application):
    """Base class of the defined mathematical functions
    (sympy.core.function.DefinedFunction)."""
    __slots__ = ()
    __module__ = 'symox.core.function'
    _matches = staticmethod(lambda e: e.is_Function and not e.is_UndefinedFunction)

    def __new__(cls, *args):
        from . import Function, sympify
        values = tuple(sympify(a) for a in args)
        evaluator = getattr(cls, 'eval', None)
        if evaluator:
            result = evaluator(*values)
            if result is not None:
                return sympify(result)
        return Function(cls.__name__)(*values)


_PUBLIC_ROOTS = {'Application': Application, 'DefinedFunction': DefinedFunction}
_abstract = {}


def _in_function_hierarchy(name):
    from ._sympy_classes import CLASSES
    if name in _PUBLIC_ROOTS:
        return True
    bases = CLASSES[name][1] if name in CLASSES else ()
    return bool(bases) and all(_in_function_hierarchy(base) for base in bases)


def _head_bases(name):
    # The head derives from the public classes of its SymPy bases, which
    # end at DefinedFunction or Application. Heads whose SymPy class is
    # also an operation (Max, Min: lattice operations, defined in
    # functions.elementary.miscellaneous) keep the expression base alone.
    from ._sympy_classes import CLASSES
    if not _in_function_hierarchy(name):
        return (Expr,)
    return tuple(abstract_class(base) for base in CLASSES[name][1]) + (Expr,)


def abstract_class(name):
    """The public class of SymPy's abstract function class ``name`` (such as
    TrigonometricFunction).

    It has no native constructor: native nodes are its instances through
    their nominal classes, and the head classes derive from it as in SymPy.
    """
    cls = _PUBLIC_ROOTS.get(name) or _abstract.get(name)
    if cls is None:
        from ._sympy_classes import CLASSES
        module, bases = CLASSES[name]
        cls = _abstract[name] = _ExprType(name, tuple(abstract_class(base) for base in bases), {
            '__module__': 'symox' + module[len('sympy'):], '__qualname__': name,
            '__slots__': (), '_matches': None})
    return cls


def define_abstract(namespace, names):
    """Export SymPy's abstract function classes ``names`` from the module
    ``namespace``."""
    exported = namespace.setdefault('__all__', [])
    for name in names.split():
        namespace[name] = abstract_class(name)
        if name not in exported:
            exported.append(name)
