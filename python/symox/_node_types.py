"""SymPy's nominal classes for native expression nodes.

Every native scalar node is an instance of a heap subtype of ``_core.Expr``
named after the SymPy class of that node (``Add``, ``Integer``, ``sin``...).
The subtypes reproduce SymPy's module, name and base-class graph from the
generated table in ``_sympy_classes`` and carry no state of their own.

Downstream protocols identify nodes through ``type(expr)``, ``__class__``
and the MRO, and monkeypatch attributes on SymPy classes after import (for
example SageMath's ``_sage_`` conversions). Public symox classes and heads
forward such assignments here, so normal attribute lookup on instances
follows SymPy's MRO and later patches are seen without any extra cache.
"""
import copyreg
import importlib
import numbers
import sys
import weakref

from . import _core
from ._basic import nominal_hash
from ._host import HostTree
from ._sympy_classes import CLASSES

_classes = {}


def _facade_active():
    root = sys.modules.get('sympy')
    return getattr(root, '__symox_version__', None) is not None


def _is_function_class(cls):
    # SymPy's FunctionClass metaclass prints applied heads by bare name.
    application = _classes.get('Application')
    return application is not None and issubclass(cls, application)


def sympy_class_repr(cls):
    """SymPy's spelling of a nominal class, or None when it has no counterpart."""
    entry = CLASSES.get(cls.__name__)
    if entry is None:
        return None
    module = entry[0] if _facade_active() else 'symox'+entry[0][len('sympy'):]
    return "<class '%s.%s'>" % (module, cls.__name__)


class NodeType(type):
    """Metaclass of the nominal node classes."""

    def __repr__(cls):
        if _is_function_class(cls):
            return cls.__name__
        return sympy_class_repr(cls) or type.__repr__(cls)

    # SymPy's public class of a node is its type: ``type(expr) == Eq``,
    # ``type(expr) in (Eq, Ne)`` and lookups of ``type(expr)`` in a dictionary
    # keyed by public classes hold. The public symox class compares equal to
    # the nominal class of its nodes, and both hash by name.
    __hash__ = nominal_hash

    def __eq__(cls, other):
        if cls is other:
            return True
        if isinstance(other, type) and not isinstance(other, NodeType) and other is public_class(cls) \
                and hash(other) == hash(cls):
            return True
        return NotImplemented

    def __ne__(cls, other):
        equal = NodeType.__eq__(cls, other)
        return equal if equal is NotImplemented else not equal

    def __call__(cls, *args, **kwargs):
        # ``type(expr)(*args)`` rebuilds a node in SymPy. Native nodes are
        # created by their public constructors, never by allocating a class.
        constructor = _public_constructor(cls.__name__)
        if constructor is None or constructor is cls:
            raise TypeError("cannot create '%s' instances directly" % cls.__name__)
        return constructor(*args, **kwargs)


def is_native_node(value):
    """Whether ``value`` is a native expression node of symox's own class:
    the base node type or the nominal class it was created with, never a
    user subclass, a host class or a wrapper with its own dispatch."""
    kind = type(value)
    return kind is _core.Expr or type(kind) is NodeType


def native_module(cls):
    """The module recorded in the class dictionary, ignoring the facade."""
    return type.__dict__['__module__'].__get__(cls)


def is_symox_class(cls):
    """Whether ``cls`` is defined by symox itself.

    Only symox's own public classes stand for the nominal class of their
    name: a user class reusing a SymPy name (Torch's ``Max``, Devito's
    ``Derivative``) must not see that node class's attributes, accept its
    node types or report its SymPy path.
    """
    module = native_module(cls) if isinstance(cls, type) else getattr(cls, '__module__', None)
    return isinstance(module, str) and (module == 'symox' or module.startswith('symox.'))


def _class_module(cls):
    module = native_module(cls)
    if isinstance(module, str) and module.startswith('symox.') and _classes.get(cls.__name__) is cls \
            and _facade_active():
        return 'sympy'+module[len('symox'):]
    return module


class _ModuleName(str):
    """A metaclass's own ``__module__`` that also serves its classes'.

    The nominal classes report SymPy's module under the routing facade, as
    their repr does: protocols test ``'sympy' in type(obj).__module__``
    (SageMath's coercion into SR). Reading ``__module__`` on a class first
    consults data descriptors of its metaclass, so this str subclass in the
    metaclass dictionary computes the class value while the metaclass itself
    still reports its own module as a plain string.
    """

    def __new__(cls, module, compute):
        self = str.__new__(cls, module)
        self._compute = compute
        return self

    def __get__(self, cls, owner=None):
        return self if cls is None else self._compute(cls)

    def __set__(self, cls, value):
        type.__dict__['__module__'].__set__(cls, value)


NodeType.__module__ = _ModuleName(NodeType.__module__, _class_module)


def _wrapper_module(cls):
    module = native_module(cls)
    entry = CLASSES.get(cls.__name__)
    if entry is not None and isinstance(module, str) and (module == 'symox' or module.startswith('symox.')) \
            and _facade_active():
        return entry[0]
    return module


class WrapperType(type):
    """Metaclass of public symox classes wrapping native values (sets).

    Under the routing facade they report SymPy's class path, like NodeType.
    """

    def __repr__(cls):
        module = cls.__module__
        if isinstance(module, str) and module.startswith('sympy') and cls.__name__ in CLASSES:
            return "<class '%s.%s'>" % (module, cls.__name__)
        return type.__repr__(cls)


WrapperType.__module__ = _ModuleName(WrapperType.__module__, _wrapper_module)


def _public_constructor(name):
    entry = CLASSES.get(name)
    candidates = []
    if entry is not None:
        candidates.append('symox'+entry[0][len('sympy'):])
    candidates.append('symox')
    for module_name in candidates:
        module = sys.modules.get(module_name)
        if module is None:
            try:
                module = importlib.import_module(module_name)
            except ImportError:
                continue
        value = module.__dict__.get(name)
        if value is not None and not isinstance(value, NodeType):
            return value
    return None


def public_class(cls):
    """The public symox class of a nominal class, from loaded modules only."""
    name = type.__dict__['__name__'].__get__(cls)
    entry = CLASSES.get(name)
    if entry is None:
        return None
    for module_name in ('symox'+entry[0][len('sympy'):], 'symox'):
        module = sys.modules.get(module_name)
        value = getattr(module, '__dict__', {}).get(name)
        if isinstance(value, type) and not isinstance(value, NodeType):
            return value
    return None


def _restore_class(name):
    cls = node_class(name)
    if cls is None:
        raise TypeError('unknown native node class %r' % name)
    return cls


copyreg.pickle(NodeType, lambda cls: (_restore_class, (cls.__name__,)))


def node_class(name):
    """Return the nominal class of a SymPy class name, or None."""
    cls = _classes.get(name)
    if cls is not None or name not in CLASSES:
        return cls
    module, base_names = CLASSES[name]
    # HostTree is the root: a tree containing user subclass instances keeps
    # them alive in its slot while reporting its SymPy class.
    bases = tuple(node_class(base) for base in base_names) or (HostTree,)
    namespace = {'__module__': 'symox'+module[len('sympy'):], '__qualname__': name,
                 '__slots__': ()}
    if name == 'Basic':
        namespace['_sage_'] = _sage_placeholder
    cls = NodeType(name, bases, namespace)
    _classes[name] = cls
    return cls


def _sage_placeholder(self):
    # SymPy's Basic._sage_: SageMath installs the real conversions lazily.
    old_method = self._sage_
    from sage.interfaces.sympy import sympy_init
    sympy_init()
    if old_method == self._sage_:
        raise NotImplementedError('conversion to SageMath is not implemented')
    else:
        return self._sage_()


_sage_placeholder.__name__ = _sage_placeholder.__qualname__ = '_sage_'


def _patchable(name):
    # ABCMeta stores its registry on the class it initialises.
    return not (name.startswith('__') and name.endswith('__')) and not name.startswith('_symox') \
        and not name.startswith('_abc_') and name != '_matches'


def class_attribute(name, attribute):
    """Look ``attribute`` up as SymPy would on the class ``name`` itself.

    Only attributes stored on the nominal classes (SymPy placeholders and
    forwarded patches) are visible. Descriptors bind to the class.
    """
    cls = node_class(name)
    if cls is None or not _patchable(attribute):
        raise AttributeError(attribute)
    for base in cls.__mro__:
        if base is _core.Expr:
            break
        if attribute in base.__dict__:
            value = base.__dict__[attribute]
            getter = getattr(type(value), '__get__', None)
            return value if getter is None else getter(value, None, cls)
    raise AttributeError(attribute)


_anchors = weakref.WeakKeyDictionary()


def nominal_anchor(cls):
    """The nominal class standing for ``cls`` in SymPy's MRO, or None.

    A symox class stands for the nominal class of its own name. A class
    defined outside symox (a user subclass of ``Symbol``, ``Function`` or
    ``Expr``) inherits, as in SymPy, the attributes of its nearest symox
    ancestor that has a SymPy counterpart: patches of ``Basic`` or
    ``Symbol`` reach it, those of an unrelated class of the same name do not.
    """
    try:
        return _anchors[cls]
    except (KeyError, TypeError):
        pass
    anchor = None
    if is_symox_class(cls):
        anchor = node_class(cls.__name__)
    else:
        for base in cls.__mro__:
            if base is object or isinstance(base, NodeType) or not is_symox_class(base):
                continue
            name = base.__dict__.get('_symox_nominal', base.__name__)
            if name in CLASSES:
                anchor = node_class(name)
                break
    try:
        _anchors[cls] = anchor
    except TypeError:
        pass
    return anchor


def _anchor_attribute(cls, attribute, instance):
    anchor = nominal_anchor(cls) if _patchable(attribute) else None
    if anchor is None:
        raise AttributeError(attribute)
    for base in anchor.__mro__:
        if base is HostTree:
            break
        if attribute in base.__dict__:
            value = base.__dict__[attribute]
            getter = getattr(type(value), '__get__', None)
            return value if getter is None else getter(value, instance, cls)
    raise AttributeError(attribute)


def anchored_class_attribute(cls, attribute):
    """Class-level lookup of ``attribute`` on the nominal anchor of ``cls``."""
    return _anchor_attribute(cls, attribute, None)


def anchored_instance_attribute(obj, attribute):
    """Instance-level lookup of ``attribute`` for instances of classes
    defined outside symox, whose MRO does not contain nominal classes."""
    try:
        return _anchor_attribute(type(obj), attribute, obj)
    except AttributeError:
        raise AttributeError("'%s' object has no attribute '%s'" % (type(obj).__name__, attribute)) from None


def has_class_attribute(name, attribute):
    try:
        class_attribute(name, attribute)
    except AttributeError:
        return False
    return True


def forward_setattr(owner, name, attribute, value):
    """Forward a monkeypatch of a public symox class or head to its node class."""
    target = _forward_target(owner, name, attribute)
    if target is not None:
        type.__setattr__(target, attribute, value)


def forward_delattr(owner, name, attribute):
    target = _forward_target(owner, name, attribute)
    if target is not None and attribute in target.__dict__:
        type.__delattr__(target, attribute)


def _forward_target(owner, name, attribute):
    if not _patchable(attribute) or name not in CLASSES:
        return None
    if not is_symox_class(owner):
        return None
    # Native attributes stay authoritative: symox installs its own Python
    # methods on several facades and those must not shadow Rust methods.
    if hasattr(_core.Expr, attribute):
        return None
    return node_class(name)


def _register_numbers():
    numbers.Number.register(node_class('Number'))
    numbers.Real.register(node_class('Float'))
    numbers.Rational.register(node_class('Rational'))
    numbers.Integral.register(node_class('Integer'))


_register_numbers()


try:
    _core._register_node_types(node_class)
except ValueError:
    pass  # The native module outlives a reload of the facade.
