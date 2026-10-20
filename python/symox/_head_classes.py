"""Rebinding of the native gate and lattice constructors to their classes.

SymPy's logic gates (``logic.boolalg``) and ``Max``/``Min``
(``functions.elementary.miscellaneous``) are classes. The root namespace
binds the native constructors while it initializes, before the modules
defining the classes can be imported, and a few modules capture them then.
Once initialization is done every symox module sees the classes.
"""
import sys
import types

_GATES = ('And', 'Or', 'Not', 'Xor', 'Nand', 'Nor', 'Xnor', 'Implies', 'Equivalent', 'ITE')
_LATTICE_FUNCTIONS = ('Max', 'Min')


# Function classes of SymPy kept by their own constructors: Piecewise
# normalizes its conditions before the kernel sees them.
_OWN_CONSTRUCTORS = frozenset(('Piecewise',))
_HEAD_OPTIONS = frozenset(('_abstract', '_heads', '_arity'))


def _function_classes(namespace, replacements):
    # The remaining heads SymPy defines as function classes are built in
    # Python (plain constructors, or classes outside the function hierarchy):
    # bind head classes over them, and the public abstract classes over
    # their abstract bases, so subclassing and issubclass follow SymPy.
    from ._function_heads import _in_function_hierarchy, abstract_class, builtin_head, is_builtin_head
    from ._sympy_classes import CLASSES
    for name, old in list(namespace.items()):
        if (name not in CLASSES or name in _OWN_CONSTRUCTORS or getattr(old, '__name__', None) != name
                or not _in_function_hierarchy(name) or is_builtin_head(old)):
            continue
        if isinstance(old, type):
            if old.__module__ != 'symox' and not old.__module__.startswith('symox.'):
                continue
            if old.__dict__.get('_abstract', False):
                new = abstract_class(name)
                for key, value in old.__dict__.items():
                    if not key.startswith('__') and key not in _HEAD_OPTIONS:
                        type.__setattr__(new, key, value)
            elif old.__bases__ != (object,) and not old.__dict__.get('_heads'):
                continue
            else:
                new = builtin_head(name, old)
        elif isinstance(old, types.FunctionType):
            new = builtin_head(name, old)
        else:
            continue
        replacements[id(old)] = (old, new)
        namespace[name] = new


def install(namespace):
    from .logic import boolalg
    from .functions.elementary import miscellaneous
    replacements = {}
    _function_classes(namespace, replacements)
    for module, names in ((boolalg, _GATES), (miscellaneous, _LATTICE_FUNCTIONS)):
        for name in names:
            cls = getattr(module, name)
            old = namespace.get(name)
            if old is not None and old is not cls:
                replacements[id(old)] = (old, cls)
            namespace[name] = cls
    root = namespace['__name__']
    for module_name, module in list(sys.modules.items()):
        if module is None or module_name == root + '._core' or not (
                module_name == root or module_name.startswith(root + '.')):
            continue
        if type(module).__name__ == '_LazyModule':
            # Not executed yet (importlib.util.LazyLoader): it captures the
            # classes when it runs. Touching it would execute it.
            continue
        variables = getattr(module, '__dict__', None)
        if variables is None:
            continue
        for key, value in list(variables.items()):
            entry = replacements.get(id(value))
            if entry is not None and entry[0] is value:
                variables[key] = entry[1]
