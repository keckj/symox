"""Substitution of function classes, as in SymPy's ``f(x).subs(f, g)``.

SymPy treats a function class ``old`` in ``subs`` through
``Application._eval_subs``: every application ``old(*args)`` whose
argument count is accepted by the callable ``new`` becomes
``new(*args)``, the arguments being substituted first.  A few classes
override that rule (``exp``, ``cot``, ``CRootOf`` and the expressions
with limits) and are reproduced here.  ``xreplace`` never matches a
function class because it only visits expression nodes.
"""
import inspect
from collections.abc import Mapping

from . import _core

# Function classes whose symox spelling accepts ``*args``; the values
# are SymPy's ``nargs``.
_NARGS = {
    'KroneckerDelta': (2, 3), 'Mod': (2,), 'Rem': (2,), 'acoth': (1,),
    'acsc': (1,), 'acsch': (1,), 'asec': (1,), 'asech': (1,), 'csch': (1,),
    'frac': (1,), 'sech': (1,), 'sinc': (1,),
}
# SymPy spells these as plain Python functions rather than classes.
_PLAIN = frozenset({'E1', 'Eijk', 'cbrt', 'sqrt', 'udivisor_sigma'})
_LIMITS = frozenset({'Integral', 'Sum', 'Product'})
# Classes whose own _eval_subs never replaces their head by a class.
_KEEP_HEAD = frozenset({'cot', 'CRootOf', 'ComplexRootOf', 'And', 'Or', 'Xor', 'Piecewise'})
# SymPy's nargs for the native Application classes that are not
# Naturals0.
_CLASS_NARGS = {
    'Not': (1,), 'BesselBase': (2,), 'SphericalBesselBase': (2,),
    'SphericalHankelBase': (2,),
}


def _Application():
    from ._function_heads import Application
    return Application


def function_class(value):
    """Whether ``value`` is a SymPy function class such as ``sin`` or ``Function('f')``."""
    from . import _is_function_head
    from ._function_heads import is_builtin_head
    from ._function_heads import Application
    if is_builtin_head(value):
        return value.__name__ not in _PLAIN
    if _is_function_head(value):
        return True
    # SymPy's FunctionClass is the metaclass of every Application
    # subclass, such as Max, And or Piecewise.
    return isinstance(value, type) and issubclass(value, Application)


def function_class_name(value):
    """SymPy's metaclass name for a function class, else None."""
    from . import _is_function_head
    if not function_class(value):
        return None
    return 'UndefinedFunction' if _is_function_head(value) else 'FunctionClass'


def _function_like(value):
    # SymPy requires ``new.is_Function`` and ``callable(new)``.
    if function_class(value):
        return True
    return (isinstance(value, _core.Expr) and callable(value)
            and getattr(value, 'is_Function', False) is True)


def _accepts(new, count):
    from ._function_heads import is_builtin_head
    if is_builtin_head(new):
        if new.__name__ in _NARGS:
            return count in _NARGS[new.__name__]
        try:
            inspect.signature(new).bind(*range(count))
        except TypeError:
            return False
        return True
    nargs = getattr(new, 'nargs', None)
    if isinstance(new, type) and (nargs is None or isinstance(nargs, property)):
        # A native class such as Max or Not has no class-level nargs;
        # SymPy's FunctionClass derives it from ``eval`` (Naturals0 for
        # ``*args``).
        if not issubclass(new, _Application()):
            return False
        for base in new.__mro__:
            if base.__name__ in _CLASS_NARGS:
                return count in _CLASS_NARGS[base.__name__]
        return True
    if nargs is None:
        return False
    try:
        return count in nargs
    except TypeError:
        # Lambda.nargs is a bare native FiniteSet node.
        return any(item == count for item in nargs.args)


def _same_head(head, old):
    return head is old or (function_class(head) and head == old)


def substitute(expr, old, new):
    """``expr.subs(old, new)`` for a function class ``old``."""
    from . import Pow, exp, preorder_traversal
    if not isinstance(expr, _core.Expr):
        return expr
    if any(getattr(node, 'class_name', None) in _LIMITS for node in preorder_traversal(expr)):
        # ExprWithLimits._eval_subs intersects the free symbols of the
        # limits with ``old.free_symbols``, a property of the class.
        raise TypeError('function-class substitution into Integral, Sum or Product is not supported')
    bases = getattr(old, '__mro__', None) or (old,)
    if function_class(old) and any(getattr(base, '__name__', None) in _KEEP_HEAD for base in bases):
        # These classes override _eval_subs and keep their head.
        return expr
    if _function_like(new):
        def query(node):
            return (_same_head(getattr(node, 'func', None), old)
                    and _accepts(new, len(node.args)))

        def value(node):
            return new(*node.args)
    elif old is exp:
        # exp._eval_subs: exp(a).subs(exp, y) is y**a.
        def query(node):
            return getattr(node, 'func', None) is exp

        def value(node):
            return Pow(new, node.args[0])
    else:
        return expr
    return expr.replace(query, value, simultaneous=False)


def _sequence(arg1, arg2):
    from . import Dict
    if arg2 is not None:
        return [(arg1, arg2)], False
    if isinstance(arg1, set):
        return list(arg1), True
    if isinstance(arg1, (Dict, Mapping)):
        return list(arg1.items()), True
    return list(arg1), False


def _one_shot(value):
    try:
        return iter(value) is value
    except TypeError:
        return False


def prepared(args):
    """``args`` with one-shot pair iterators materialized once."""
    from . import Dict
    if len(args) != 1 or isinstance(args[0], (Mapping, set, Dict, str, bytes)):
        return args
    try:
        items = list(args[0])
    except TypeError:
        return args
    inner = [tuple(item) if _one_shot(item) else item for item in items]
    if _one_shot(args[0]) or any(a is not b for a, b in zip(inner, items)):
        return (inner,)
    return args


def _involves_function_class(args):
    if not args or len(args) > 2:
        return False
    try:
        items, _ = _sequence(args[0], args[1] if len(args) == 2 else None)
        return any(function_class(old) or function_class(new) for old, new in items)
    except (TypeError, ValueError, IndexError):
        return False


def _sympify_pair(old, new):
    from . import Symbol, sympify
    if not function_class(old):
        old = Symbol(old) if isinstance(old, str) else sympify(old)
    if not _function_like(new) and not function_class(new):
        new = sympify(new)
    return old, new


def subs(original, self, args, kwargs):
    """SymPy's ``Basic.subs`` for a sequence containing function classes."""
    from . import nan, oo, zoo
    from .core.sorting import _nodes, default_sort_key, ordered
    items, unordered = _sequence(args[0], args[1] if len(args) == 2 else None)
    sequence = [_sympify_pair(old, new) for old, new in items]
    # Identical pairs are skipped; the native substitution skips the
    # remaining equal ones itself.
    sequence = [(old, new) for old, new in sequence if old is not new]
    options = dict(kwargs)
    simultaneous = options.pop('simultaneous', False)
    if unordered:
        table = dict(sequence)
        keys = list(ordered(table, default=False, keys=(lambda item: -_nodes(item), default_sort_key)))
        sequence = [(key, table[key]) for key in keys]
        if not simultaneous:
            illegal = (nan, zoo, oo, -oo)
            redo = [index for index, (_, new) in enumerate(sequence)
                    if not function_class(new) and new in illegal]
            for index in reversed(redo):
                sequence.insert(0, sequence.pop(index))
    if simultaneous:
        if any(function_class(new) for _, new in sequence):
            # SymPy builds Dummy(commutative=new.is_commutative), a
            # property of the class.
            raise ValueError('Dummy commutativity must be True or False.')
        # A function class replaced by a dummy product matches nothing.
        rest = [(old, new) for old, new in sequence if not function_class(old)]
        return original(self, rest, simultaneous=True, **options) if rest else self
    result = self
    pending = []
    for old, new in sequence + [(None, None)]:
        if old is not None and not function_class(old):
            # An absent expression replaced by a function class is a
            # no-op; a present one is left to the native error.
            if not function_class(new) or result.has(old):
                pending.append((old, new))
            continue
        if pending:
            result = original(result, pending, **options)
            pending = []
            if not isinstance(result, _core.Expr):
                return result
        if old is not None:
            result = substitute(result, old, new)
    return result


def xreplace_rule(rule):
    """``rule`` without function-class keys, which xreplace never matches."""
    if isinstance(rule, Mapping) and any(function_class(key) for key in rule):
        return {key: value for key, value in rule.items() if not function_class(key)}
    return rule


def install():
    from functools import wraps
    from . import _core as core
    expr_subs, expr_xreplace = core.Expr.subs, core.Expr.xreplace

    @wraps(expr_subs)
    def subs_method(self, *args, **kwargs):
        args = prepared(args)
        if _involves_function_class(args):
            return subs(expr_subs, self, args, kwargs)
        return expr_subs(self, *args, **kwargs)

    @wraps(expr_xreplace)
    def xreplace_method(self, rule, *args, **kwargs):
        return expr_xreplace(self, xreplace_rule(rule), *args, **kwargs)

    core.Expr.subs, core.Expr.xreplace = subs_method, xreplace_method

    matrix_subs = core.Matrix.subs

    @wraps(matrix_subs)
    def matrix_subs_method(self, *args, **kwargs):
        args = prepared(args)
        if _involves_function_class(args):
            # SymPy's MatrixBase.subs maps Expr.subs over the entries.
            return self.applyfunc(lambda entry: entry.subs(*args, **kwargs))
        return matrix_subs(self, *args, **kwargs)

    core.Matrix.subs = matrix_subs_method
