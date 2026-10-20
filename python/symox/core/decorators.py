"""Argument-sympifying and operator-priority decorators.

These follow sympy/core/decorators.py (SymPy 1.14). They are used by Python
classes that define binary special methods next to native expressions.
"""
from functools import wraps

from .sympify import SympifyError, sympify


def _sympifyit(arg, retval=None):
    """Decorator to smartly sympify the second argument of a binary method.

    If ``retval`` is given, a failed strict sympification returns it
    (typically ``NotImplemented``); otherwise the SympifyError propagates.
    """
    def deco(func):
        return __sympifyit(func, arg, retval)
    return deco


def __sympifyit(func, arg, retval=None):
    """Decorator to sympify ``arg`` argument for function ``func``."""
    # we support f(a,b) only
    if not func.__code__.co_argcount:
        raise LookupError("func not found")
    # only b is sympified
    assert func.__code__.co_varnames[1] == arg
    if retval is None:
        @wraps(func)
        def __sympifyit_wrapper(a, b):
            return func(a, sympify(b, strict=True))
    else:
        @wraps(func)
        def __sympifyit_wrapper(a, b):
            try:
                # An external class with _op_priority knows how to deal with
                # symbolic objects. Otherwise it must be converted.
                if not hasattr(b, '_op_priority'):
                    b = sympify(b, strict=True)
                return func(a, b)
            except SympifyError:
                return retval
    return __sympifyit_wrapper


def call_highest_priority(method_name):
    """Dispatch a binary special method to the operand with higher _op_priority.

    ``method_name`` is the name of the other operand's method to call, e.g.
    ``@call_highest_priority('__rmul__')`` on ``__mul__``.
    """
    def priority_decorator(func):
        @wraps(func)
        def binary_op_wrapper(self, other):
            if hasattr(other, '_op_priority'):
                if other._op_priority > self._op_priority:
                    f = getattr(other, method_name, None)
                    if f is not None:
                        return f(self)
            return func(self, other)
        return binary_op_wrapper
    return priority_decorator


def sympify_method_args(cls):
    """Class decorator that installs methods marked with ``sympify_return``."""
    for attrname, obj in cls.__dict__.items():
        if isinstance(obj, _SympifyWrapper):
            setattr(cls, attrname, obj.make_wrapped(cls))
    return cls


def sympify_return(*args):
    """Method decorator to sympify its argument and check its type.

    Used with ``sympify_method_args``; ``args`` are ``(parameter, expected
    type)`` and the return value used when sympification or the check fails.
    """
    def wrapper(func):
        return _SympifyWrapper(func, args)
    return wrapper


class _SympifyWrapper:
    """Internal class used by sympify_return and sympify_method_args"""

    def __init__(self, func, args):
        self.func = func
        self.args = args

    def make_wrapped(self, cls):
        func = self.func
        parameters, retval = self.args

        [(parameter, expectedcls)] = parameters

        # Handle forward references to the current class using strings
        if expectedcls == cls.__name__:
            expectedcls = cls

        # A failure at import time is not recoverable.
        nargs = func.__code__.co_argcount
        if nargs != 2:
            raise RuntimeError('sympify_return can only be used with 2 argument functions')
        if func.__code__.co_varnames[1] != parameter:
            raise RuntimeError('parameter name mismatch "%s" in %s' %
                    (parameter, func.__name__))

        @wraps(func)
        def _func(self, other):
            if not hasattr(other, '_op_priority'):
                try:
                    other = sympify(other, strict=True)
                except SympifyError:
                    return retval
            if not isinstance(other, expectedcls):
                return retval
            return func(self, other)

        return _func
