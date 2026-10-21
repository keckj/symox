"""Integration using Meijer G-functions: the native port of the reference
``meijerint_indefinite``, ``meijerint_definite`` and ``meijerint_inversion``.

Each function returns ``None`` when the method does not apply, like the
reference; ``meijerint_indefinite`` also raises the reference errors.
"""
from symox import _core
from symox import sympify


def _first(result):
    return None if result is None else result[0]


def meijerint_indefinite(f, x):
    """Antiderivative of ``f`` by Meijer G-functions, or None.

    The errors the reference raises propagate with its classes (as
    ``NotImplementedError('Gammas partially over the strip.')`` for
    ``sin(x)**3/x**3``); ``integrate`` itself goes on with its other
    methods.
    """
    result = _core._meijerint_indefinite(sympify(f), sympify(x))
    if isinstance(result, tuple):
        class_name, message = result
        raise _error_class(class_name)(message)
    return result


def _error_class(name):
    if name == 'MellinTransformStripError':
        from symox._transforms import MellinTransformStripError
        return MellinTransformStripError
    if name == 'PolynomialError':
        from symox import PolynomialError
        return PolynomialError
    import builtins
    cls = getattr(builtins, name, None)
    if isinstance(cls, type) and issubclass(cls, Exception):
        return cls
    return NotImplementedError


def meijerint_definite(f, x, a, b):
    """``(value, condition)`` for the integral of ``f`` over ``(a, b)``, or None."""
    return _core._meijerint_debug('definite', sympify(f), sympify(x), sympify(a), sympify(b))


def meijerint_inversion(f, x, t):
    """Inverse Laplace transform of ``f(x)`` at ``t`` by Meijer G-functions, or None."""
    return _first(_core._meijerint_debug('inversion', sympify(f), sympify(x), sympify(t)))


__all__ = ['meijerint_indefinite', 'meijerint_definite', 'meijerint_inversion']
