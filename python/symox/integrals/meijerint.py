"""Integration using Meijer G-functions: the native port of the reference
``meijerint_indefinite``, ``meijerint_definite`` and ``meijerint_inversion``.

Each function returns ``None`` when the method does not apply, like the
reference.
"""
from symox import _core
from symox import sympify


def _first(result):
    return None if result is None else result[0]


def meijerint_indefinite(f, x):
    """Antiderivative of ``f`` by Meijer G-functions, or None."""
    return _first(_core._meijerint_debug('indefinite', sympify(f), sympify(x)))


def meijerint_definite(f, x, a, b):
    """``(value, condition)`` for the integral of ``f`` over ``(a, b)``, or None."""
    return _core._meijerint_debug('definite', sympify(f), sympify(x), sympify(a), sympify(b))


def meijerint_inversion(f, x, t):
    """Inverse Laplace transform of ``f(x)`` at ``t`` by Meijer G-functions, or None."""
    return _first(_core._meijerint_debug('inversion', sympify(f), sympify(x), sympify(t)))


__all__ = ['meijerint_indefinite', 'meijerint_definite', 'meijerint_inversion']
