"""Context-local global parameters.

``global_parameters`` follows sympy/core/parameters.py (SymPy 1.14). Its
``evaluate`` attribute reads and sets the same context-local state as the
``evaluate`` context manager, so native constructors and Python code agree.
"""
from contextlib import contextmanager

from .. import evaluate
from .._public_core import _evaluation
__all__ = ['evaluate']


class _global_parameters:
    """Global parameters seen by the native constructors.

    ``evaluate`` is stored per thread and per context. ``distribute`` and
    ``exp_is_pow`` are fixed by the native canonicalization rules (numeric
    coefficients distribute over sums and ``exp(x)`` is a function), so only
    their reference default values can be set.
    """
    _fixed = {'distribute': True, 'exp_is_pow': False}

    @property
    def evaluate(self):
        return _evaluation.get()

    @evaluate.setter
    def evaluate(self, value):
        _evaluation.set(value)

    def __getattr__(self, name):
        try:
            return self._fixed[name]
        except KeyError:
            raise AttributeError(
                "'_global_parameters' object has no attribute %r" % name) from None

    def __setattr__(self, name, value):
        if name in self._fixed:
            if value != self._fixed[name]:
                raise NotImplementedError(
                    'native canonicalization does not support %s=%r' % (name, value))
            return
        super().__setattr__(name, value)

global_parameters = _global_parameters()


@contextmanager
def distribute(x):
    """Control automatic distribution of Number over Add (fixed natively)."""
    old = global_parameters.distribute
    try:
        global_parameters.distribute = x
        yield
    finally:
        global_parameters.distribute = old


@contextmanager
def _exp_is_pow(x):
    """Control whether exp(x) is a Pow (fixed natively to the default)."""
    old = global_parameters.exp_is_pow
    try:
        global_parameters.exp_is_pow = x
        yield
    finally:
        global_parameters.exp_is_pow = old

__all__ += ['global_parameters', 'distribute']
