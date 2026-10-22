"""Algebraic number fields, as sympy.polys.numberfields.

Minimal polynomials, primitive elements and field isomorphisms are native.
Orders, maximal orders (Round Two), prime decomposition and the module
machinery are ports of the reference pure-Python algorithms onto the native
exact domains. The submodules that need the full package are loaded on
first access, since this package is imported while symox initialises.
"""
from ..._core import minimal_polynomial, minpoly, primitive_element, to_number_field
from .galoisgroups import galois_group

__all__ = ['minimal_polynomial', 'minpoly', 'primitive_element', 'to_number_field',
           'galois_group', 'field_isomorphism', 'isolate', 'round_two',
           'prime_decomp', 'prime_valuation', 'galois_group']

from .subfield import field_isomorphism

_LAZY = {
    'isolate': 'utilities',
    'round_two': 'basis',
    'prime_decomp': 'primes',
    'prime_valuation': 'primes',
}


def __getattr__(name):
    module = _LAZY.get(name)
    if module is None:
        raise AttributeError("module %r has no attribute %r" % (__name__, name))
    import importlib
    value = getattr(importlib.import_module('.' + module, __name__), name)
    globals()[name] = value
    return value
