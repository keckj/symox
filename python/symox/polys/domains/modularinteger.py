"""Module of the reference modular integer element factory.

The reference implementation (sympy/polys/domains/modularinteger.py) builds
one ``ModularInteger`` subclass per (modulus, domain, symmetric, parent)
and memoizes those classes in ``_modular_integer_cache``. Native ``GF(p)``
elements are instances of the single native ``DomainElement`` type and need
no generated classes, so the cache exists and is always empty: clearing it
(as downstream cache managers do) is exact and has no other effect.
"""
from ..._core import CoercionFailed, DomainElement

_modular_integer_cache: dict = {}

__all__ = []
