"""Native exact coefficient domains and their closed arithmetic elements.

Supported coefficient fields use the native ZZ, QQ, GF(p), rational-function,
and certified algebraic-field backends. Domain identity includes the actual
symbol assumptions, generator order, and selected algebraic embedding.
"""
from ..._core import Domain, DomainElement, QQ, ZZ, GF, DomainError, CoercionFailed, ExactQuotientFailed

__all__ = ["Domain", "DomainElement", "QQ", "ZZ", "GF", "DomainError", "CoercionFailed", "ExactQuotientFailed"]
