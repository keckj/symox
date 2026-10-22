"""Native exact coefficient domains and their closed arithmetic elements.

Supported coefficient fields use the native ZZ, QQ, GF(p), rational-function,
and certified algebraic-field backends. Domain identity includes the actual
symbol assumptions, generator order, and selected algebraic embedding.
"""
from ..._core import Domain, DomainElement, QQ, ZZ, GF, DomainError, CoercionFailed, ExactQuotientFailed
from ..._modular_ring import GF

__all__ = ["Domain", "DomainElement", "QQ", "ZZ", "GF", "DomainError", "CoercionFailed", "ExactQuotientFailed"]

from ... import _core
if hasattr(_core, 'EX'):
    from ..._core import EX, EXRAW, RR, CC, ZZ_I, QQ_I
    from ._types import *
    from ._types import __all__ as _type_names
    FF = GF
    __all__ += ['EX', 'EXRAW', 'RR', 'CC', 'ZZ_I', 'QQ_I', 'FF'] + _type_names

from ._compat import *
from ._compat import __all__ as _compat_names
__all__ += _compat_names
