from .. import sympify,SympifyError
from .._sympify import converter
__all__ = ['sympify','SympifyError']

class CantSympify:
    """Marker for objects which must not be implicitly converted to scalars."""

def kernS(s):
    from symox import parse_expr
    return parse_expr(s,evaluate=False) if isinstance(s,str) else sympify(s)
__all__ += ['CantSympify','kernS']
__all__ += ['converter']


def _sympify(a):
    """Short version of sympify for internal usage: ``sympify(a, strict=True)``.

    Only objects with a registered or implicit symbolic conversion are
    accepted; strings and other non-Basic objects raise SympifyError.
    """
    return sympify(a, strict=True)
__all__ += ['_sympify']
