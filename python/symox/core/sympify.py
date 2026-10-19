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
