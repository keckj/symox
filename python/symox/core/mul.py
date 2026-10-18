"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox', 'Mul prod')

class NC_Marker:
    """Private sentinel used to mark a noncommutative factor group."""

def expand_2arg(expr):
    from symox import Mul,Add
    if not getattr(expr,'is_Mul',False) or len(expr.args)!=2:return expr
    a,b=expr.args
    if a.is_Add:return Add(*(term*b for term in a.args))
    if b.is_Add:return Add(*(a*term for term in b.args))
    return expr
__all__ += ['NC_Marker','expand_2arg']
