"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox', 'Atom Basic')


def as_Basic(expr):
    from symox import sympify,Basic
    result=sympify(expr)
    if not isinstance(result,Basic):raise TypeError('expected a symbolic Basic object')
    return result
__all__.append('as_Basic')
