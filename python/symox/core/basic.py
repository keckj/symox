"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox', 'Atom Basic')


def as_Basic(expr):
    from symox import sympify, SympifyError
    try:
        return sympify(expr, strict=True)
    except SympifyError:
        raise TypeError('Argument must be a Basic object, not `%s`' % type(expr).__name__) from None
__all__.append('as_Basic')
