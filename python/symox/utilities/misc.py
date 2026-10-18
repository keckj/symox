from .. import capture,filldedent
__all__ = ['capture','filldedent']

# Compatible imports reuse the existing native objects.
from symox._api_exports import export as _export
_export(globals(), 'symox.ntheory._intfunc', 'as_int')
