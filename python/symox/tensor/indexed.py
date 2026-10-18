"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox', 'Idx Indexed IndexedBase')

class IndexException(Exception):
    """An indexed expression has inconsistent or insufficient index information."""
__all__.append('IndexException')
