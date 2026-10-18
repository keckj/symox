"""Native scalar arrays and concrete tensor operations."""
from .array import *
from .array import __all__

# Compatible imports reuse the existing native objects.
from symox._api_exports import export as _export
_export(globals(), 'symox', 'Idx Indexed IndexedBase get_contraction_structure get_indices shape')
