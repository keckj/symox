"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox.utilities.cache', 'cacheit clear_cache')
_export(globals(), 'functools', 'cached_property')
