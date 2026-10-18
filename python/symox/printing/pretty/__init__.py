"""Pretty printer compatibility exports."""

from .pretty import PrettyPrinter, pretty, pretty_print, pprint

__all__ = ["PrettyPrinter", "pretty", "pretty_print", "pprint"]

# Compatible imports reuse the existing native objects.
from symox._api_exports import export as _export
_export(globals(), 'symox.printing', 'pager_print pprint_try_use_unicode pprint_use_unicode')
