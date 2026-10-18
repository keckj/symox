"""Native fortran code generation entry point."""
from symox._core import fcode
__all__ = ['fcode']

# Compatible imports reuse the existing native objects.
from symox._api_exports import export as _export
_export(globals(), 'symox.printing', 'print_fcode')
