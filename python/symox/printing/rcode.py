"""Native symbolic expression output for rcode."""
from ._languages import RCodePrinter

def rcode(expr, assign_to=None, **settings):
    return RCodePrinter(settings).doprint(expr, assign_to)

# Compatible imports reuse the existing native objects.
from symox._api_exports import export as _export
_export(globals(), 'symox.printing', 'print_rcode')
