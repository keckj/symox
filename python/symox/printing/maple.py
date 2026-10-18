"""Native symbolic expression output for maple."""
from ._languages import MapleCodePrinter

def maple_code(expr, assign_to=None, **settings):
    return MapleCodePrinter(settings).doprint(expr, assign_to)

# Compatible imports reuse the existing native objects.
from symox._api_exports import export as _export
_export(globals(), 'symox.printing', 'print_maple_code')
