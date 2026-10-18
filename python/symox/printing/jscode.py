"""Native symbolic expression output for jscode."""
from ._languages import JavascriptCodePrinter

def jscode(expr, assign_to=None, **settings):
    return JavascriptCodePrinter(settings).doprint(expr, assign_to)

# Compatible imports reuse the existing native objects.
from symox._api_exports import export as _export
_export(globals(), 'symox.printing', 'print_jscode')
