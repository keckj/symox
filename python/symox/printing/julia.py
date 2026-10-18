"""Native symbolic expression output for julia."""
from ._languages import JuliaCodePrinter

def julia_code(expr, assign_to=None, **settings):
    return JuliaCodePrinter(settings).doprint(expr, assign_to)
