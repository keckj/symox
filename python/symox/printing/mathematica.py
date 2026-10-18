"""Native symbolic expression output for mathematica."""
from ._languages import MCodePrinter

def mathematica_code(expr, assign_to=None, **settings):
    return MCodePrinter(settings).doprint(expr, assign_to)
