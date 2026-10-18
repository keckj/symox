"""Native symbolic expression output for octave."""
from ._languages import OctaveCodePrinter

def octave_code(expr, assign_to=None, **settings):
    return OctaveCodePrinter(settings).doprint(expr, assign_to)
