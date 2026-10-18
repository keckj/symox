"""Native symbolic expression output for maple."""
from ._languages import MapleCodePrinter

def maple_code(expr, assign_to=None, **settings):
    return MapleCodePrinter(settings).doprint(expr, assign_to)
