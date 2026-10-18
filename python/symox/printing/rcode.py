"""Native symbolic expression output for rcode."""
from ._languages import RCodePrinter

def rcode(expr, assign_to=None, **settings):
    return RCodePrinter(settings).doprint(expr, assign_to)
