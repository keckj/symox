"""Native symbolic expression output for jscode."""
from ._languages import JavascriptCodePrinter

def jscode(expr, assign_to=None, **settings):
    return JavascriptCodePrinter(settings).doprint(expr, assign_to)
