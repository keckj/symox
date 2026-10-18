"""Native symbolic expression output for glsl."""
from ._languages import GLSLPrinter

def glsl_code(expr, assign_to=None, **settings):
    return GLSLPrinter(settings).doprint(expr, assign_to)
