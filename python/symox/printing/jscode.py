"""Native symbolic expression output for jscode."""
from ._languages import JavascriptCodePrinter

def jscode(expr, assign_to=None, **settings):
    return JavascriptCodePrinter(settings).doprint(expr, assign_to)

# Compatible imports reuse the existing native objects.
from symox._api_exports import export as _export
_export(globals(), 'symox.printing', 'print_jscode')


known_functions = {'Abs': 'Math.abs',
 'Max': 'Math.max',
 'Min': 'Math.min',
 'acos': 'Math.acos',
 'acosh': 'Math.acosh',
 'asin': 'Math.asin',
 'asinh': 'Math.asinh',
 'atan': 'Math.atan',
 'atan2': 'Math.atan2',
 'atanh': 'Math.atanh',
 'ceiling': 'Math.ceil',
 'cos': 'Math.cos',
 'cosh': 'Math.cosh',
 'exp': 'Math.exp',
 'floor': 'Math.floor',
 'log': 'Math.log',
 'sign': 'Math.sign',
 'sin': 'Math.sin',
 'sinh': 'Math.sinh',
 'tan': 'Math.tan',
 'tanh': 'Math.tanh'}
