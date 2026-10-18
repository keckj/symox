"""Native symbolic expression output for rcode."""
from ._languages import RCodePrinter

def rcode(expr, assign_to=None, **settings):
    return RCodePrinter(settings).doprint(expr, assign_to)

# Compatible imports reuse the existing native objects.
from symox._api_exports import export as _export
_export(globals(), 'symox.printing', 'print_rcode')


known_functions = {'Abs': 'abs',
 'Max': 'max',
 'Min': 'min',
 'acos': 'acos',
 'acosh': 'acosh',
 'asin': 'asin',
 'asinh': 'asinh',
 'atan': 'atan',
 'atan2': 'atan2',
 'atanh': 'atanh',
 'beta': 'beta',
 'ceiling': 'ceiling',
 'cos': 'cos',
 'cosh': 'cosh',
 'digamma': 'digamma',
 'erf': 'erf',
 'exp': 'exp',
 'factorial': 'factorial',
 'floor': 'floor',
 'gamma': 'gamma',
 'log': 'log',
 'sign': 'sign',
 'sin': 'sin',
 'sinh': 'sinh',
 'sqrt': 'sqrt',
 'tan': 'tan',
 'tanh': 'tanh',
 'trigamma': 'trigamma'}
