"""Native symbolic expression output for mathematica."""
from ._languages import MCodePrinter

def mathematica_code(expr, assign_to=None, **settings):
    return MCodePrinter(settings).doprint(expr, assign_to)

from ._languages import _COMMON,_EXTRA,_MMA
known_functions={name:[(lambda *args:True,_MMA.get(name,name[0].upper()+name[1:]))] for name in _COMMON|_EXTRA['mathematica']}
