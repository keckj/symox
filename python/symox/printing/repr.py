"""Native constructor representation and recursive container formatting."""
from .printer import Printer

class ReprPrinter(Printer):
    printmethod='_sympyrepr'
    def emptyPrinter(self,expr):
        from symox import _core
        if isinstance(expr,(list,tuple,set,frozenset)):
            items=', '.join(self._print(item) for item in expr)
            if isinstance(expr,list):return '['+items+']'
            if isinstance(expr,tuple):return '('+items+(',' if len(expr)==1 else '')+')'
            if isinstance(expr,frozenset):return 'frozenset({'+items+'})'
            return '{'+items+'}' if expr else 'set()'
        if isinstance(expr,dict):return '{'+', '.join(self._print(k)+': '+self._print(v) for k,v in expr.items())+'}'
        if isinstance(expr,(str,bytes)):return repr(expr)
        if hasattr(expr,'srepr'):return expr.srepr()
        try:return _core.srepr(expr)
        except (TypeError,ValueError):return repr(expr)

def srepr(expr,**settings):return ReprPrinter(settings).doprint(expr)
