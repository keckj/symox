"""Subclassable textual expression printer."""

from .precedence import precedence
from .printer import Printer


class StrPrinter(Printer):
    printmethod = "_sympystr"

    def parenthesize(self, item, level, strict=False):
        value = self._print(item)
        p = precedence(item)
        return f"({value})" if p < level or (not strict and p == level) else value

    def _print_UnevaluatedExpr(self, expr):
        return self._print(expr.args[0])


def sstr(expr, **settings):
    return StrPrinter(settings).doprint(expr)

# Compatible imports reuse the existing native objects.
from symox._api_exports import export as _export
_export(globals(), 'symox.printing', 'sstrrepr')


class StrReprPrinter(StrPrinter):
    def emptyPrinter(self,expr):
        if isinstance(expr,str):return repr(expr)
        if isinstance(expr,list):return '['+', '.join(self._print(value) for value in expr)+']'
        if isinstance(expr,tuple):return '('+', '.join(self._print(value) for value in expr)+(',' if len(expr)==1 else '')+')'
        if isinstance(expr,dict):return '{'+', '.join(self._print(key)+': '+self._print(value) for key,value in expr.items())+'}'
        return super().emptyPrinter(expr)

__all__ += ['StrReprPrinter']
