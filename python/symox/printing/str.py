"""Subclassable textual expression printer."""

from .precedence import precedence
from .printer import Printer


import threading

# The ``abbrev`` setting of the printer currently delegating to ``str``;
# None outside a printer, where StrPrinter's defaults apply (as for SymPy's
# ``Basic.__str__``, which prints through ``sstr``).
_active = threading.local()


def abbreviate_units():
    """Whether quantities print by their abbreviation in the current context."""
    abbrev = getattr(_active, 'abbrev', None)
    return bool(StrPrinter._default_settings.get('abbrev') if abbrev is None else abbrev)


class StrPrinter(Printer):
    printmethod = "_sympystr"
    _default_settings = {
        "order": None, "full_prec": "auto", "sympy_integers": False, "abbrev": False,
        "perm_cyclic": True, "min": None, "max": None, "dps": None,
    }

    def emptyPrinter(self, expr):
        previous = getattr(_active, 'abbrev', None)
        _active.abbrev = bool(self._settings.get('abbrev'))
        try:
            return super().emptyPrinter(expr)
        finally:
            _active.abbrev = previous

    def parenthesize(self, item, level, strict=False):
        value = self._print(item)
        p = precedence(item)
        return f"({value})" if p < level or (not strict and p == level) else value

    def _print_UnevaluatedExpr(self, expr):
        return self._print(expr.args[0])

    def _print_Symbol(self, expr):
        return expr.name

    def _print_Dummy(self, expr):
        return '_' + expr.name

    def _print_Wild(self, expr):
        return expr.name + '_'

    _print_MatrixSymbol = _print_Symbol


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
