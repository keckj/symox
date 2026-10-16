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
