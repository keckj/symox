"""Debug pretty printing using the existing native layout engine."""

from ..printer import Printer


class PrettyPrinter(Printer):
    def doprint(self, expr):
        from symox import _core
        return _core.pretty(expr, **self._settings)


def pretty(expr, **settings):
    return PrettyPrinter(settings).doprint(expr)


def pretty_print(expr, **settings):
    print(pretty(expr, **settings))


pprint = pretty_print
