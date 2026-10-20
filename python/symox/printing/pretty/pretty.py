"""Debug pretty printing using the existing native layout engine."""

from .._legacy import LegacyPrinter as Printer


_use_unicode=True


class PrettyPrinter(Printer):
    def doprint(self, expr):
        from symox import _core
        settings=dict(self._settings)
        if settings.get("use_unicode") is None:settings["use_unicode"]=_use_unicode
        return _core.pretty(expr, **settings)


def pretty(expr, **settings):
    return PrettyPrinter(settings).doprint(expr)


def pretty_print(expr, **settings):
    print(pretty(expr, **settings))


pprint = pretty_print


def pprint_use_unicode(flag=None):
    global _use_unicode
    previous=_use_unicode
    if flag is not None:_use_unicode=bool(flag)
    return previous


def pprint_try_use_unicode():
    import sys
    encoding=getattr(sys.stdout,'encoding',None)
    if encoding:
        try:'α√∫'.encode(encoding)
        except (UnicodeError,LookupError):return
        pprint_use_unicode(True)

# Compatible imports reuse the existing native objects.
from symox._api_exports import export as _export
_export(globals(), 'symox.printing', 'pager_print')

from .pretty_symbology import (U,annotated,center_pad,hobj,is_subscriptable_in_unicode,
 pretty_atom,pretty_symbol,pretty_try_use_unicode,pretty_use_unicode,vobj,xobj,xsym)
