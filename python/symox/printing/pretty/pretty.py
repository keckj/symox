"""Debug pretty printing using the existing native layout engine."""

from .._legacy import LegacyPrinter as Printer


_use_unicode=True


class PrettyPrinter(Printer):
    def doprint(self, expr):
        global _use_unicode
        from symox import _core
        if getattr(type(expr), '_symox_pretty_as_str', False):
            # The reference PrettyPrinter has no handler for the class
            # (Quaternion): its emptyPrinter renders the str form.
            return str(expr)
        settings=dict(self._settings)
        if settings.get("use_unicode") is None:settings["use_unicode"]=_use_unicode
        try:
            return _core.pretty(expr, **settings)
        except TypeError as exc:
            if 'cannot convert object' not in str(exc):raise
        # A Python-level object the native layout engine cannot represent:
        # the reference dispatches to its _pretty method, otherwise its
        # emptyPrinter renders the str form.
        if not hasattr(expr, '_pretty'):return str(expr)
        previous=_use_unicode
        _use_unicode=bool(settings["use_unicode"])
        try:
            self._settings=settings
            form=expr._pretty(self)
        finally:
            _use_unicode=previous
        if hasattr(form, 'render'):
            return form.render(wrap_line=settings.get('wrap_line', True),
                               num_columns=settings.get('num_columns'))
        return str(form)

    def _print(self, expr, **kwargs):
        from .stringpict import prettyForm
        if hasattr(expr, 'picture'):return expr
        if isinstance(expr, str):return prettyForm(expr)
        if hasattr(expr, '_pretty'):
            try:return expr._pretty(self)
            except TypeError:pass
        text=PrettyPrinter(self._settings).doprint(expr)
        return prettyForm(text)


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
