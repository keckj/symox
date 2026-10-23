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
        if _composed(self, expr):
            return self._render_composed(expr, settings)
        try:
            return _core.pretty(expr, **settings)
        except TypeError as exc:
            if not _unconvertible(exc):raise
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

    # ------------------------------------------------------------------
    # Composition of pictures (reference PrettyPrinter methods for
    # containers and Python-level objects). Native expressions keep the
    # native layout and contribute their picture with its baseline.
    # ------------------------------------------------------------------

    @property
    def _use_unicode(self):
        return bool(self._settings.get('use_unicode', _use_unicode))

    def _render_composed(self, expr, settings):
        global _use_unicode
        from . import _stringpict_ref as sp
        previous=_use_unicode
        _use_unicode=bool(settings["use_unicode"])
        old_wrap=sp._GLOBAL_WRAP_LINE
        try:
            self._settings=settings
            form=self._compose(expr)
            return form.render(wrap_line=settings.get('wrap_line', True),
                               num_columns=settings.get('num_columns'),
                               use_unicode=_use_unicode)
        finally:
            _use_unicode=previous
            sp._GLOBAL_WRAP_LINE=old_wrap

    def _compose(self, expr):
        from symox import _core
        from ._stringpict_ref import stringPict, prettyForm
        from .stringpict import stringPict as _legacyPict
        if isinstance(expr, stringPict):return expr
        if isinstance(expr, _legacyPict):return _from_legacy(expr)
        if isinstance(expr, str):return prettyForm(expr)
        # The reference Printer tries the object's print method first.
        if hasattr(expr, '_pretty') and not isinstance(expr, (_core.Expr, type)):
            form=expr._pretty(self)
            if isinstance(form, stringPict):return form
            if isinstance(form, _legacyPict):return _from_legacy(form)
            return prettyForm(str(form))
        method=_composer(self, expr)
        if method is not None:return method(expr)
        try:
            text,baseline,binding=_core._pretty_form(expr, use_unicode=self._use_unicode,
                root_notation=self._settings.get('root_notation', True),
                use_unicode_sqrt_char=self._settings.get('use_unicode_sqrt_char', True))
        except TypeError as exc:
            if not _unconvertible(exc):raise
            return prettyForm(str(expr))
        return prettyForm(text, baseline, binding)

    def _print_seq(self, seq, left=None, right=None, delimiter=', ',
            parenthesize=lambda x: False, ifascii_nougly=True):
        from ._stringpict_ref import stringPict, prettyForm
        pforms = []
        for item in seq:
            pform = self._compose(item)
            if parenthesize(item):
                pform = prettyForm(*pform.parens())
            if pforms:
                pforms.append(delimiter)
            pforms.append(pform)
        if not pforms:
            s = stringPict('')
        else:
            s = prettyForm(*stringPict.next(*pforms))
        s = prettyForm(*s.parens(left, right, ifascii_nougly=ifascii_nougly))
        return s

    def _print_list(self, l):
        return self._print_seq(l, '[', ']')

    def _print_tuple(self, t):
        from ._stringpict_ref import stringPict, prettyForm
        if len(t) == 1:
            ptuple = prettyForm(*stringPict.next(self._compose(t[0]), ','))
            return prettyForm(*ptuple.parens('(', ')', ifascii_nougly=True))
        else:
            return self._print_seq(t, '(', ')')

    def _print_dict(self, d):
        from ._stringpict_ref import stringPict, prettyForm
        from symox import default_sort_key
        keys = sorted(d.keys(), key=default_sort_key)
        items = []
        for k in keys:
            K = self._compose(k)
            V = self._compose(d[k])
            s = prettyForm(*stringPict.next(K, ': ', V))
            items.append(s)
        return self._print_seq(items, '{', '}')

    def _print_set(self, s):
        from ._stringpict_ref import prettyForm
        from symox import default_sort_key
        if not s:
            return prettyForm('set()')
        items = sorted(s, key=default_sort_key)
        pretty = self._print_seq(items)
        pretty = prettyForm(*pretty.parens('{', '}', ifascii_nougly=True))
        return pretty

    def _print_frozenset(self, s):
        from ._stringpict_ref import stringPict, prettyForm
        from symox import default_sort_key
        if not s:
            return prettyForm('frozenset()')
        items = sorted(s, key=default_sort_key)
        pretty = self._print_seq(items)
        pretty = prettyForm(*pretty.parens('{', '}', ifascii_nougly=True))
        pretty = prettyForm(*pretty.parens('(', ')', ifascii_nougly=True))
        pretty = prettyForm(*stringPict.next(type(s).__name__, pretty))
        return pretty

    _print_mappingproxy = _print_dict

    # Control systems (sympy.physics.control), as in the reference.
    def _print_TransferFunction(self, expr):
        from symox import Mul, Pow
        if not expr.num == 1:
            num, den = expr.num, expr.den
            res = Mul(num, Pow(den, -1, evaluate=False), evaluate=False)
            return self._compose(res)
        else:
            return self._compose(1)/self._compose(expr.den)

    def _print_Series(self, expr):
        from ._stringpict_ref import prettyForm
        args = list(expr.args)
        for i, a in enumerate(expr.args):
            args[i] = prettyForm(*self._compose(a).parens())
        return prettyForm.__mul__(*args)

    def _print_MIMOSeries(self, expr):
        from ._stringpict_ref import prettyForm
        from symox.physics.control.lti import MIMOParallel
        args = list(expr.args)
        pretty_args = []
        for a in reversed(args):
            if (isinstance(a, MIMOParallel) and len(expr.args) > 1):
                expression = self._compose(a)
                expression.baseline = expression.height()//2
                pretty_args.append(prettyForm(*expression.parens()))
            else:
                expression = self._compose(a)
                expression.baseline = expression.height()//2
                pretty_args.append(expression)
        return prettyForm.__mul__(*pretty_args)

    def _print_Parallel(self, expr):
        from ._stringpict_ref import stringPict, prettyForm
        s = None
        for item in expr.args:
            pform = self._compose(item)
            if s is None:
                s = pform     # First element
            else:
                s = prettyForm(*stringPict.next(s))
                s.baseline = s.height()//2
                s = prettyForm(*stringPict.next(s, ' + '))
                s = prettyForm(*stringPict.next(s, pform))
        return s

    def _print_MIMOParallel(self, expr):
        from ._stringpict_ref import stringPict, prettyForm
        from symox.physics.control.lti import TransferFunctionMatrix
        s = None
        for item in expr.args:
            pform = self._compose(item)
            if s is None:
                s = pform     # First element
            else:
                s = prettyForm(*stringPict.next(s))
                s.baseline = s.height()//2
                s = prettyForm(*stringPict.next(s, ' + '))
                if isinstance(item, TransferFunctionMatrix):
                    s.baseline = s.height() - 1
                s = prettyForm(*stringPict.next(s, pform))
        return s

    def _print_Feedback(self, expr):
        from ._stringpict_ref import stringPict, prettyForm
        from symox.physics.control import TransferFunction, Series
        num, tf = expr.sys1, TransferFunction(1, 1, expr.var)
        num_arg_list = list(num.args) if isinstance(num, Series) else [num]
        den_arg_list = list(expr.sys2.args) if \
            isinstance(expr.sys2, Series) else [expr.sys2]
        if isinstance(num, Series) and isinstance(expr.sys2, Series):
            den = Series(*num_arg_list, *den_arg_list)
        elif isinstance(num, Series) and isinstance(expr.sys2, TransferFunction):
            if expr.sys2 == tf:
                den = Series(*num_arg_list)
            else:
                den = Series(*num_arg_list, expr.sys2)
        elif isinstance(num, TransferFunction) and isinstance(expr.sys2, Series):
            if num == tf:
                den = Series(*den_arg_list)
            else:
                den = Series(num, *den_arg_list)
        else:
            if num == tf:
                den = Series(*den_arg_list)
            elif expr.sys2 == tf:
                den = Series(*num_arg_list)
            else:
                den = Series(*num_arg_list, *den_arg_list)
        denom = prettyForm(*stringPict.next(self._compose(tf)))
        denom.baseline = denom.height()//2
        denom = prettyForm(*stringPict.next(denom, ' + ')) if expr.sign == -1 \
            else prettyForm(*stringPict.next(denom, ' - '))
        denom = prettyForm(*stringPict.next(denom, self._compose(den)))
        return self._compose(num)/denom

    def _print_MIMOFeedback(self, expr):
        from ._stringpict_ref import stringPict, prettyForm
        from symox.physics.control import MIMOSeries, TransferFunctionMatrix
        inv_mat = self._compose(MIMOSeries(expr.sys2, expr.sys1))
        plant = self._compose(expr.sys1)
        _feedback = prettyForm(*stringPict.next(inv_mat))
        _feedback = prettyForm(*stringPict.right("I + ", _feedback)) if expr.sign == -1 \
            else prettyForm(*stringPict.right("I - ", _feedback))
        _feedback = prettyForm(*stringPict.parens(_feedback))
        _feedback.baseline = 0
        _feedback = prettyForm(*stringPict.right(_feedback, '-1 '))
        _feedback.baseline = _feedback.height()//2
        _feedback = prettyForm.__mul__(_feedback, prettyForm(" "))
        if isinstance(expr.sys1, TransferFunctionMatrix):
            _feedback.baseline = _feedback.height() - 1
        _feedback = prettyForm(*stringPict.next(_feedback, plant))
        return _feedback

    def _print_TransferFunctionMatrix(self, expr):
        from ._stringpict_ref import prettyForm
        mat = self._compose(expr._expr_mat)
        mat.baseline = mat.height() - 1
        subscript = '\N{GREEK SMALL LETTER TAU}' if self._use_unicode else r'{t}'
        mat = prettyForm(*mat.right(subscript))
        return mat

    def _print_matrix_contents_rows(self, rows):
        """Reference ``_print_matrix_contents`` over a list of rows."""
        from ._stringpict_ref import prettyForm
        Ms = [[self._compose(c) for c in row] for row in rows]
        ncols = len(rows[0]) if rows else 0
        hsep, vsep = 2, 1
        maxw = [max([Ms[i][j].width() for i in range(len(rows))] or [0])
                for j in range(ncols)]
        D = None
        for i in range(len(rows)):
            D_row = None
            for j in range(ncols):
                s = Ms[i][j]
                left, right = center_pad(s.width(), maxw[j])
                s = prettyForm(*s.right(right))
                s = prettyForm(*s.left(left))
                if D_row is None:
                    D_row = s
                    continue
                D_row = prettyForm(*D_row.right(' '*hsep))
                D_row = prettyForm(*D_row.right(s))
            if D is None:
                D = D_row
                continue
            for _ in range(vsep):
                D = prettyForm(*D.below(' '))
            D = prettyForm(*D.below(D_row))
        if D is None:
            D = prettyForm('')
        return D

    def _print_DiagramGrid(self, grid):
        from symox import Symbol
        rows = [[grid[i, j] if grid[i, j] else Symbol(" ")
                 for j in range(grid.width)] for i in range(grid.height)]
        return self._print_matrix_contents_rows(rows)

    def _print(self, expr, **kwargs):
        from .stringpict import prettyForm
        if hasattr(expr, 'picture'):return expr
        if isinstance(expr, str):return prettyForm(expr)
        if hasattr(expr, '_pretty'):
            try:return expr._pretty(self)
            except TypeError:pass
        text=PrettyPrinter(self._settings).doprint(expr)
        return prettyForm(text)


_COMPOSERS=('list','tuple','dict','set','frozenset','mappingproxy','TransferFunction',
 'Series','MIMOSeries','Parallel','MIMOParallel','Feedback','MIMOFeedback',
 'TransferFunctionMatrix','DiagramGrid')


def _composer(printer, expr):
    """Reference ``_print_<Class>`` handler of a Python-level value."""
    for cls in type(expr).__mro__:
        if cls.__name__ in _COMPOSERS:
            if cls.__module__ not in ('builtins','mappingproxy') and \
                    not cls.__module__.endswith(('physics.control.lti',
                                                 'categories.diagram_drawing')):
                continue
            return getattr(printer, '_print_'+cls.__name__)
    return None


def _composed(printer, expr):
    return _composer(printer, expr) is not None


def _unconvertible(exc):
    return 'cannot convert object' in str(exc) or \
        getattr(exc, '_symox_unsupported_conversion', False)


def _from_legacy(form):
    from ._stringpict_ref import prettyForm
    return prettyForm('\n'.join(form.picture), form.baseline, getattr(form, 'binding', 0) or 0)


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
