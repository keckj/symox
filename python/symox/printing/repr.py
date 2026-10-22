"""Native constructor representation and recursive container formatting."""
from ._legacy import LegacyPrinter as Printer

class ReprPrinter(Printer):
    printmethod='_sympyrepr'

    def _print(self, expr, **kwargs):
        if isinstance(expr, type):
            # function classes (sin, user Function subclasses) print by name
            return expr.__name__
        from symox import Function
        if type(expr) is Function:
            # an undefined function head, in the form of its applied srepr
            parts = [repr(expr.name)]
            parts += ['%s=%s' % item for item in sorted(expr._assumptions.items())]
            if expr._nargs is not None:
                parts.append('nargs=(%s)' % ''.join('%d,' % n for n in expr._nargs))
            return 'Function(%s)' % ', '.join(parts)
        return super()._print(expr, **kwargs)

    def _print_Matrix(self, expr):
        if expr.shape == (0, 0):
            arguments = '[]'
        elif not expr.rows or not expr.cols:
            arguments = '%s, %s, []' % expr.shape
        else:
            arguments = '[' + ', '.join('[' + ', '.join(
                self._print(expr[i, j]) for j in range(expr.cols)) + ']'
                for i in range(expr.rows)) + ']'
        name = 'ImmutableDenseMatrix' if expr._immutable else 'MutableDenseMatrix'
        return name + '(' + arguments + ')'

    def emptyPrinter(self,expr):
        from symox import _core
        if expr is None or isinstance(expr, (bool, int, float, complex)):
            return repr(expr)
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
