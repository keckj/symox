"""Python-specific code generation nodes (port of ``sympy.codegen.pynodes``)."""
from .abstract_nodes import List as AbstractList
from .ast import Token


class List(AbstractList):
    # The native Tuple rebuilds plain tuples; SymPy keeps the List class
    # through ``func`` and the tree transformations, and prints it as a list.

    @property
    def func(self):
        return type(self)

    def _rewrap(self, result):
        from ..core.containers import Tuple
        if type(result) is Tuple:
            return type(self)(*result.args)
        return result

    def replace(self, *args, **kwargs):
        return self._rewrap(super().replace(*args, **kwargs))

    def subs(self, *args, **kwargs):
        return self._rewrap(super().subs(*args, **kwargs))

    def xreplace(self, *args, **kwargs):
        return self._rewrap(super().xreplace(*args, **kwargs))

    def _sympystr(self, printer, *args, **kwargs):
        # StrPrinter._print_List
        return '[%s]' % ', '.join(printer._print(a) for a in self.args)

    def _sympyrepr(self, printer, *args, **kwargs):
        # ReprPrinter: the class name and the srepr of the arguments.
        return 'List(%s)' % ', '.join(printer._print(a) for a in self.args)

    def __str__(self):
        from ..printing.str import sstr
        return sstr(self)

    def __repr__(self):
        return str(self)


class NumExprEvaluate(Token):
    """represents a call to :class:`numexpr`s :func:`evaluate`"""
    _fields = ('expr',)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # SymPy's Token keeps a field without a constructor as given.
        if args:
            object.__setattr__(self, 'expr', args[0])
        elif 'expr' in kwargs:
            object.__setattr__(self, 'expr', kwargs['expr'])

    def __str__(self):
        from ..printing.str import sstr
        return 'NumExprEvaluate(%s)' % sstr(self.expr)

    def __repr__(self):
        from ..printing.repr import srepr
        return 'NumExprEvaluate(%s)' % srepr(self.expr)


__all__ = ['List', 'NumExprEvaluate']
