"""Default printing mixin for non-expression classes (SymPy's ``Printable``)."""


class Printable:
    """Mix-in giving ``str``/``repr`` through ``sstr`` and Jupyter LaTeX output."""
    __slots__ = ()

    def __str__(self):
        from .str import sstr
        return sstr(self, order=None)

    __repr__ = __str__

    def _repr_disabled(self):
        return None

    _repr_png_ = _repr_disabled
    _repr_svg_ = _repr_disabled

    def _repr_latex_(self):
        from .latex import latex
        s = latex(self, mode='plain')
        return "$\\displaystyle %s$" % s


DefaultPrinting = Printable
__all__ = ['Printable', 'DefaultPrinting']
