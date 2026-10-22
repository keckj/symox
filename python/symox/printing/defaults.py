"""Default str/repr printing mixin for non-expression objects."""


class DefaultPrinting:
    """Print through ``sstr`` for both ``str`` and ``repr``, as the reference mixin."""
    __slots__ = ()

    def __str__(self):
        from .str import sstr
        return sstr(self, order=None)

    __repr__ = __str__

    def _repr_disabled(self):
        return None

    def _repr_latex_(self):
        from .latex import latex
        s = latex(self, mode='plain')
        return "$\\displaystyle %s$" % s


Printable = DefaultPrinting
