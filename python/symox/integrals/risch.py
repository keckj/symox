"""Formal non-elementary integral marker, using native held integrals.

This constructor does not prove non-elementarity. The native integration
engine does not yet implement the complete Risch decision procedure and must
not use this marker to label an unsupported integral as proved non-elementary.
"""
from symox import Integral, _core


class NonElementaryIntegral(Integral, _core.Expr):
    """An explicitly constructed held integral with its marker type retained."""
    __slots__ = ()
    is_NonElementaryIntegral = True

    def __new__(cls, function, *symbols, **assumptions) -> 'Integral':
        node = Integral(function, *symbols, **assumptions)
        return _core.Expr.__new__(cls, node)

    @property
    def func(self):
        return type(self)

    def __getnewargs__(self):
        return (self.function, *self.limits)

    def __reduce_ex__(self, protocol):
        return type(self), self.__getnewargs__()

    def _retain_marker(self, result):
        return _core.Expr.__new__(type(self), result) if result.is_Integral else result

    def subs(self, *args, **kwargs):
        return self._retain_marker(super().subs(*args, **kwargs))

    def xreplace(self, rule):
        return self._retain_marker(super().xreplace(rule))

    def doit(self, **hints):
        return self._retain_marker(super().doit(**hints))

    def simplify(self, **kwargs):
        return self._retain_marker(super().simplify(**kwargs))

    def expand(self, *args, **kwargs):
        return self._retain_marker(super().expand(*args, **kwargs))

    def copy(self):
        return self

    def srepr(self):
        from symox import srepr, Tuple
        return type(self).__name__ + '(' + ', '.join(
            srepr(arg) for arg in (self.function, *(Tuple(*limit) for limit in self.limits))) + ')'

    def to_sympy(self):
        from sympy.integrals.risch import NonElementaryIntegral as ReferenceIntegral
        return ReferenceIntegral(self.function.to_sympy(),
            *(tuple(item.to_sympy() for item in limit) for limit in self.limits))


__all__ = ['NonElementaryIntegral']
