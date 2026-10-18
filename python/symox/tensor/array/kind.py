"""Kinds for concrete and deferred multidimensional arrays."""
from ..._kind import Kind,NumberKind


class ArrayKind(Kind):
    __slots__ = ('element_kind',)
    _instances = {}
    def __new__(cls,element_kind=NumberKind):
        if not isinstance(element_kind,Kind):
            raise TypeError('An array element kind must be a Kind')
        if element_kind not in cls._instances:
            result = object.__new__(cls)
            result.element_kind = element_kind
            cls._instances[element_kind] = result
        return cls._instances[element_kind]
    def __repr__(self):return 'ArrayKind(%r)' % self.element_kind
    def __reduce__(self):return type(self),(self.element_kind,)
    def to_sympy(self):
        from sympy.tensor.array import ArrayKind
        return ArrayKind(self.element_kind.to_sympy())


__all__ = ['ArrayKind']
