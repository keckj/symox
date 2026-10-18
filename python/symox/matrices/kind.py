"""Interned kinds of matrix-valued objects."""
from .._kind import Kind,NumberKind

class MatrixKind(Kind):
    __slots__=('element_kind',)
    _instances={}
    def __new__(cls,element_kind=NumberKind):
        if element_kind not in cls._instances:
            obj=super().__new__(cls);obj.element_kind=element_kind;cls._instances[element_kind]=obj
        return cls._instances[element_kind]
    def __repr__(self):return 'MatrixKind(%r)'%self.element_kind
    def __reduce__(self):return type(self),(self.element_kind,)

__all__=['MatrixKind']
