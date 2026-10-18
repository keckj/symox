"""Interned scalar kinds used by typed array metadata."""


class Kind:
    __slots__ = ()


class _ScalarKind(Kind):
    __slots__ = ('name',)
    def __init__(self,name):self.name=name
    def __repr__(self):return self.name
    def __reduce__(self):return _scalar_kind,(self.name,)
    def to_sympy(self):
        from sympy.core import kind
        return getattr(kind,self.name)


NumberKind = _ScalarKind('NumberKind')
BooleanKind = _ScalarKind('BooleanKind')
UndefinedKind = _ScalarKind('UndefinedKind')


def _scalar_kind(name):
    return {'NumberKind':NumberKind,'BooleanKind':BooleanKind,
            'UndefinedKind':UndefinedKind}[name]


def scalar_kind(value,_depth=0):
    from . import true,false
    if _depth>128:return UndefinedKind
    if value in (true,false) or getattr(value,'is_Relational',False) or getattr(value,'is_Boolean',False):
        return BooleanKind
    if (not value.args or value.is_Number) and value.is_commutative is True:
        return NumberKind
    if getattr(value.func,'__name__','')=='exp' and value.is_commutative is True:
        return NumberKind
    if value.is_Add or value.is_Mul or value.is_Pow:
        kinds={scalar_kind(arg,_depth+1) for arg in value.args}
        return next(iter(kinds)) if len(kinds)==1 else UndefinedKind
    if value.is_Integral or value.is_Sum or value.is_Product:
        return scalar_kind(value.function,_depth+1)
    return UndefinedKind


__all__ = ['Kind','NumberKind','BooleanKind','UndefinedKind']
