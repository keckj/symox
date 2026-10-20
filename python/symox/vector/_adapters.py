"""Dispatch native scalar arithmetic and public transforms for vector objects."""
from functools import wraps
from .._dispatch_arguments import FirstArgument


def install():
    import symox
    from ._basis import _BasisDependent
    Expr = symox._core.Expr
    if getattr(Expr,'_coordinate_vector_dispatch_installed',False):
        return
    for name in ('__add__','__sub__','__mul__','__truediv__'):
        original = getattr(Expr,name)
        def arithmetic(original,name):
            @wraps(original)
            def operation(self,other):
                if isinstance(other,_BasisDependent):
                    if name == '__add__':
                        return other.__radd__(self)
                    if name == '__sub__':
                        return other.__rsub__(self)
                    if name == '__mul__':
                        return other.__rmul__(self)
                    raise TypeError('division by a vector or dyadic is undefined')
                return original(self,other)
            return operation
        setattr(Expr,name,arithmetic(original,name))
    for name in ('diff','integrate','simplify','expand','trigsimp'):
        original = getattr(symox,name)
        def transform(original,name):
            first=FirstArgument(original)
            @wraps(original)
            def operation(*args,**kwargs):
                expression=first.get(args,kwargs)
                if isinstance(expression,_BasisDependent):
                    expression,rest,options=first.split(args,kwargs)
                    return getattr(expression,name)(*rest,**options)
                return original(*args,**kwargs)
            return operation
        setattr(symox,name,transform(original,name))
    Expr._coordinate_vector_dispatch_installed = True
