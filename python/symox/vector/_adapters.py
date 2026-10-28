"""Dispatch native scalar arithmetic and public transforms for vector objects."""
from functools import wraps
from .._dispatch_arguments import FirstArgument


_SCALAR_DISPATCH_PROOFS = ()


def _dispatch_states(operation):
    from types import FunctionType
    pending=[operation];seen=set();states=[]
    while pending:
        function=pending.pop()
        if type(function) is not FunctionType or id(function) in seen:continue
        seen.add(id(function));closure=function.__closure__
        contents=tuple(cell.cell_contents for cell in closure or ())
        states.append((function,function.__code__,closure,contents,function.__defaults__,function.__kwdefaults__))
        pending.extend(value for value in contents if type(value) is FunctionType)
    return tuple(states)


def install():
    import symox
    import builtins
    global _SCALAR_DISPATCH_PROOFS
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
        operation=arithmetic(original,name)
        proof=(operation,original,_dispatch_states(operation),operation.__globals__,
               _BasisDependent,_BasisDependent.__mro__,builtins.isinstance,name)
        operation._symox_vector_proof=proof
        _SCALAR_DISPATCH_PROOFS+=(proof,)
        setattr(Expr,name,operation)
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
    install._stock_proofs=_SCALAR_DISPATCH_PROOFS
    Expr._coordinate_vector_dispatch_installed = True


# Canonical templates, fixed before the package calls install(). Runtime
# replacement of a wrapper or the installer cannot create a new certificate.
from types import CodeType as _CodeType
_SCALAR_INSTALL_CODE = install.__code__
_SCALAR_FACTORY_CODE = next(code for code in _SCALAR_INSTALL_CODE.co_consts
                            if type(code) is _CodeType and code.co_name=='arithmetic')
_SCALAR_OPERATION_CODE = next(code for code in _SCALAR_FACTORY_CODE.co_consts
                              if type(code) is _CodeType and code.co_name=='operation')
