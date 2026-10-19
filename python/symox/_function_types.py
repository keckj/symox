"""Per-head native application classes for nominal undefined-function queries."""
from weakref import WeakValueDictionary
from . import Function,_core
from .core.function import UndefinedFunction


class _NativeFunctionType(UndefinedFunction):
    def __new__(metaclass,name,bases,namespace):
        return type.__new__(metaclass,name,bases,namespace)

    def __call__(cls,*args,**kwargs):
        return cls._symox_function_head(*args,**kwargs)

    @property
    def _symox_head_name(cls):
        return cls.__name__


_types=WeakValueDictionary()


def application_type(name,*,nargs=None,**assumptions):
    key=(name,nargs,tuple(sorted(assumptions.items())))
    result=_types.get(key)
    if result is None:
        head=Function(name,nargs=nargs,**assumptions)
        result=_NativeFunctionType(name,(_core.Expr,),{
            '__module__':'symox','__slots__':(),
            '_symox_function_head':head,'nargs':head.nargs,
            'name':name,'_assumptions':head._assumptions,'_nargs':head._nargs,
        })
        _types[key]=result
    return result
