"""Per-head native application classes for nominal undefined-function queries."""
import copyreg
from weakref import WeakValueDictionary
from . import Function
from . import core
from .core.function import UndefinedFunction
from ._node_types import NodeType,node_class


class _NativeFunctionType(UndefinedFunction,NodeType):
    def __new__(metaclass,name,bases,namespace):
        return type.__new__(metaclass,name,bases,namespace)

    def __call__(cls,*args,**kwargs):
        return cls._symox_function_head(*args,**kwargs)

    @property
    def _symox_head_name(cls):
        return cls.__name__


_types=WeakValueDictionary()
_NODE_HEADS=frozenset({'ComplexRootOf','RootSum'})


def application_type(name,*,nargs=None,**assumptions):
    key=(name,nargs,tuple(sorted(assumptions.items())))
    result=_types.get(key)
    if result is None:
        head=Function(name,nargs=nargs,**assumptions)
        if name in _NODE_HEADS:
            # Public node classes that symox represents as applied heads.
            bases,namespace=(node_class(name),),{'__module__':node_class(name).__module__}
        else:
            bases=(node_class('AppliedUndef'),)
            namespace={'__module__':None,'_sage_':core.function._undef_sage_helper}
        result=_NativeFunctionType(name,bases,{
            **namespace,'__slots__':(),
            '_symox_function_head':head,'nargs':head.nargs,
            'name':name,'_assumptions':head._assumptions,'_nargs':head._nargs,
        })
        _types[key]=result
    return result


def _restore_application_type(name,nargs,assumptions):
    return application_type(name,nargs=nargs,**assumptions)


copyreg.pickle(_NativeFunctionType,lambda cls:(_restore_application_type,(cls.__name__,cls._nargs,dict(cls._assumptions))))
