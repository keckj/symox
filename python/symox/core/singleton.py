"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox', 'S')

from .._compat import _SingletonRegistry
from .. import S

SingletonRegistry=type(S)


def _register(cls):
    setattr(S,cls.__name__,cls())
def _unregister(cls):
    if cls.__name__ in vars(S):delattr(S,cls.__name__)

# The public registry is shared by the constructors and native scalar facade.
if not hasattr(type(S),'register'):type(S).register=staticmethod(_register)
if not hasattr(type(S),'unregister'):type(S).unregister=staticmethod(_unregister)

class Singleton(type):
    def __new__(mcls,name,bases,namespace,**kwargs):
        cls=super().__new__(mcls,name,bases,namespace,**kwargs);cls._singleton_instance=None
        return cls
    def __call__(cls,*args,**kwargs):
        if args or kwargs:raise TypeError('singletons do not take arguments')
        if cls._singleton_instance is None:
            instance=super().__call__();cls._singleton_instance=instance;setattr(S,cls.__name__,instance)
        return cls._singleton_instance
    def __instancecheck__(cls,value):return type(value) is cls

__all__ += ['Singleton','SingletonRegistry']
