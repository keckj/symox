"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox', 'S')

from .._compat import _SingletonRegistry, _ExprType
from .._basic import Basic, _BasicMeta
from .. import S

SingletonRegistry=type(S)


# SymPy's registry installs singleton classes lazily: the instance is created
# and stored on S at the first access of S.<name>.
S.__dict__.setdefault('_classes_to_install',{})

def _register(cls):
    # Make sure a duplicate class overwrites the old one
    if cls.__name__ in vars(S):delattr(S,cls.__name__)
    S._classes_to_install[cls.__name__]=cls
def _unregister(cls):
    S._classes_to_install.pop(cls.__name__,None)
    if cls.__name__ in vars(S):delattr(S,cls.__name__)

_registry_getattr=SingletonRegistry.__getattr__
def _getattr(self,name):
    pending=self.__dict__.get('_classes_to_install',{})
    if name not in pending:return _registry_getattr(self,name)
    value=pending[name]()
    self.__setattr__(name,value)
    del pending[name]
    return value

# The public registry is shared by the constructors and native scalar facade.
if not hasattr(type(S),'register'):
    type(S).register=staticmethod(_register)
    type(S).__getattr__=_getattr
if not hasattr(type(S),'unregister'):type(S).unregister=staticmethod(_unregister)

class Singleton(_ExprType,_BasicMeta):
    """
    Metaclass for singleton classes.

    A singleton class has only one instance which is returned every time the
    class is instantiated. Additionally, this instance can be accessed through
    the global registry object ``S`` as ``S.<class_name>``. The instance is
    created when the class is created, and installed on ``S`` at its first
    access, as in SymPy.
    """
    def __init__(cls,*args,**kwargs):
        cls._instance=obj=Basic.__new__(cls)
        cls.__new__=lambda cls:obj
        cls.__getnewargs__=lambda obj:()
        cls.__getstate__=lambda obj:None
        S.register(cls)

__all__ += ['Singleton','SingletonRegistry']
