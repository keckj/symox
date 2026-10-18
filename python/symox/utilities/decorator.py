from .. import public,threaded,xthreaded
__all__ = ['public','threaded','xthreaded']

# Compatible imports reuse the existing native objects.
from symox._api_exports import export as _export
_export(globals(), 'symox.utilities.cache', 'memoize_property')

from contextlib import contextmanager as _contextmanager
from functools import wraps as _wraps
from typing import TypeVar as _TypeVar
from .._public_tools import _threaded
T=_TypeVar('T')

def threaded_factory(func,use_add):return _threaded(func,use_add)

def conserve_mpmath_dps(func):
    @_wraps(func)
    def wrapped(*args,**kwargs):
        import mpmath
        precision=mpmath.mp.dps
        try:return func(*args,**kwargs)
        finally:mpmath.mp.dps=precision
    return wrapped

def deprecated(message,*,deprecated_since_version,active_deprecations_target,stacklevel=3):
    def decorate(function):
        @_wraps(function)
        def wrapped(*args,**kwargs):
            from .exceptions import sympy_deprecation_warning
            sympy_deprecation_warning(message,deprecated_since_version=deprecated_since_version,
                active_deprecations_target=active_deprecations_target,stacklevel=stacklevel)
            return function(*args,**kwargs)
        return wrapped
    return decorate

def doctest_depends_on(exe=None,modules=None,disable_viewers=None,python_version=None,ground_types=None):
    def decorate(function):
        function._doctest_depends_on={'executables':exe or [],'modules':modules or [],
            'disable_viewers':disable_viewers or [],'python_version':python_version,'ground_types':ground_types}
        return function
    return decorate

class no_attrs_in_subclass:
    def __init__(self,cls,f):self.cls,self.f=cls,f
    def __get__(self,instance,owner):
        if owner is not self.cls:raise AttributeError(getattr(self.f,'__name__','attribute'))
        return self.f.__get__(instance,owner) if hasattr(self.f,'__get__') else self.f
__all__ += ['T','threaded_factory','conserve_mpmath_dps','deprecated','doctest_depends_on','no_attrs_in_subclass']
