"""Expose real instance methods on constructor classes backed by native nodes."""
from inspect import getattr_static


def instance_method(native_type,name):
    if name.startswith('_'):
        raise AttributeError(name)
    descriptor=getattr_static(native_type,name)
    # Factory methods have to construct the requested public class. Delegating
    # those to another native class would silently change mutability or kind.
    if isinstance(descriptor,(staticmethod,classmethod)):
        raise AttributeError(name)
    method=getattr(native_type,name)
    if callable(method) and getattr(method,'__self__',None) is None:
        return method
    raise AttributeError(name)
