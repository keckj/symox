"""Bounded typed memoization with explicit cache invalidation."""
from functools import lru_cache, wraps

CACHE=[]


def cacheit(function):
    cached=lru_cache(maxsize=1000,typed=True)(function)
    @wraps(function)
    def wrapper(*args,**kwargs):
        try:
            hash((args,tuple(kwargs.items())))
        except TypeError:
            return function(*args,**kwargs)
        return cached(*args,**kwargs)
    wrapper.cache_info=cached.cache_info
    wrapper.cache_clear=cached.cache_clear
    CACHE.append(wrapper)
    return wrapper


def clear_cache():
    for function in CACHE:
        function.cache_clear()


def memoize_property(propfunc):
    attribute='_'+propfunc.__name__
    sentinel=object()
    @wraps(propfunc)
    def accessor(self):
        result=getattr(self,attribute,sentinel)
        if result is sentinel:
            result=propfunc(self)
            setattr(self,attribute,result)
        return result
    return property(accessor)


__all__=['cacheit','memoize_property']
