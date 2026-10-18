"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox.utilities.cache', 'cacheit clear_cache')
_export(globals(), 'functools', 'cached_property')


def print_cache():
    from ..utilities.cache import CACHE
    for function in CACHE:print(function.__name__,function.cache_info())
def lazy_function(module,name):
    import importlib
    from functools import wraps
    target=None
    def wrapped(*args,**kwargs):
        nonlocal target
        if target is None:target=getattr(importlib.import_module(module),name)
        return target(*args,**kwargs)
    wrapped.__name__=name
    return wrapped
__all__ += ['print_cache','lazy_function']
