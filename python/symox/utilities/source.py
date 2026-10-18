"""Resolve explicitly requested dotted Python objects."""
from importlib import import_module

def get_mod_func(callback):
    if not isinstance(callback,str):return callback.__module__,callback.__name__
    return tuple(callback.rsplit('.',1)) if '.' in callback else (callback,'')
def get_class(lookup_view):
    if callable(lookup_view):return lookup_view
    module,attribute=get_mod_func(lookup_view);value=import_module(module)
    return getattr(value,attribute) if attribute else value
__all__=['get_mod_func','get_class']
