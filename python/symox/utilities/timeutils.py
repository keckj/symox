from .. import timed
__all__ = ['timed']

from functools import wraps as _wraps
from time import perf_counter as _clock
from .misc import debug as _debug

def timethis(name):
    def decorate(function):
        @_wraps(function)
        def wrapped(*args,**kwargs):
            start=_clock()
            try:return function(*args,**kwargs)
            finally:_debug(name,function.__name__,_clock()-start)
        return wrapped
    return decorate
__all__ += ['timethis']
