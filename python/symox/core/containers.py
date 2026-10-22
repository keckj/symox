from .. import Dict,Tuple
__all__ = ['Dict','Tuple']

from collections.abc import MutableSet
from .._kind import Kind

class OrderedSet(MutableSet):
    def __init__(self,iterable=()):self._data=dict.fromkeys(iterable)
    def __contains__(self,value):return value in self._data
    def __iter__(self):return iter(self._data)
    def __len__(self):return len(self._data)
    def add(self,key):self._data[key]=None
    def discard(self,key):self._data.pop(key,None)
    def update(self,iterable):self._data.update(dict.fromkeys(iterable))
    def pop(self,last=True):
        if not self._data:raise KeyError('pop from empty OrderedSet')
        value=next(reversed(self._data)) if last else next(iter(self._data));del self._data[value];return value
    def __repr__(self):return 'OrderedSet(%r)'%list(self)
    def intersection(self,other):return self.__class__([val for val in self if val in other])
    def difference(self,other):return self.__class__([val for val in self if val not in other])

class TupleKind(Kind):
    _instances={}
    def __new__(cls,*element_kind):
        if element_kind not in cls._instances:
            obj=super().__new__(cls);obj.element_kind=element_kind;cls._instances[element_kind]=obj
        return cls._instances[element_kind]
    def __repr__(self):return 'TupleKind('+', '.join(map(str,self.element_kind))+')'

def tuple_wrapper(method):
    from functools import wraps
    @wraps(method)
    def wrapped(*args,**kwargs):return method(*(Tuple(*arg) if isinstance(arg,(tuple,list)) else arg for arg in args),**kwargs)
    return wrapped
__all__ += ['OrderedSet','TupleKind','tuple_wrapper']
