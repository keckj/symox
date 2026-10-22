"""Multiple dispatch with subclass specificity and explicit ambiguity errors."""
from functools import update_wrapper


class Dispatcher:
    def __init__(self,name,doc=None):
        self.__name__,self.__doc__ = name,doc
        self.funcs = {}
        self._cache = {}
    def add(self,signature,func):
        signature = tuple(signature)
        if any(not isinstance(t,(type,tuple)) for t in signature):raise TypeError('dispatch signatures contain types or tuples of types')
        self.funcs[signature] = func;self._cache.clear()
    def register(self,*types,**kwargs):
        def decorate(func):self.add(types,func);return func
        return decorate
    def dispatch(self,*types):
        if types in self._cache:return self._cache[types]
        candidates = [s for s in self.funcs if len(s) == len(types) and all(issubclass(t,spec) for t,spec in zip(types,s))]
        def more_specific(first,second):
            def subclass(a,b):return all(issubclass(t,b) for t in (a if isinstance(a,tuple) else (a,)))
            return all(subclass(a,b) for a,b in zip(first,second)) and not all(subclass(b,a) for a,b in zip(first,second))
        best = [s for s in candidates if not any(more_specific(other,s) for other in candidates)]
        if len(best) > 1:raise TypeError('ambiguous dispatch of %s for %r' % (self.__name__,types))
        function = None if not best else self.funcs[best[0]]
        self._cache[types] = function
        return function
    def __call__(self,*args,**kwargs):
        function = self.dispatch(*(type(a) for a in args))
        if function is None:raise NotImplementedError('no matching signature for '+self.__name__)
        return function(*args,**kwargs)
    def __get__(self,instance,owner):
        if instance is None:return self
        from functools import partial
        return partial(self,instance)


_namespace = {}


def dispatch(*types,**kwargs):
    namespace = kwargs.pop('namespace',_namespace)
    if kwargs:raise TypeError('unknown dispatch options')
    def decorate(func):
        name = func.__name__
        dispatcher = namespace.get(name)
        if dispatcher is None:dispatcher = namespace[name] = Dispatcher(name,func.__doc__)
        if not isinstance(dispatcher,Dispatcher):raise TypeError('dispatch name is already occupied')
        dispatcher.add(types,func)
        return dispatcher
    return decorate


class MDNotImplementedError(NotImplementedError):
    """A dispatched implementation declines; the next match is tried."""


__all__ = ['dispatch','Dispatcher','MDNotImplementedError']
