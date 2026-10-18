"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox._kind', 'Kind')

from .._kind import UndefinedKind

from .._kind import KindMeta

class KindDispatcher:
    def __init__(self,name,commutative=False,doc=None):self.name,self.commutative,self.doc=name,commutative,doc;self._rules=[]
    def register(self,*types,**kwargs):
        if len(types)!=2:raise RuntimeError('kind dispatch requires two classes')
        def decorator(function):self._rules.append((types,function));return function
        return decorator
    def _pair(self,left,right):
        matches=[]
        for (first,second),function in self._rules:
            if isinstance(left,first) and isinstance(right,second):matches.append((len(type(left).__mro__)-type(left).__mro__.index(first)+len(type(right).__mro__)-type(right).__mro__.index(second),function,left,right))
            if self.commutative and isinstance(right,first) and isinstance(left,second):matches.append((0,function,right,left))
        if not matches:return left if left is right else UndefinedKind
        _,function,left,right=max(matches,key=lambda entry:entry[0]);return function(left,right)
    def __call__(self,*kinds):
        if not kinds:return UndefinedKind
        result=kinds[0]
        for kind in kinds[1:]:result=self._pair(result,kind)
        return result
    def __repr__(self):return '<dispatched '+self.name+'>'
__all__ += ['KindMeta','KindDispatcher']
