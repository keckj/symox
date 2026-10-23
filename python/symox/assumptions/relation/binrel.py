"""User-defined binary relations with native Boolean lowering and dispatch."""
from ..assume import Predicate,AppliedPredicate,_PredicateType,_AppliedType
from ... import _core

class _BinaryType(_PredicateType):
    def __instancecheck__(cls,value):
        if cls in type(value).__mro__:return True
        return cls is BinaryRelation and isinstance(value,_core.Expr) and value.class_name in ('EqualityPredicate','UnequalityPredicate','StrictGreaterThanPredicate','StrictLessThanPredicate','GreaterThanPredicate','LessThanPredicate')

class BinaryRelation(Predicate,metaclass=_BinaryType):
    is_reflexive=None
    is_symmetric=None
    def __new__(cls,*args,**kwargs):return object.__new__(cls)
    def __init__(self,name=None):
        self.name=str(name if name is not None else getattr(type(self),'name',type(self).__name__))
        self._native=_core._predicate(self.name,undefined=True);self._handlers=[]
        from ..._public_tools import register_handler
        register_handler(self.name,self._eval_handler)
    def _symox_(self):return self._native
    def _eval_handler(self,*arguments,assumptions=True):return self.eval(arguments,assumptions)
    def __call__(self,*args):
        if len(args)!=2:raise ValueError('a binary relation takes exactly two arguments')
        return AppliedBinaryRelation(self,*args)
    def register(self,*types):
        if len(types)!=2:raise ValueError('binary relation dispatch needs two types')
        def decorate(function):self._handlers.append((types,function));return function
        return decorate
    def eval(self,args,assumptions=True):
        if args[0]==args[1] and self.is_reflexive is not None:return self.is_reflexive
        for types,handler in reversed(self._handlers):
            if all(isinstance(arg,typ) for arg,typ in zip(args,types)):
                value=handler(*args,assumptions=assumptions)
                if value is not None:return value
        return None
    reversed=property(lambda self:self if self.is_symmetric else None)
    negated=None
    def __repr__(self):return 'Q.'+self.name
    __str__=__repr__
    def __hash__(self):return hash(self._native)
    def __eq__(self,other):return isinstance(other,BinaryRelation) and self.name==other.name


class _AppliedBinaryType(_AppliedType):
    def __instancecheck__(cls,value):return cls in type(value).__mro__ or isinstance(value,_core.Expr) and value.class_name=='AppliedBinaryRelation'
class AppliedBinaryRelation(AppliedPredicate,metaclass=_AppliedBinaryType):
    def __new__(cls,predicate,*args):
        if len(args)!=2:raise ValueError('a binary relation takes exactly two arguments')
        if isinstance(predicate,_core.Expr):return _core._applied_predicate((predicate,*args))
        value=object.__new__(cls)
        from symox import sympify
        value.function=predicate;value.arguments=tuple(sympify(arg) for arg in args);return value
    class_name='AppliedBinaryRelation'
    is_Boolean=True
    args=property(lambda self:(self.function,*self.arguments))
    lhs=property(lambda self:self.arguments[0])
    rhs=property(lambda self:self.arguments[1])
    free_symbols=property(lambda self:self.lhs.free_symbols|self.rhs.free_symbols)
    def _symox_(self):return _core._applied_predicate((self.function._symox_(),*self.arguments))
    @property
    def reversed(self):
        relation=self.function.reversed
        return relation(self.rhs,self.lhs) if relation is not None else self
    @property
    def reversedsign(self):
        relation=self.function.reversed
        return relation(-self.lhs,-self.rhs) if relation is not None else self
    @property
    def negated(self):
        from symox import Not
        return self.function.negated(*self.arguments) if self.function.negated is not None else Not(self._symox_())
    def __bool__(self):
        from symox import ask
        result=ask(self)
        if result is None:raise TypeError('truth value of relation is undecidable')
        return result
    def __repr__(self):return str(self.function)+'('+', '.join(map(str,self.arguments))+')'
    def __hash__(self):return hash(self._symox_())
    def __eq__(self,other):return isinstance(other,AppliedBinaryRelation) and self.args==other.args
    def subs(self,*args,**kwargs):return type(self)(self.function,*(arg.subs(*args,**kwargs) for arg in self.arguments))

__all__=['BinaryRelation','AppliedBinaryRelation']
