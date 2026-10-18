"""Deferred rectangular comprehensions, expanded with native scalar entries."""
import inspect
import itertools
import math

from ._array import ImmutableDenseNDimArray,MAX_ENTRIES,MAX_RANK


class ArrayComprehension:
    __slots__ = ('_function','_limits')

    def __new__(cls,function,*symbols,**assumptions):
        from ... import sympify,Expr
        if assumptions:
            raise TypeError('Array comprehensions do not accept assumptions')
        result = object.__new__(cls)
        result._function = function if cls is ArrayComprehensionMap else sympify(function)
        if cls is not ArrayComprehensionMap and not isinstance(result._function,Expr):
            raise TypeError('A comprehension requires a scalar expression')
        if len(symbols)>MAX_RANK:
            raise ValueError('Array rank exceeds 32')
        limits = []
        for source in symbols:
            source = source.args if isinstance(source,Expr) else source
            if not isinstance(source,(tuple,list)) or len(source)!=3:
                raise ValueError('Comprehensions require (index, lower, upper) limits')
            index,lower,upper = map(sympify,source)
            if not index.is_Symbol or any(index==limit[0] for limit in limits):
                raise ValueError('Comprehension indices must be distinct symbols')
            if any(not isinstance(bound,Expr) or bound.is_integer is False
                   or bound.is_finite is False for bound in (lower,upper)):
                raise TypeError('Comprehension bounds must be integer expressions')
            if index in lower.free_symbols|upper.free_symbols:
                raise ValueError('A comprehension index cannot occur in its own bounds')
            if (upper-lower).is_negative is True:
                raise ValueError('Lower bound must not exceed upper bound')
            limits.append((index,lower,upper))
        result._limits = tuple(limits)
        return result

    @property
    def function(self):return self._function
    @property
    def limits(self):return self._limits
    @property
    def args(self):return (self.function,*self.limits)
    @property
    def func(self):return type(self)
    @property
    def variables(self):return [limit[0] for limit in self.limits]
    bound_symbols = variables
    @property
    def shape(self):return tuple(upper-lower+1 for _,lower,upper in self.limits)
    @property
    def is_shape_numeric(self):
        return all(bound.is_Integer for limit in self.limits for bound in limit[1:])
    @property
    def free_symbols(self):
        symbols = set() if isinstance(self,ArrayComprehensionMap) else self.function.free_symbols
        for index,lower,upper in self.limits:
            symbols = (symbols-{index})|lower.free_symbols|upper.free_symbols
        return symbols
    @property
    def kind(self):
        from ..._kind import scalar_kind,UndefinedKind
        from .kind import ArrayKind
        return ArrayKind(UndefinedKind if isinstance(self,ArrayComprehensionMap)
                         else scalar_kind(self.function))
    def rank(self):return len(self.limits)
    def __len__(self):
        if not self.is_shape_numeric:
            raise ValueError('Symbolic comprehension length cannot be evaluated')
        return math.prod(map(int,self.shape)) if self.shape else 0
    def _get_element(self,values):
        from ... import sympify
        return self.function.xreplace(dict(zip(self.variables,map(sympify,values))))
    def doit(self,**hints):
        if not self.is_shape_numeric:
            return self
        if len(self)>MAX_ENTRIES:
            raise ValueError('Comprehension expansion exceeds 100000 positions')
        ranges = [range(int(lower),int(upper)+1) for _,lower,upper in self.limits]
        values = [self._get_element(values) for values in itertools.product(*ranges)]
        return ImmutableDenseNDimArray(values,tuple(map(int,self.shape)))
    def tolist(self):
        if not self.is_shape_numeric:
            raise ValueError('A symbolic comprehension cannot be expanded to a list')
        return self.doit().tolist()
    def tomatrix(self):
        if not self.is_shape_numeric or self.rank()!=2:
            raise ValueError('Matrix conversion requires a concrete rank-two comprehension')
        return self.doit().tomatrix()
    def xreplace(self,mapping):
        if self in mapping:return mapping[self]
        function = self.function if isinstance(self,ArrayComprehensionMap) else self.function.xreplace(mapping)
        limits = [tuple(value.xreplace(mapping) for value in limit) for limit in self.limits]
        return type(self)(function,*limits)
    def subs(self,*args,**options):
        from ... import Dummy,Symbol,sympify
        from collections.abc import Mapping
        pairs = (list(args[0].items()) if len(args)==1 and isinstance(args[0],Mapping)
                 else list(args[0]) if len(args)==1 else [args])
        pairs = [(sympify(old),sympify(new)) for old,new in pairs]
        replacements = set().union(*(new.free_symbols for _,new in pairs)) if pairs else set()
        rename = {index:Dummy(index.name,integer=True) for index in self.variables if index in replacements}
        protected = [(old,new) for old,new in pairs if old not in self.variables]
        function = self.function if isinstance(self,ArrayComprehensionMap) else self.function.xreplace(rename).subs(protected,**options)
        limits = [(rename.get(index,index),lower.subs(pairs,**options),upper.subs(pairs,**options))
                  for index,lower,upper in self.limits]
        return type(self)(function,*limits)
    def __eq__(self,other):return type(self) is type(other) and self.args==other.args
    def __hash__(self):return hash((type(self),self.args))
    def __repr__(self):return '%s(%s)' % (type(self).__name__,', '.join(map(str,self.args)))
    def __reduce__(self):return type(self),self.args
    def to_sympy(self):
        from sympy.tensor.array import ArrayComprehension
        return ArrayComprehension(self.function.to_sympy(),
            *(tuple(value.to_sympy() for value in limit) for limit in self.limits))


class ArrayComprehensionMap(ArrayComprehension):
    __slots__ = ()
    def __new__(cls,function,*symbols,**assumptions):
        if not inspect.isfunction(function) or function.__name__!='<lambda>':
            raise ValueError('ArrayComprehensionMap requires a lambda')
        parameters = tuple(inspect.signature(function).parameters.values())
        if len(parameters)>1 or any(parameter.kind not in (
                inspect.Parameter.POSITIONAL_ONLY,inspect.Parameter.POSITIONAL_OR_KEYWORD)
                for parameter in parameters):
            raise ValueError('A comprehension lambda must accept zero or one argument')
        return super().__new__(cls,function,*symbols,**assumptions)
    def _get_element(self,values):
        return self.function() if not inspect.signature(self.function).parameters else self.function(math.prod(values))
    def to_sympy(self):
        from sympy.tensor.array import ArrayComprehensionMap
        return ArrayComprehensionMap(self.function,
            *(tuple(value.to_sympy() for value in limit) for limit in self.limits))


__all__ = ['ArrayComprehension','ArrayComprehensionMap']
