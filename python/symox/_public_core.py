"""Remaining core constructors and context-local Python evaluation controls."""
from collections.abc import Mapping
from contextvars import ContextVar
from functools import wraps
import math
import operator

from . import _core
from ._compat import _ExprType

_evaluation = ContextVar('symox_evaluation',default=True)
_divide_error = ContextVar('symox_divide_error',default=False)
_wild_functions = {}


class SympifyError(ValueError,TypeError):
    def __init__(self,expr,base_exc=None):
        self.expr,self.base_exc = expr,base_exc
        super().__init__('cannot sympify %r%s' % (expr,': '+str(base_exc) if base_exc else ''))


class PrecisionExhausted(ArithmeticError):
    """A requested numerical accuracy could not be attained."""


class Derivative(metaclass=_ExprType):
    """Held native derivative, with the reference constructor keyword schema."""
    _matches=staticmethod(lambda expr:expr.is_Derivative)

    def __new__(cls,expr,*variables,**kwargs):
        evaluate=kwargs.get('evaluate',False)
        if evaluate and not kwargs.get('simplify',True):
            raise NotImplementedError('evaluated derivatives use native canonical simplification')
        return _core.Derivative(expr,*variables,evaluate=bool(evaluate))


class evaluate:
    def __init__(self,x):
        if not isinstance(x,bool):
            raise TypeError('evaluation state must be boolean')
        self.x,self._tokens = x,[]

    def __enter__(self):
        self._tokens.append(_evaluation.set(self.x))
        return self

    def __exit__(self,*exc):
        _evaluation.reset(self._tokens.pop())


def seterr(divide=False):
    if not isinstance(divide,bool):
        raise TypeError('divide must be boolean')
    _divide_error.set(divide)


class _FunctionClassType(type):
    def __instancecheck__(cls,value):
        from . import Function,_is_function_head
        from ._host import is_host_function_class
        return (value is Function or _is_function_head(value) or is_host_function_class(value)
                or hasattr(value,'_symox_head_name') or hasattr(value,'_heads'))


class FunctionClass(type,metaclass=_FunctionClassType):
    """Metaclass interface for native function-head factories."""


class WildFunction(metaclass=_ExprType):
    _matches = staticmethod(lambda expr:expr in _wild_functions)

    def __new__(cls,name,**assumptions):
        counts = assumptions.pop('nargs',None)
        if counts is not None:
            try:
                counts = (operator.index(counts),)
            except TypeError:
                counts = tuple(operator.index(n) for n in counts)
            if not counts or any(n < 0 for n in counts):
                raise ValueError('function argument counts must be nonnegative')
            counts = tuple(sorted(set(counts)))
        node = _core.Wild(name,(),(),**assumptions)
        _wild_functions[node] = counts
        return node


class Dict(Mapping):
    """An immutable mapping with native symbolic keys and values."""
    def __init__(self,*args):
        if len(args) == 1 and isinstance(args[0],Mapping):
            entries = args[0].items()
        else:
            entries = args
        from . import sympify
        data = {}
        for pair in entries:
            if getattr(pair,'is_Tuple',False):
                pair = pair.args
            if len(pair) != 2:
                raise ValueError('dictionary entries must be key/value pairs')
            key,value = pair
            data[sympify(key)] = sympify(value)
        from types import MappingProxyType
        self._data = MappingProxyType(data)

    def __getitem__(self,key):
        from . import sympify
        return self._data[sympify(key)]

    def __iter__(self):
        return iter(self._data)

    def __len__(self):
        return len(self._data)

    def __hash__(self):
        return hash(frozenset(self._data.items()))

    def __repr__(self):
        return '{'+', '.join('%s: %s' % pair for pair in sorted(self._data.items(),key=lambda p:str(p[0])))+'}'

    @property
    def args(self):
        from . import Tuple
        return tuple(Tuple(k,v) for k,v in sorted(self._data.items(),key=lambda p:str(p[0])))

    @property
    def free_symbols(self):
        return set().union(*(k.free_symbols | v.free_symbols for k,v in self._data.items()))

    def atoms(self,*types):
        from ._expr_methods import _structural_nodes
        from ._structural_queries import matcher
        selectors=[matcher(t if isinstance(t,type) else type(t),matching=False) for t in types]
        return {node for node in _structural_nodes(self)
                if (any(select(node) for select in selectors) if types else not getattr(node,'args',()))}

    def _symox_(self):
        return _core._apply('Dict',*self.args)

    def subs(self,*args,**kwargs):
        return Dict(*((k.subs(*args,**kwargs),v.subs(*args,**kwargs)) for k,v in self._data.items()))

    def xreplace(self,rule):
        return Dict(*((k.xreplace(rule),v.xreplace(rule)) for k,v in self._data.items()))

    def has(self,*patterns):
        return any(pattern is Dict or pattern == self for pattern in patterns) or any(
            key.has(*patterns) or value.has(*patterns) for key,value in self._data.items())

    def to_sympy(self):
        from sympy import Dict as ReferenceDict
        return ReferenceDict({key.to_sympy():value.to_sympy() for key,value in self._data.items()})

    _sympy_=to_sympy


def evalf(x,prec,options):
    """Low-level binary-precision evaluation returning mpmath mpf tuples."""
    import mpmath
    mp = mpmath.mp.clone()
    from . import sympify
    prec = operator.index(prec)
    if prec < 1:
        raise ValueError('precision must be positive')
    options = dict(options)
    if prec > options.get('maxprec',prec):
        raise PrecisionExhausted('requested precision exceeds maxprec')
    digits = math.ceil(prec*math.log10(2))+5
    value = sympify(x).evalf(digits,subs=options.get('subs'),chop=options.get('chop',False))
    if value.free_symbols:
        raise ValueError('numeric evaluation has unresolved symbols')
    real,imag = value.as_real_imag()
    with mp.workprec(prec):
        re = mp.mpf(str(real.evalf(digits)))
        im = mp.mpf(str(imag.evalf(digits)))
        return (None if not re else re._mpf_,None if not im else im._mpf_,
                None if not re else prec,None if not im else prec)


def install(namespace):
    Expr = _core.Expr
    if getattr(Expr,'_public_core_installed',False):
        return
    original_sympify = namespace['sympify']
    from ._sympify import make_sympify
    namespace['sympify'] = make_sympify(original_sympify)
    namespace['Derivative'] = Derivative
    Expr._args=property(lambda self:self.args)
    Expr.is_Equality=property(lambda self:self.is_Relational and self.rel_op=='==')
    operations = ('__add__','__radd__','__sub__','__rsub__','__mul__','__rmul__',
                  '__truediv__','__rtruediv__','__pow__','__rpow__')
    Expr._op_priority = 10.0
    for name in operations:
        original = getattr(Expr,name)
        def arithmetic(original,name):
            @wraps(original)
            def operation(self,other,*args):
                priority = getattr(other,'_op_priority',None)
                if priority is not None and priority > self._op_priority:
                    reflected = '__'+name[3:] if name.startswith('__r') else '__r'+name[2:]
                    handler = getattr(other,reflected,None)
                    if handler is not None:
                        return handler(self,*args)
                scalar = isinstance(other,(Expr,int,float,complex))
                if scalar and name in ('__truediv__','__rtruediv__') and _divide_error.get() and self == 0 and other == 0:
                    raise ValueError('indeterminate 0/0')
                if scalar and not _evaluation.get() and not args:
                    value = _core.sympify(other)
                    a,b = (value,self) if name.startswith('__r') else (self,value)
                    if 'sub' in name:
                        return _core._node('Add',(a,_core._node('Mul',(-1,b),False)),False)
                    if 'add' in name:
                        return _core._node('Add',(a,b),False)
                    if 'mul' in name:
                        return _core._node('Mul',(a,b),False)
                    if 'div' in name:
                        return _core._node('Mul',(a,_core._node('Pow',(b,-1),False)),False)
                    return _core._node('Pow',(a,b),False)
                try:
                    return original(self,other,*args)
                except TypeError as error:
                    if getattr(error,'_symox_unsupported_conversion',False):
                        return NotImplemented
                    raise
            return operation
        setattr(Expr,name,arithmetic(original,name))
    original_match = Expr.match
    @wraps(original_match)
    def match(self,pattern,*args,**kwargs):
        result = original_match(self,pattern,*args,**kwargs)
        if result is not None:
            for key,value in result.items():
                if key in _wild_functions and value != key:
                    counts = _wild_functions[key]
                    if not value.is_Function or counts is not None and len(value.args) not in counts:
                        return None
        return result
    Expr.match = match
    previous = Expr.nargs
    def nargs(self):
        if self in _wild_functions:
            from . import Naturals0,FiniteSet
            counts = _wild_functions[self]
            return Naturals0 if counts is None else FiniteSet(*counts)
        return previous.__get__(self,Expr)
    Expr.nargs = property(nargs)
    Expr._public_core_installed = True
