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
# Equal WildFunctions may carry different nargs, as in the reference: each
# created object keeps its own counts, keyed by identity.
_wild_function_instances = {}


def _wild_nargs(node):
    entry = _wild_function_instances.get(id(node))
    if entry is not None and entry[0] is node:
        return entry[1]
    return _wild_functions[node]


class SympifyError(ValueError,TypeError):
    def __init__(self,expr,base_exc=None):
        self.expr,self.base_exc = expr,base_exc
        super().__init__('cannot sympify %r%s' % (expr,': '+str(base_exc) if base_exc else ''))

    def __str__(self):
        # SymPy's SympifyError text.
        if self.base_exc is None:
            return "SympifyError: %r" % (self.expr,)
        return ("Sympify of expression '%s' failed, because of exception being "
                "raised:\n%s: %s" % (self.expr, self.base_exc.__class__.__name__,
                                     str(self.base_exc)))


class PrecisionExhausted(ArithmeticError):
    """A requested numerical accuracy could not be attained."""


class Derivative(metaclass=_ExprType):
    """Held native derivative, with the reference constructor keyword schema."""
    _matches=staticmethod(lambda expr:expr.is_Derivative)
    # ArrayDerivative (and the held matrix Derivative) subclass Derivative.
    _foreign_instance=staticmethod(lambda value:getattr(type(value),'_symox_held_matrix_derivative',False))

    def __new__(cls,expr,*variables,**kwargs):
        evaluate=kwargs.get('evaluate',False)
        if evaluate and not kwargs.get('simplify',True):
            from ._differentiation import unsimplified_diff
            return unsimplified_diff(expr,variables)
        if isinstance(expr,_core._MatrixExpr) and variables:
            # Derivative(X, v) of a matrix expression stays held, as in SymPy
            # (doit: ZeroMatrix(*X.shape) when X does not contain v, else
            # X.diff(v)).
            from .tensor.array import array_derivatives
            held=array_derivatives.held_matrix_derivative(expr,variables,array_derivatives.Derivative,contained=True)
            if held is not None:
                return held.doit() if evaluate else held
        result=_core.Derivative(expr,*variables,evaluate=bool(evaluate))
        if evaluate and not getattr(result,'is_Derivative',False):
            from ._differentiation import _finish_multiple
            return _finish_multiple(result,_core.sympify(expr),variables)
        return result


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


def _wild_function_value(expr):
    # Marks a WildFunction: a distinct object from the Wild of the same name.
    # WildFunction.matches accepts any Function, a WildFunction included.
    return bool(getattr(expr,'is_Function',False)) or _is_wild_function_node(expr)


def _is_wild_function_node(expr):
    try:
        return expr in _wild_functions
    except TypeError:
        return False


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
        node = _core.Wild(name,(),(_wild_function_value,),**assumptions)
        _wild_functions[node] = counts
        _wild_function_instances[id(node)] = (node,counts)
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


def _modular_pow(self, other, mod):
    """``Expr.__pow__(other, mod)`` of the reference: the integer modular
    power (with ``mod_inverse`` for a negative exponent), else
    ``self**other % mod``."""
    from . import sympify, mod_inverse
    from .ntheory._intfunc import as_int
    try:
        _self, other, mod = as_int(self), as_int(other), as_int(mod)
        if other >= 0:
            return sympify(pow(_self, other, mod))
        return sympify(mod_inverse(pow(_self, -other, mod), mod))
    except ValueError:
        power = self**other
        try:
            return power % mod
        except TypeError:
            return NotImplemented


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
            operand_dispatch=getattr(original,'_symox_operand_dispatch',None)
            if operand_dispatch is None:
                tensor_type,tensor_operation,normal_operation=None,None,original
            else:
                tensor_type,tensor_operation,normal_operation=operand_dispatch
            reflected = '__'+name[3:] if name.startswith('__r') else '__r'+name[2:]
            division = name in ('__truediv__','__rtruediv__')
            scalars = (Expr,int,float,complex)
            @wraps(original)
            def operation(self,other,*args):
                if name == '__pow__' and args and args[0] is not None:
                    return _modular_pow(self,other,args[0])
                # Exact built-in numeric types have no symbolic priority or
                # foreign conversion hooks. Subclasses keep full dispatch.
                other_type = type(other)
                plain_number = other_type is int or other_type is float or other_type is complex
                priority = None if plain_number else getattr(other,'_op_priority',None)
                if priority is not None and priority > self._op_priority:
                    handler = getattr(other,reflected,None)
                    if handler is not None:
                        return handler(self,*args)
                scalar = plain_number or isinstance(other,scalars)
                if scalar and division and _divide_error.get() and self == 0 and other == 0:
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
                    if tensor_operation is not None and isinstance(other,tensor_type):
                        return tensor_operation(self,other,*args)
                    return normal_operation(self,other,*args)
                except TypeError as error:
                    if getattr(error,'_symox_unsupported_conversion',False):
                        return NotImplemented
                    raise
            # This wrapper owns priority/context logic and cannot advertise a
            # bypass of that logic to another consumer of dispatch metadata.
            operation.__dict__.pop('_symox_operand_dispatch',None)
            return operation
        setattr(Expr,name,arithmetic(original,name))
    original_match = Expr.match
    @wraps(original_match)
    def match(self,pattern,*args,**kwargs):
        result = original_match(self,pattern,*args,**kwargs)
        if result is None and _is_wild_function_node(pattern) and _is_wild_function_node(self):
            # WildFunction.matches accepts any Function, a WildFunction
            # included, whose single argument is its name in the reference.
            counts = _wild_nargs(pattern)
            return {pattern:self} if counts is None or 1 in counts else None
        if result is not None:
            for key,value in result.items():
                if key in _wild_functions:
                    counts = _wild_nargs(key)
                    function = value.is_Function or _is_wild_function_node(value)
                    if not function or counts is not None and len(value.args) not in counts:
                        return None
        return result
    Expr.match = match
    original_matches = Expr.matches
    @wraps(original_matches)
    def matches(self,expr,repl_dict=None,old=False):
        if _is_wild_function_node(self) and _is_wild_function_node(expr):
            counts = _wild_nargs(self)
            if counts is not None and 1 not in counts:
                return None
            repl_dict = {} if repl_dict is None else dict(repl_dict)
            repl_dict[self] = expr
            return repl_dict
        return original_matches(self,expr,repl_dict,old)
    Expr.matches = matches
    previous = Expr.nargs
    def nargs(self):
        if self in _wild_functions:
            from . import Naturals0,FiniteSet
            counts = _wild_nargs(self)
            return Naturals0 if counts is None else FiniteSet(*counts)
        return previous.__get__(self,Expr)
    Expr.nargs = property(nargs)
    from ._integer_head_cache import install as install_integer_head_cache
    install_integer_head_cache(namespace,_evaluation)
    Expr._public_core_installed = True
