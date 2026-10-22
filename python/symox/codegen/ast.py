"""Immutable program nodes shared by the native expression code generators."""
import symox as _S


def _immutable(value):
    if isinstance(value,(list,tuple)):return tuple(_immutable(v) for v in value)
    return value


class CodegenAST:
    _fields=()
    defaults={}
    def __init__(self,*args,**kwargs):
        if len(args)>len(self._fields):raise TypeError('too many arguments for '+type(self).__name__)
        values=dict(zip(self._fields,args))
        if set(values)&set(kwargs):raise TypeError('argument supplied twice')
        values.update(kwargs)
        if set(values)-set(self._fields):raise TypeError('unknown node fields: '+str(set(values)-set(self._fields)))
        for field in self._fields:
            if field in values:value=values[field]
            elif field in self.defaults:value=self.defaults[field]
            else:raise TypeError('missing field '+field)
            value=self._coerce(field,value);object.__setattr__(self,field,_immutable(value))
    def _coerce(self,field,value):
        if field in ('body',) and isinstance(value,(tuple,list)):return CodeBlock(*value)
        if field in ('name','text','message','label','namespace') and value is not None:return str(value)
        if field in ('symbol','target','condition','counter','first','last','step','offset') and value is not None and not isinstance(value,(CodegenAST,range)):
            return _S.Symbol(value) if isinstance(value,str) else _S.sympify(value)
        return value
    def __setattr__(self,name,value):raise AttributeError('code generation nodes are immutable')
    @property
    def args(self):return tuple(getattr(self,field) for field in self._fields)
    @property
    def func(self):return type(self)
    @property
    def free_symbols(self):
        def free(v):
            if isinstance(v,tuple):return set().union(*(free(x) for x in v))
            return getattr(v,'free_symbols',set())
        return set().union(*(free(v) for v in self.args))
    def __hash__(self):return hash((type(self),self.args))
    def __eq__(self,other):return type(self) is type(other) and self.args==other.args
    def __repr__(self):return type(self).__name__+'('+', '.join(repr(v) for v in self.args)+')'
    def subs(self,*args,**kwargs):
        def replace(v):
            if isinstance(v,tuple):return tuple(replace(x) for x in v)
            return v.subs(*args,**kwargs) if hasattr(v,'subs') else v
        return type(self)(*(replace(v) for v in self.args))
    def xreplace(self,rule):
        if self in rule:return rule[self]
        def replace(v):
            if isinstance(v,tuple):return tuple(replace(x) for x in v)
            return v.xreplace(rule) if hasattr(v,'xreplace') else rule.get(v,v)
        return type(self)(*(replace(v) for v in self.args))
    def replace(self,query,value,**kwargs):
        def replace(v):
            if isinstance(v,tuple):return tuple(replace(x) for x in v)
            return v.replace(query,value,**kwargs) if hasattr(v,'replace') and not isinstance(v,str) else v
        return type(self)(*(replace(v) for v in self.args))
    def kwargs(self,exclude=(),apply=None):
        """The fields as a dictionary, SymPy's ``Token.kwargs``."""
        kwargs={k:getattr(self,k) for k in self._fields if k not in exclude}
        if apply is not None:return {k:apply(v) for k,v in kwargs.items()}
        return kwargs
    def has(self,*patterns):return any(self==p or isinstance(p,type) and isinstance(self,p) for p in patterns) or any(getattr(v,'has',lambda *p:False)(*patterns) for v in self.args)


class Node(CodegenAST):
    _fields=('attrs',);defaults={'attrs':()}
    def attr_params(self,looking_for):
        for attr in self.attrs:
            if str(attr.name)==str(looking_for):return attr.parameters
        # Like the reference, a missing attribute gives None, not ().
        return None


class _Expression:
    def __add__(self,v):return self._symox_()+v
    __radd__=__add__
    def __sub__(self,v):return self._symox_()-v
    def __rsub__(self,v):return v-self._symox_()
    def __mul__(self,v):return self._symox_()*v
    __rmul__=__mul__
    def __truediv__(self,v):return self._symox_()/v
    def __rtruediv__(self,v):return v/self._symox_()
    def __pow__(self,v):return self._symox_()**v
    def __neg__(self):return -self._symox_()


class Token(CodegenAST):pass
class NoneToken(Token):
    # SymPy's NoneToken compares equal to None.
    def __bool__(self):return False
    def __eq__(self,other):return other is None or isinstance(other,NoneToken)
    def __hash__(self):return hash(())
    def __str__(self):return 'None'
class BreakToken(Token):pass
class ContinueToken(Token):pass
none=NoneToken();break_=BreakToken();continue_=ContinueToken()


class String(Token):
    _fields=('text',)
    def __str__(self):return self.text
class QuotedString(String):pass
class Attribute(Token):
    _fields=('name','parameters');defaults={'parameters':()}
class Type(Token):
    _fields=('name',)
    def __str__(self):return str(self.name)
    @classmethod
    def from_expr(cls,expr):
        value=_S.sympify(expr)
        return integer if value.is_integer else real if value.is_real else complex_ if value.is_complex else untyped
    def cast_check(self,value,rtol=None,atol=0,precision_targets=None):
        import numpy as np
        value=_S.sympify(value);name=str(self.name)
        if name in ('integer','intc') or isinstance(self,IntBaseType):
            if value.is_integer is not True:raise ValueError('value is not an exact integer')
            result=int(value)
            if hasattr(self,'nbits'):
                low=0 if isinstance(self,UnsignedIntType) else -(1<<(int(self.nbits)-1));high=(1<<int(self.nbits))-1 if isinstance(self,UnsignedIntType) else (1<<(int(self.nbits)-1))-1
                if not low<=result<=high:raise ValueError('integer does not fit target type')
            return _S.Integer(result)
        if isinstance(self,ComplexBaseType):return _S.sympify(complex(value.evalf()))
        bits=int(getattr(self,'nbits',64));dtype=np.float32 if bits<=32 else np.float64 if bits<=64 else np.longdouble
        numeric=float(value.evalf());converted=dtype(numeric)
        if np.isfinite(numeric) and not np.isfinite(converted):raise ValueError('floating conversion overflows')
        tolerance=float(rtol) if rtol is not None else float(np.finfo(dtype).eps)
        if abs(float(converted)-numeric)>float(atol)+tolerance*abs(numeric):raise ValueError('floating conversion exceeds tolerance')
        return _S.Float(str(converted))
class IntBaseType(Type):pass
class SignedIntType(IntBaseType):_fields=('name','nbits')
class UnsignedIntType(IntBaseType):_fields=('name','nbits')
class FloatBaseType(Type):pass


class _FloatLimits:
    """The derived limits of SymPy's ``FloatType``."""
    @property
    def max_exponent(self):
        return _S.Integer(2)**(int(self.nexp) - 1)

    @property
    def min_exponent(self):
        return 3 - self.max_exponent

    @property
    def max(self):
        two = _S.Integer(2)
        return (1 - two**-(int(self.nmant) + 1))*two**self.max_exponent

    @property
    def tiny(self):
        return _S.Integer(2)**(self.min_exponent - 1)

    @property
    def eps(self):
        return _S.Integer(2)**(-int(self.nmant))

    @property
    def dig(self):
        return _S.floor(int(self.nmant)*_S.log(2)/_S.log(10))

    @property
    def decimal_dig(self):
        return _S.ceiling((int(self.nmant) + 1)*_S.log(2)/_S.log(10) + 1)

    def cast_nocheck(self, value):
        if value == _S.oo:
            return float('inf')
        elif value == -_S.oo:
            return float('-inf')
        return _S.Float(str(_S.sympify(value).evalf(self.decimal_dig)), self.decimal_dig)


class FloatType(_FloatLimits,FloatBaseType):_fields=('name','nbits','nmant','nexp')
class ComplexBaseType(Type):pass
class ComplexType(_FloatLimits,ComplexBaseType):_fields=('name','nbits','nmant','nexp')

untyped=Type('untyped');integer=IntBaseType('integer');intc=SignedIntType('intc',32)
real=FloatBaseType('real');float32=FloatType('float32',32,23,8);float64=FloatType('float64',64,52,11);float80=FloatType('float80',80,63,15)
complex_=ComplexBaseType('complex');complex64=ComplexType('complex64',64,23,8);complex128=ComplexType('complex128',128,52,11)
bool_=Type('bool');void=Type('void');int8=SignedIntType('int8',8);int16=SignedIntType('int16',16);int32=SignedIntType('int32',32);int64=SignedIntType('int64',64)
uint8=UnsignedIntType('uint8',8);uint16=UnsignedIntType('uint16',16);uint32=UnsignedIntType('uint32',32);uint64=UnsignedIntType('uint64',64)
value_const=Attribute('value_const');pointer_const=Attribute('pointer_const');restrict=Attribute('restrict');volatile=Attribute('volatile')


class Variable(Node,_Expression):
    _fields=('symbol','type','value','attrs');defaults={'type':untyped,'value':none,'attrs':()}
    def _symox_(self):return self.symbol
    @property
    def free_symbols(self):return self.symbol.free_symbols
    @classmethod
    def deduced(cls,symbol,value=None,attrs=(),**kwargs):
        if value is None:value=symbol
        return cls(symbol,kwargs.get('type',Type.from_expr(value)),value,attrs)
    def as_Declaration(self,**kwargs):
        values=dict(zip(self._fields,self.args));values.update(kwargs);return Declaration(type(self)(**values))
class Pointer(Variable):
    def __getitem__(self,indices):return Element(self.symbol,(indices,) if not isinstance(indices,tuple) else indices)
class Declaration(CodegenAST):_fields=('variable',)
class Element(CodegenAST):
    """Element of a (possibly N-dimensional) array.

    As in the reference, the symbol and offset are sympified and the indices
    and strides become Tuples (a string such as 'ijk' gives one symbol per
    character). An Element is a code node, not an arithmetic expression: it
    is printed with its strides and offset by the code printers."""
    _fields=('symbol','indices','strides','offset');defaults={'strides':none,'offset':none}
    def _coerce(self,field,value):
        if value is None or isinstance(value,NoneToken):return self.defaults.get(field,none)
        if field in ('symbol','offset'):return _S.sympify(value)
        return _S.Tuple(*value)
    def _repr_args(self,printer):
        parts=[]
        for position,(field,value) in enumerate(zip(self._fields,self.args)):
            if field in self.defaults and value==self.defaults[field]:continue
            if isinstance(value,_S.Tuple):
                joined=', '.join(printer(v) for v in value.args)
                text=('({0},)' if len(value.args)==1 else '({0})').format(joined)
            else:text=printer(value)
            parts.append(text if position==0 else field+'='+text)
        return type(self).__name__+'('+', '.join(parts)+')'
    def _sympyrepr(self,printer,*args,**kwargs):return self._repr_args(printer._print)
    def _sympystr(self,printer,*args,**kwargs):return self._repr_args(printer._print)
    def __repr__(self):
        from symox.printing.repr import srepr
        return self._repr_args(srepr)
    def __str__(self):
        from symox.printing.str import sstr
        return self._repr_args(sstr)
class FunctionCall(CodegenAST,_Expression):
    _fields=('name','function_args')
    def _symox_(self):return _S.Function(str(self.name))(*self.function_args)
class FunctionPrototype(Node):
    _fields=('return_type','name','parameters','attrs');defaults={'attrs':()}
    @classmethod
    def from_FunctionDefinition(cls,definition):return cls(definition.return_type,definition.name,definition.parameters,definition.attrs)
class FunctionDefinition(FunctionPrototype):_fields=('return_type','name','parameters','body','attrs')
class For(CodegenAST):_fields=('target','iterable','body')
class While(CodegenAST):_fields=('condition','body')
class Scope(CodegenAST):_fields=('body',)
class Return(CodegenAST):
    _fields=('return',)
    return_=property(lambda self:getattr(self,'return'))
class RuntimeError_(CodegenAST):_fields=('message',)
class Raise(CodegenAST):_fields=('exception',)
class Stream(CodegenAST):_fields=('name',)
stdout=Stream('stdout');stderr=Stream('stderr')
class Print(CodegenAST):
    _fields=('print_args','format_string','file');defaults={'format_string':none,'file':none}
def _check_assignment(lhs,rhs):
    """The matrix structure checks of SymPy's ``AssignmentBase._check_args``."""
    lhs_is_mat=hasattr(lhs,'shape') and not isinstance(lhs,Indexed)
    rhs_is_mat=hasattr(rhs,'shape') and not isinstance(rhs,Indexed)
    if lhs_is_mat:
        if not rhs_is_mat:raise ValueError("Cannot assign a scalar to a matrix.")
        elif lhs.shape!=rhs.shape:raise ValueError("Dimensions of lhs and rhs do not align.")
    elif rhs_is_mat and not lhs_is_mat:raise ValueError("Cannot assign a matrix to a scalar.")
class AssignmentBase(CodegenAST):
    _fields=('lhs','rhs')
    def __init__(self,lhs,rhs):
        lhs,rhs=_expression(lhs),_expression(rhs);_check_assignment(lhs,rhs)
        super().__init__(lhs,rhs)
class AugmentedAssignment(AssignmentBase):
    binop=''
    @property
    def op(self):return self.binop+'='
class AddAugmentedAssignment(AugmentedAssignment):binop='+'
class SubAugmentedAssignment(AugmentedAssignment):binop='-'
class MulAugmentedAssignment(AugmentedAssignment):binop='*'
class DivAugmentedAssignment(AugmentedAssignment):binop='/'
class ModAugmentedAssignment(AugmentedAssignment):binop='%'
def aug_assign(lhs,op,rhs):
    operators={'+':AddAugmentedAssignment,'-':SubAugmentedAssignment,'*':MulAugmentedAssignment,'/':DivAugmentedAssignment,'%':ModAugmentedAssignment}
    if op not in operators:raise ValueError('unknown augmented assignment operator')
    return operators[op](lhs,rhs)

from dataclasses import dataclass


def _expression(value):
    from symox import sympify
    return value if isinstance(value, (Indexed, IndexedBase, CodegenAST)) else sympify(value)


@dataclass(frozen=True, init=False)
class Assignment(AssignmentBase):
    lhs: object
    rhs: object
    op = ':='

    def __init__(self, lhs, rhs):
        lhs, rhs = _expression(lhs), _expression(rhs)
        _check_assignment(lhs, rhs)
        object.__setattr__(self, "lhs", lhs)
        object.__setattr__(self, "rhs", rhs)

    @property
    def args(self):
        return self.lhs, self.rhs

    @property
    def free_symbols(self):
        return self.lhs.free_symbols | self.rhs.free_symbols


@dataclass(frozen=True, init=False)
class CodeBlock(CodegenAST):
    _statements: tuple

    def __init__(self, *args):
        object.__setattr__(self, "_statements", tuple(args))

    @property
    def args(self):return self._statements

    @property
    def free_symbols(self):
        return set().union(*(getattr(arg, "free_symbols", set()) for arg in self.args))


@dataclass(frozen=True)
class Comment(Token):
    text: str

    def __post_init__(self):
        object.__setattr__(self, "text", str(self.text))

    def __str__(self):
        return self.text

    @property
    def args(self):
        return (self.text,)


class _IndexedArithmetic:
    """Allow indexed references in scalar symbolic expressions."""

    is_Symbol = False
    is_Number = False
    is_negative = None
    is_Pow = False
    is_Matrix = False
    is_Atom = True
    is_commutative = True

    def _symox_(self):
        from symox import Symbol
        return Symbol(str(self))

    def has(self, *patterns):
        return any(self == p or isinstance(p, type) and isinstance(self, p) for p in patterns) or any(
            getattr(a, 'has', lambda *p: False)(*patterns) for a in self.args)

    def atoms(self, *types):
        found = {self} if not types or isinstance(self, types) else set()
        for a in self.args:
            if hasattr(a, 'atoms'):
                found |= a.atoms(*types)
        return found

    def __add__(self, other): return self._symox_() + other
    def __radd__(self, other): return other + self._symox_()
    def __sub__(self, other): return self._symox_() - other
    def __rsub__(self, other): return other - self._symox_()
    def __mul__(self, other): return self._symox_() * other
    def __rmul__(self, other): return other * self._symox_()
    def __truediv__(self, other): return self._symox_() / other
    def __rtruediv__(self, other): return other / self._symox_()
    def __pow__(self, other): return self._symox_() ** other
    def __neg__(self): return -self._symox_()


@dataclass(frozen=True, init=False)
class IndexedBase(_IndexedArithmetic):
    label: object
    shape: object

    def __init__(self, label, shape=None, **assumptions):
        from symox import Symbol
        label = Symbol(label, **assumptions) if isinstance(label, str) else _expression(label)
        object.__setattr__(self, "label", label)
        object.__setattr__(self, "shape", None if shape is None else tuple(_expression(d) for d in shape))

    @property
    def args(self):
        return (self.label,) if self.shape is None else (self.label, self.shape)

    @property
    def free_symbols(self):
        return self.label.free_symbols

    def __getitem__(self, indices):
        return Indexed(self, *(indices if isinstance(indices, tuple) else (indices,)))

    def __str__(self):
        return str(self.label)


@dataclass(frozen=True, init=False)
class Indexed(_IndexedArithmetic):
    base: IndexedBase
    indices: tuple

    def __init__(self, base, *indices):
        if not indices:
            raise ValueError("Indexed needs at least one index")
        if not isinstance(base, IndexedBase):
            base = IndexedBase(base)
        if base.shape is not None and len(indices) != len(base.shape):
            raise ValueError("Rank mismatch")
        object.__setattr__(self, "base", base)
        object.__setattr__(self, "indices", tuple(_expression(i) for i in indices))

    @property
    def args(self):
        return (self.base,) + self.indices

    @property
    def free_symbols(self):
        return self.base.free_symbols | set().union(*(i.free_symbols for i in self.indices)) | {self}

    def __str__(self):
        return f"{self.base}[{', '.join(map(str, self.indices))}]"

    @property
    def rank(self):
        return len(self.indices)

    @property
    def shape(self):
        """The base shape, or else the ranges of the indices, as SymPy's
        ``Indexed.shape``."""
        if self.base.shape:
            return self.base.shape
        from symox import Tuple
        from symox.tensor.indexed import IndexException
        sizes = []
        for i in self.indices:
            upper = getattr(i, 'upper', None)
            lower = getattr(i, 'lower', None)
            if None in (upper, lower):
                raise IndexException("Range is not defined for all indices in: %s" % self)
            try:
                size = upper - lower + 1
            except TypeError:
                raise IndexException("Shape cannot be inferred from Idx with undefined range: %s" % self)
            sizes.append(size)
        return Tuple(*sizes)

    def _symox_(self):
        # Flatten multidimensional references consistently when they enter
        # native scalar arithmetic (whose atoms are scalar symbols).
        from symox import Symbol
        from symox.printing._legacy import LegacyCXX17CodePrinter
        return Symbol(LegacyCXX17CodePrinter().doprint(self))

__all__ = [name for name,value in globals().items() if not name.startswith("_") and (isinstance(value,type) and issubclass(value,CodegenAST) or name in ("aug_assign",))]
__all__ += 'untyped integer intc real float32 float64 float80 complex_ complex64 complex128 bool_ void int8 int16 int32 int64 uint8 uint16 uint32 uint64 value_const pointer_const restrict volatile none break_ continue_ stdout stderr Indexed IndexedBase'.split()
