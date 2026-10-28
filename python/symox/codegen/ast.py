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
    def has(self,*patterns):
        if any(self==p or isinstance(p,type) and isinstance(self,p) for p in patterns):return True
        # A native expression cannot contain a code-generation node, so only
        # the expression patterns are searched below it.
        plain=tuple(p for p in patterns if not isinstance(p,CodegenAST) and not (isinstance(p,type) and issubclass(p,CodegenAST)))
        for v in self.args:
            if isinstance(v,CodegenAST):
                if v.has(*patterns):return True
            elif plain and getattr(v,'has',None) is not None and v.has(*plain):return True
        return False


class Node(CodegenAST):
    _fields=('attrs',);defaults={'attrs':()}
    def attr_params(self,looking_for):
        name=str(looking_for.name) if isinstance(looking_for,Attribute) else str(looking_for)
        for attr in self.attrs:
            if str(attr.name)==name:return attr.parameters
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
    def __str__(self):return 'NoneToken()'
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
        # SymPy's ``Type.from_expr``.
        if isinstance(expr,(float,_S.Float)):return real
        if isinstance(expr,(int,_S.Integer)) or getattr(expr,'is_integer',False):return integer
        if getattr(expr,'is_real',False):return real
        if isinstance(expr,complex) or getattr(expr,'is_complex',False):return complex_
        if isinstance(expr,bool) or getattr(expr,'is_Relational',False):return bool_
        raise ValueError('Could not deduce type from expr.')
    def cast_nocheck(self,value):return value
    def _check(self,value):pass
    def cast_check(self,value,rtol=None,atol=0,precision_targets=None):
        """SymPy's ``Type.cast_check``: cast, check the bounds, then the rounding."""
        val=_S.sympify(value)
        ten=_S.Integer(10)
        exp10=getattr(self,'decimal_dig',None)
        if rtol is None:rtol=1e-15 if exp10 is None else 2.0*ten**(-exp10)
        def tol(num):return atol+rtol*abs(num)
        new_val=self.cast_nocheck(value)
        self._check(new_val)
        delta=new_val-val
        if abs(delta)>tol(val):raise ValueError("Casting gives a significantly different value.")
        return new_val
class IntBaseType(Type):
    def cast_nocheck(self,i):return _S.Integer(int(i))
class _SizedIntType(IntBaseType):
    _fields=('name','nbits')
    def _check(self,value):
        if value<self.min:raise ValueError("Value is too small: %d < %d"%(value,self.min))
        if value>self.max:raise ValueError("Value is too big: %d > %d"%(value,self.max))
class SignedIntType(_SizedIntType):
    min=property(lambda self:-2**(_S.Integer(self.nbits)-1))
    max=property(lambda self:2**(_S.Integer(self.nbits)-1)-1)
class UnsignedIntType(_SizedIntType):
    min=property(lambda self:_S.Integer(0))
    max=property(lambda self:2**_S.Integer(self.nbits)-1)
class FloatBaseType(Type):
    def cast_nocheck(self,value):return _S.Float(value)


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

    def _check(self, value):
        if value < -self.max:
            raise ValueError("Value is too small: %d < %d" % (value, -self.max))
        if value > self.max:
            raise ValueError("Value is too big: %d > %d" % (value, self.max))
        if abs(value) < self.tiny:
            raise ValueError("Smallest (absolute) value for data type bigger than new value.")


class FloatType(_FloatLimits,FloatBaseType):_fields=('name','nbits','nmant','nexp')
class ComplexBaseType(Type):
    # SymPy's ``ComplexBaseType`` derives from ``FloatBaseType``: the real
    # and imaginary parts are cast and checked as floats.
    def _part_cast(self,value):return FloatBaseType.cast_nocheck(self,value)
    def _part_check(self,value):pass
    def cast_nocheck(self,value):
        return self._part_cast(_S.re(value))+self._part_cast(_S.im(value))*_S.I
    def _check(self,value):
        self._part_check(_S.re(value));self._part_check(_S.im(value))
class ComplexType(_FloatLimits,ComplexBaseType):
    _fields=('name','nbits','nmant','nexp')
    def _part_cast(self,value):return _FloatLimits.cast_nocheck(self,value)
    def _part_check(self,value):_FloatLimits._check(self,value)
    cast_nocheck=ComplexBaseType.cast_nocheck
    _check=ComplexBaseType._check

untyped=Type('untyped');integer=IntBaseType('integer');intc=IntBaseType('intc')
real=FloatBaseType('real');float32=FloatType('float32',32,23,8);float64=FloatType('float64',64,52,11);float80=FloatType('float80',80,63,15)
complex_=ComplexBaseType('complex');complex64=ComplexType('complex64',64,23,8);complex128=ComplexType('complex128',128,52,11)
bool_=Type('bool');void=Type('void');int8=SignedIntType('int8',8);int16=SignedIntType('int16',16);int32=SignedIntType('int32',32);int64=SignedIntType('int64',64)
uint8=UnsignedIntType('uint8',8);uint16=UnsignedIntType('uint16',16);uint32=UnsignedIntType('uint32',32);uint64=UnsignedIntType('uint64',64)
value_const=Attribute('value_const');pointer_const=Attribute('pointer_const');restrict=Attribute('restrict');volatile=Attribute('volatile')


class Variable(Node,_Expression):
    _fields=('symbol','type','value','attrs');defaults={'type':untyped,'value':none,'attrs':()}
    def _coerce(self,field,value):
        # SymPy's ``_construct_value = staticmethod(sympify)``.
        if field=='value' and value is not None and not isinstance(value,(CodegenAST,_S.Basic)):
            try:return _S.sympify(value)
            except Exception:return value
        return super()._coerce(field,value)
    def _symox_(self):return self.symbol
    @property
    def free_symbols(self):return self.symbol.free_symbols
    @classmethod
    def deduced(cls,symbol,value=None,attrs=(),**kwargs):
        # SymPy's ``Variable.deduced``.
        if isinstance(symbol,Variable):return symbol
        cast_check=kwargs.get('cast_check',True)
        try:type_=Type.from_expr(symbol)
        except ValueError:type_=Type.from_expr(value)
        if value is not None and cast_check:value=type_.cast_check(value)
        return cls(symbol,type_,none if value is None else value,attrs)
    def as_Declaration(self,**kwargs):
        values=dict(zip(self._fields,self.args));values.update(kwargs);return Declaration(type(self)(**values))
class Pointer(Variable):
    def __getitem__(self,indices):return Element(self.symbol,(indices,) if not isinstance(indices,tuple) else indices)
class Declaration(CodegenAST):
    _fields=('variable',)
    def _coerce(self,field,value):
        # SymPy's ``_construct_variable = Variable``.
        return value if isinstance(value,Variable) else Variable(value)
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
    # SymPy's FunctionCall has is_commutative None, which Mul treats as
    # noncommutative: keep its factor order and Pow(-1) printing.
    def _symox_(self):return _S.Function(str(self.name),commutative=False)(*self.function_args)
def _function_parameter(arg):
    # SymPy's ``FunctionPrototype._construct_parameters``.
    if isinstance(arg,Declaration):return arg.variable
    if isinstance(arg,Variable):return arg
    return Variable.deduced(arg)
class FunctionPrototype(Node):
    _fields=('return_type','name','parameters','attrs');defaults={'attrs':()}
    def _coerce(self,field,value):
        if field=='parameters':return tuple(_function_parameter(a) for a in value)
        return super()._coerce(field,value)
    @classmethod
    def from_FunctionDefinition(cls,func_def):
        if not isinstance(func_def,FunctionDefinition):raise TypeError("func_def is not an instance of FunctionDefinition")
        return cls(**func_def.kwargs(exclude=('body',)))
class FunctionDefinition(FunctionPrototype):
    _fields=('return_type','name','parameters','body','attrs')
    def _coerce(self,field,value):
        if field=='body' and not isinstance(value,CodeBlock):return CodeBlock(*value)
        return super()._coerce(field,value)
    @classmethod
    def from_FunctionPrototype(cls,func_proto,body):
        if not isinstance(func_proto,FunctionPrototype):raise TypeError("func_proto is not an instance of FunctionPrototype")
        return cls(body=body,**func_proto.kwargs())
def _code_block(body):
    return body if isinstance(body,CodeBlock) else CodeBlock(*body)
class For(CodegenAST):
    _fields=('target','iterable','body')
    def _coerce(self,field,value):
        # SymPy's For._construct_iterable and _construct_body.
        if field=='iterable':
            if isinstance(value,range):return _S.Range(value.start,value.stop,value.step)
            if isinstance(value,list):value=tuple(value)
            if not hasattr(value,'__iter__'):raise TypeError("iterable must be an iterable")
            return value if isinstance(value,(CodegenAST,_S.Basic,_S.Set)) else _S.sympify(value)
        if field=='body':return _code_block(value)
        return super()._coerce(field,value)
class While(CodegenAST):
    _fields=('condition','body')
    def _coerce(self,field,value):
        if field=='condition':return _S.sympify(value)
        if field=='body':return _code_block(value)
        return super()._coerce(field,value)
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
    def _coerce(self,field,value):
        if field=='print_args':return tuple(a if isinstance(a,CodegenAST) else _S.sympify(a) for a in value)
        if field=='format_string' and isinstance(value,str):return QuotedString(value)
        return super()._coerce(field,value)
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

    def __repr__(self):
        return '%s(%s, %s)' % (type(self).__name__, self.lhs, self.rhs)


@dataclass(frozen=True, init=False)
class CodeBlock(CodegenAST):
    _statements: tuple

    def __init__(self, *args):
        object.__setattr__(self, "_statements", tuple(args))

    @property
    def args(self):return self._statements
    def __iter__(self):return iter(self.args)

    def __repr__(self):
        return 'CodeBlock(\n%s\n)' % ',\n'.join(repr(arg) for arg in self.args)

    @property
    def left_hand_sides(self):
        return tuple(a.lhs for a in self.args if isinstance(a, AssignmentBase))

    @property
    def right_hand_sides(self):
        return tuple(a.rhs for a in self.args if isinstance(a, AssignmentBase))

    def atoms(self, *types):
        found = set()
        for arg in self.args:
            for part in getattr(arg, 'args', ()):
                if hasattr(part, 'atoms'):
                    found |= part.atoms(*types)
        return found

    @classmethod
    def topological_sort(cls, assignments):
        """SymPy's ``CodeBlock.topological_sort``."""
        from collections import defaultdict
        from symox.utilities.iterables import topological_sort
        if not all(isinstance(i, Assignment) for i in assignments):
            raise NotImplementedError("CodeBlock.topological_sort only supports Assignments")
        if any(isinstance(i, AugmentedAssignment) for i in assignments):
            raise NotImplementedError("CodeBlock.topological_sort does not yet work with AugmentedAssignments")
        A = list(enumerate(assignments))
        var_map = defaultdict(list)
        for node in A:
            var_map[node[1].lhs].append(node)
        E = []
        for dst_node in A:
            for s in dst_node[1].rhs.free_symbols:
                for src_node in var_map[s]:
                    E.append((src_node, dst_node))
        ordered_assignments = topological_sort([A, E])
        return cls(*[a for i, a in ordered_assignments])

    def cse(self, symbols=None, optimizations=None, postprocess=None, order='canonical'):
        """SymPy's ``CodeBlock.cse``."""
        from symox.simplify.cse_main import cse
        from symox.utilities.iterables import numbered_symbols, filter_symbols
        if not all(isinstance(i, Assignment) for i in self.args):
            raise NotImplementedError("CodeBlock.cse only supports Assignments")
        if any(isinstance(i, AugmentedAssignment) for i in self.args):
            raise NotImplementedError("CodeBlock.cse does not yet work with AugmentedAssignments")
        for i, lhs in enumerate(self.left_hand_sides):
            if lhs in self.left_hand_sides[:i]:
                raise NotImplementedError("Duplicate assignments to the same "
                    "variable are not yet supported (%s)" % lhs)
        existing_symbols = self.atoms(_S.Symbol)
        if symbols is None:
            symbols = numbered_symbols()
        symbols = filter_symbols(symbols, existing_symbols)
        replacements, reduced_exprs = cse(list(self.right_hand_sides), symbols=symbols,
            optimizations=optimizations, postprocess=postprocess, order=order)
        new_block = [Assignment(var, expr) for var, expr in zip(self.left_hand_sides, reduced_exprs)]
        new_assignments = [Assignment(var, expr) for var, expr in replacements]
        return self.topological_sort(new_assignments + new_block)

    @property
    def free_symbols(self):
        return set().union(*(getattr(arg, "free_symbols", set()) for arg in self.args))


@dataclass(frozen=True)
class Comment(Token):
    text: str
    _fields=('text',)

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

    # Reference printing: IndexedBase('A') is shown as A, not a dataclass.
    __repr__ = __str__

    def _latex(self, printer):
        return printer._print(self.label)


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

    __repr__ = __str__

    def _latex(self, printer):
        # Reference LatexPrinter._print_Indexed.
        return '{%s}_{%s}' % (printer._print(self.base), ','.join(printer._print(i) for i in self.indices))

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

    def __hash__(self):
        # The hash of the symbol it becomes in arithmetic, so that results
        # keyed by that symbol (``solve``) are found with ``A[1]`` itself.
        try:
            return self.__dict__['_symbol_hash']
        except KeyError:
            value = hash(self._symox_())
            object.__setattr__(self, '_symbol_hash', value)
            return value

    def _sympy_(self):
        # Explicit interoperability only; native arithmetic never imports SymPy.
        import sympy
        from symox import to_sympy
        from symox._public_matrix import Idx
        shape = None if self.base.shape is None else tuple(map(to_sympy, self.base.shape))
        base = sympy.IndexedBase(to_sympy(self.base.label), shape=shape)
        indices = []
        for index in self.indices:
            if isinstance(index, Idx):
                bounds = None if index.lower is None else (to_sympy(index.lower), to_sympy(index.upper))
                indices.append(sympy.Idx(to_sympy(index.label), bounds))
            else:
                indices.append(to_sympy(index))
        return base[tuple(indices)]

    def _symox_(self):
        # Flatten multidimensional references consistently when they enter
        # native scalar arithmetic (whose atoms are scalar symbols).
        from symox.printing._legacy import LegacyCXX17CodePrinter
        symbol = _indexed_symbol(self, LegacyCXX17CodePrinter().doprint(self))
        _INDEXED_SYMBOLS.setdefault(symbol, self)
        return symbol


# The reference classes live in sympy.tensor.indexed.
IndexedBase.__module__ = Indexed.__module__ = 'symox.tensor.indexed'


# The scalar symbols that Indexed references became in native arithmetic,
# so that ``expr.atoms(Indexed)`` finds them (SymPy's
# ``Tuple(A[1] + A[2] - 3).atoms(Indexed)`` is ``{A[1], A[2]}``).
_INDEXED_SYMBOLS = {}


def _indexed_symbol(indexed, name):
    """Preserve structural identity when an Indexed enters scalar arithmetic."""
    from symox import Tuple, S, _core
    from symox._public_matrix import Idx
    entries = []
    for index in indexed.indices:
        if isinstance(index, Idx):
            fields = (S.true, index.label)
            if index.lower is not None:
                fields += (index.lower, index.upper)
        else:
            fields = (S.false, index)
        entries.append(Tuple(*fields))
    shape = S.false if indexed.base.shape is None else Tuple(*indexed.base.shape)
    return _core._indexed_symbol(name, Tuple(indexed.base.label, shape, Tuple(*entries)))


def _indexed_from_sympy(indexed):
    """Adapt an explicitly supplied foreign Indexed without losing bounds."""
    from symox import sympify
    from symox._public_matrix import Idx
    shape = indexed.base.shape
    base = IndexedBase(sympify(indexed.base.label), None if shape is None else tuple(map(sympify, shape)))
    indices = []
    for index in indexed.indices:
        if type(index).__name__ == 'Idx':
            bounds = None if index.lower is None else (sympify(index.lower), sympify(index.upper))
            indices.append(Idx(sympify(index.label), bounds))
        else:
            indices.append(sympify(index))
    return Indexed(base, *indices)._symox_()


def _indexed_reference(symbol):
    """Recover metadata after native arithmetic or serialization."""
    try:
        indexed = _INDEXED_SYMBOLS.get(symbol)
    except TypeError:
        return None
    if indexed is not None:
        return indexed
    identity = getattr(symbol, '_scalar_identity', None)
    if identity is None:
        return None
    from symox import S
    from symox._public_matrix import Idx
    label, shape, entries = identity.args
    base = IndexedBase(label, None if shape == S.false else shape.args)
    indices = []
    for entry in entries.args:
        fields = entry.args
        if fields[0] == S.true:
            indices.append(Idx(fields[1], None if len(fields) == 2 else fields[2:]))
        else:
            indices.append(fields[1])
    indexed = Indexed(base, *indices)
    _INDEXED_SYMBOLS[symbol] = indexed
    return indexed



def _install_indexed_assumptions():
    """IndexedBase keeps the assumptions of its label (given as keywords or
    inherited from a Symbol label) and Indexed those of its base, as the
    reference ``IndexedBase._set_assumptions``."""
    # sympy.core.assumptions._assume_defined (the facts module cannot be
    # imported while the package initializes).
    facts = ('algebraic', 'antihermitian', 'commutative', 'complex', 'composite',
             'even', 'extended_negative', 'extended_nonnegative',
             'extended_nonpositive', 'extended_nonzero', 'extended_positive',
             'extended_real', 'finite', 'hermitian', 'imaginary', 'infinite',
             'integer', 'irrational', 'negative', 'noninteger', 'nonnegative',
             'nonpositive', 'nonzero', 'odd', 'polar', 'positive', 'prime', 'rational',
             'real', 'transcendental', 'zero',)

    def base_fact(name):
        return property(lambda self: getattr(self.label, name))

    def indexed_fact(name):
        return property(lambda self: getattr(self.base.label, name))

    for fact in facts:
        name = 'is_' + fact
        setattr(IndexedBase, name, base_fact(name))
        setattr(Indexed, name, indexed_fact(name))


_install_indexed_assumptions()


def _install_indexed_atoms():
    """``expr.atoms(Indexed)`` lists the symbols that Indexed references
    became in native arithmetic (the reference finds its Indexed nodes)."""
    from symox import _core
    native_atoms = _core.Expr.atoms

    def atoms(self, *types):
        found = native_atoms(self, *types)
        if any(t is Indexed for t in types):
            found |= {s for s in native_atoms(self, _S.Symbol) if _indexed_reference(s) is not None}
        return found

    atoms.__doc__ = native_atoms.__doc__
    setattr(_core.Expr, 'atoms', atoms)


_install_indexed_atoms()



# The reference Token._sympyrepr/_sympystr layout: the first field is
# positional, later ones are keywords, defaults are omitted and a body is
# indented one level deeper than its owner.
_STRING_FIELDS=('name','message','local','original')


def _token_value(value,repr_mode,il,joiner=', '):
    if isinstance(value,CodeBlock):
        joined=(',\n'+' '*il).join(_token_value(a,repr_mode,il) for a in value.args)
        return ' '*(il-4)+'CodeBlock(\n'+' '*il+joined+'\n'+' '*(il-4)+')'
    if isinstance(value,AssignmentBase) and repr_mode:
        from symox.printing.repr import srepr
        return '%s(%s, %s)'%(type(value).__name__,srepr(value.lhs),srepr(value.rhs))
    if isinstance(value,_TOKEN_CLASSES):return _token_text(value,repr_mode,il)
    if isinstance(value,bool):return repr(value)
    if isinstance(value,int):return 'Integer(%d)'%value if repr_mode else str(value)
    if isinstance(value,str):return repr(value) if repr_mode else value
    if repr_mode:
        from symox.printing.repr import srepr
        return srepr(value)
    if isinstance(value,_S.Float):
        # A nested Float prints stripped (the reference full_prec='auto').
        from symox.printing.str import sstr
        return sstr(value,full_prec=False)
    return str(value)


def _token_text(obj,repr_mode,il=0):
    if isinstance(obj,NoneToken):return 'NoneToken()'
    if not repr_mode:
        if isinstance(obj,String):return obj.text
        if isinstance(obj,Type):return str(obj.name)
        if isinstance(obj,Attribute):
            result=str(obj.name)
            if obj.parameters:result+='(%s)'%', '.join(_token_value(a,False,0) for a in obj.parameters)
            return result
    parts=[]
    for i,field in enumerate(obj._fields):
        value=getattr(obj,field)
        if field in obj.defaults and value==obj.defaults[field]:continue
        ilvl=il+4 if field=='body' else 0
        joiner=',\n'+' '*ilvl if field=='body' else ', '
        if isinstance(value,str) and field in _STRING_FIELDS and not isinstance(obj,String):
            text="String(%r)"%value if repr_mode else value
        elif isinstance(value,(tuple,_S.Tuple)) and not isinstance(value,CodegenAST):
            items=value.args if isinstance(value,_S.Tuple) else value
            joined=joiner.join(_token_value(a,repr_mode,ilvl) for a in items)
            if field=='body':text='(\n'+' '*ilvl+joined+',\n'+' '*(ilvl-4)+')'
            else:text=('({0},)' if len(items)==1 else '({0})').format(joined)
        else:text=_token_value(value,repr_mode,ilvl)
        parts.append(text.lstrip() if i==0 else '%s=%s'%(field,text.lstrip()))
    return '%s(%s)'%(type(obj).__name__,', '.join(parts))


_TOKEN_CLASSES=(Token,Variable,Declaration,FunctionCall,FunctionPrototype,For,While,Scope,Return,RuntimeError_,Raise,Stream,Print)
for _cls in _TOKEN_CLASSES:
    _cls.__repr__=lambda self:_token_text(self,True)
    _cls._sympyrepr=lambda self,printer,*args,**kwargs:_token_text(self,True)
    _cls._sympystr=lambda self,printer,*args,**kwargs:_token_text(self,False)
for _cls in (Token,Variable,Declaration,FunctionCall,FunctionPrototype,For,While,Scope,Return,RuntimeError_,Raise,Stream,Print):
    if '__str__' not in _cls.__dict__:_cls.__str__=lambda self:_token_text(self,False)
Comment.__repr__=lambda self:_token_text(self,True)
del _cls


__all__ = [name for name,value in globals().items() if not name.startswith("_") and (isinstance(value,type) and issubclass(value,CodegenAST) or name in ("aug_assign",))]
__all__ += 'untyped integer intc real float32 float64 float80 complex_ complex64 complex128 bool_ void int8 int16 int32 int64 uint8 uint16 uint32 uint64 value_const pointer_const restrict volatile none break_ continue_ stdout stderr Indexed IndexedBase'.split()
