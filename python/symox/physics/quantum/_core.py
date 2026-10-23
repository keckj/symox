"""Ordered quantum operators, states, duality and finite representations."""
from itertools import product
import symox as S
from .._algebra import Atom,OperatorExpr,PowerAtom,_terms

class HBar(S.NumberSymbol):
    """Reduced Planck constant: a positive real number symbol, as in the reference."""
    # A module name outside symox makes this a host number-symbol class.
    __module__='symox_quantum_constants'
    _symox_user_hooks=True
    is_real=True
    is_positive=True
    is_negative=False
    is_irrational=True
    __slots__=()
    def _eval_evalf(self,prec):
        import mpmath
        # Reference HBar._as_mpf_val: the double 1.05457162e-34 at the precision.
        return S.Float(mpmath.mpf(1.05457162e-34),mpmath.libmp.prec_to_dps(prec))
    def _sympystr(self,printer,*args):return 'hbar'
    def _latex(self,printer,*args):return r'\hbar'
    def _pretty(self,printer,*args):
        from symox.printing.pretty.stringpict import prettyForm
        if printer._use_unicode:return prettyForm('\N{PLANCK CONSTANT OVER TWO PI}')
        return prettyForm('hbar')
HBar.__module__='symox.physics.quantum.constants'
hbar=HBar()


class HilbertSpaceError(Exception):pass


class HilbertSpace:
    dimension=S.oo
    def __mul__(self,other):return TensorProductHilbertSpace(self,other)
    def __add__(self,other):return DirectSumHilbertSpace(self,other)
    def __pow__(self,power):return TensorPowerHilbertSpace(self,power)
    def __repr__(self):return 'H'
    def __eq__(self,other):return type(self) is type(other) and getattr(self,'args',())==getattr(other,'args',())
    def __hash__(self):return hash((type(self),getattr(self,'args',())))


class ComplexSpace(HilbertSpace):
    def __init__(self,dimension):
        self.dimension=S.sympify(dimension)
        if self.dimension.is_integer is False or self.dimension.is_nonpositive is True:raise HilbertSpaceError('complex-space dimension must be a positive integer')
        self.args=(self.dimension,)
    def __repr__(self):return 'C(%s)'%self.dimension


class L2(HilbertSpace):
    def __init__(self,interval):self.interval=interval;self.args=(interval,)
    def __repr__(self):return 'L2(%s)'%self.interval


class FockSpace(HilbertSpace):
    def __repr__(self):return 'F'


def _space_repr(space,compound):
    text=repr(space)
    return '(%s)'%text if isinstance(space,compound) else text


class TensorProductHilbertSpace(HilbertSpace):
    def __new__(cls,*spaces):
        # Reference TensorProductHilbertSpace.eval: flatten, then merge equal
        # neighbours into tensor powers.
        if not spaces or any(not isinstance(v,HilbertSpace) for v in spaces):raise TypeError('tensor factors must be Hilbert spaces')
        flat=[];recall=False
        for space in spaces:
            if type(space) is TensorProductHilbertSpace:flat.extend(space.args);recall=True
            else:flat.append(space)
        combined=[];previous=None
        for space in flat:
            if previous is None:previous=space;continue
            if isinstance(space,TensorPowerHilbertSpace) and isinstance(previous,TensorPowerHilbertSpace) and space.base==previous.base:previous=space.base**(space.exp+previous.exp)
            elif isinstance(space,TensorPowerHilbertSpace) and space.base==previous:previous=previous**(space.exp+1)
            elif isinstance(previous,TensorPowerHilbertSpace) and space==previous.base:previous=space**(previous.exp+1)
            elif space==previous:previous=space**2
            else:combined.append(previous);previous=space
        combined.append(previous)
        if recall:return TensorProductHilbertSpace(*combined)
        if len(combined)==1:return combined[0]
        obj=object.__new__(cls);obj.args=tuple(combined);obj.dimension=S.Mul(*(v.dimension for v in combined))
        return obj
    @property
    def spaces(self):return self.args
    def __repr__(self):return '*'.join(_space_repr(v,(DirectSumHilbertSpace,TensorProductHilbertSpace)) for v in self.args)


class DirectSumHilbertSpace(HilbertSpace):
    def __new__(cls,*spaces):
        if not spaces or any(not isinstance(v,HilbertSpace) for v in spaces):raise TypeError('direct-sum terms must be Hilbert spaces')
        flat=[]
        for space in spaces:flat.extend(space.args if type(space) is DirectSumHilbertSpace else (space,))
        obj=object.__new__(cls);obj.args=tuple(flat);obj.dimension=S.Add(*(v.dimension for v in flat))
        return obj
    @property
    def spaces(self):return self.args
    def __repr__(self):return '+'.join(_space_repr(v,(TensorProductHilbertSpace,)) for v in self.args)


class TensorPowerHilbertSpace(HilbertSpace):
    def __new__(cls,space,power):
        power=S.sympify(power)
        if power==1:return space
        obj=object.__new__(cls);obj.args=(space,power);obj.dimension=space.dimension**power
        return obj
    @property
    def base(self):return self.args[0]
    @property
    def exp(self):return self.args[1]
    def __repr__(self):return '%s**%s'%(_space_repr(self.base,(DirectSumHilbertSpace,TensorProductHilbertSpace)),self.exp)


class Operator(Atom):
    is_hermitian=None
    is_unitary=None
    def __init__(self,*labels,**kwargs):
        if not labels:labels=('A',)
        super().__init__(*(S.Symbol(v) if isinstance(v,str) else S.sympify(v) for v in labels))
        self.hilbert_space=HilbertSpace()
    def _atom_repr(self):return str(self.labels[0]) if len(self.labels)==1 else 'Operator('+', '.join(map(str,self.labels))+')'
    @property
    def label(self):return tuple(self.labels)
    def adjoint(self):return Dagger(self)
    def inverse(self):return self._eval_inverse()
    inv=inverse
    def _eval_inverse(self):return self**-1
    def _apply_to(self,state,**options):return None


class HermitianOperator(Operator):
    is_hermitian=True
    def _eval_inverse(self):return self if self.is_unitary else Operator._eval_inverse(self)
class UnitaryOperator(Operator):
    is_unitary=True
    def _multiply_atom(self,other):
        # Reference UnitaryOperator: U*Dagger(U) is the identity.
        if isinstance(other,DaggerAtom) and other.argument==self:return S.Integer(1),()
        return None


class IdentityOperator(Operator):
    is_identity=True
    is_hermitian=True
    is_unitary=True
    def __init__(self,N=S.oo):self.N=S.sympify(N);super().__init__('I')
    def _key(self):return type(self),self.N
    def _apply_to(self,state,**options):return state
    def _represent(self,**options):return S.eye(int(self.N if self.N!=S.oo else options.get('ndim',2)))


class StateBase(Atom):
    # Reference QExpr._label_separator: labels print without separator.
    _label_separator=''
    is_ket=False
    is_bra=False
    is_orthogonal=False
    def __init__(self,*labels,**kwargs):
        super().__init__(*(S.Symbol(v) if isinstance(v,str) else S.sympify(v) for v in labels));self.hilbert_space=HilbertSpace()
    @property
    def label(self):return self.labels
    @classmethod
    def dual_class(cls):
        # Reference Ket/Bra.dual_class: the matching Bra/Ket class.
        import sys
        name=cls.__name__
        swapped=name.replace('Ket','\0').replace('Bra','Ket').replace('\0','Bra')
        if swapped!=name:
            target=getattr(sys.modules.get(cls.__module__),swapped,None)
            if isinstance(target,type):return target
        for klass in cls.__mro__:
            if klass in _DUALS:return _DUALS[klass]
        raise NotImplementedError('dual_class must be implemented in a subclass')
    def _inner_with(self,ket):
        if self.is_orthogonal and ket.is_orthogonal:return S.Mul(*(S.KroneckerDelta(a,b) for a,b in zip(self.labels,ket.labels))) if len(self.labels)==len(ket.labels) else S.Integer(0)
        return InnerProduct(self,ket)
    def _multiply_atom(self,other):
        if self.is_bra and isinstance(other,StateBase) and other.is_ket:return self._inner_with(other),()
        if self.is_bra and isinstance(other,OuterProduct):return self._inner_with(other.ket),(other.bra,)
        if self.is_ket and isinstance(other,StateBase) and other.is_bra:return S.Integer(1),(OuterProduct(self,other),)
        return None
    @property
    def dual(self):return Dagger(self)


class State(StateBase):pass
class KetBase(StateBase):
    is_ket=True
    def _atom_repr(self):return '|'+self._label_separator.join(map(str,self.labels))+'>'
class BraBase(StateBase):
    is_bra=True
    def _atom_repr(self):return '<'+self._label_separator.join(map(str,self.labels))+'|'
class Ket(KetBase):pass
class Bra(BraBase):pass
class _TimeDependence:
    # The last argument is the time; the label is everything before it.
    @property
    def time(self):return self.labels[-1]
    @property
    def label(self):return self.labels[:-1]
    def _time_repr(self):return self._label_separator.join(map(str,self.label))+';'+str(self.time)
class TimeDepState(_TimeDependence,State):pass
class TimeDepKet(_TimeDependence,Ket):
    def _atom_repr(self):return '|'+self._time_repr()+'>'
class TimeDepBra(_TimeDependence,Bra):
    def _atom_repr(self):return '<'+self._time_repr()+'|'
class OrthogonalState(State):is_orthogonal=True
class OrthogonalKet(Ket):is_orthogonal=True
class OrthogonalBra(Bra):is_orthogonal=True


_DUALS={Ket:Bra,Bra:Ket,KetBase:BraBase,BraBase:KetBase,TimeDepKet:TimeDepBra,TimeDepBra:TimeDepKet,OrthogonalKet:OrthogonalBra,OrthogonalBra:OrthogonalKet,State:State}


class DaggerAtom(Atom):
    def __init__(self,argument):self.argument=argument;super().__init__(argument)
    def _atom_repr(self):return 'Dagger(%s)'%self.argument
    def _multiply_atom(self,other):
        if getattr(self.argument,'is_unitary',False) and other==self.argument:return S.Integer(1),()
        return None


def _conjugate_scalar(value):
    """``conjugate`` of a scalar, with inner products ``<b|k>`` turned into ``<k|b>``."""
    result=S.conjugate(value)
    try:products=value.atoms(InnerProduct) if isinstance(value,S.Expr) else ()
    except Exception:products=()
    for product in products:result=result.xreplace({S.conjugate(product):product._adjoint_product()})
    return result


# Reference sympy.physics.quantum.dagger: importing it renames ``adjoint``.
from ... import _nc_adjoint as _nc_adjoint_module
_nc_adjoint_module.adjoint.__name__='Dagger'


def Dagger(expr):
    if isinstance(expr,S.Matrix):return expr.conjugate().T
    if isinstance(expr,InnerProduct):return expr._adjoint_product()
    if isinstance(expr,DaggerAtom):return expr.argument
    if isinstance(expr,StateBase):
        target=_DUALS.get(type(expr))
        if target is not None:return target(*expr.labels)
        if hasattr(expr,'_dual'):return expr._dual()
        return DaggerAtom(expr)
    if isinstance(expr,OuterProduct):return OuterProduct(Dagger(expr.bra),Dagger(expr.ket))
    if isinstance(expr,TensorProductAtom):return TensorProduct(*(Dagger(v) for v in expr.factors))
    if isinstance(expr,Atom) and hasattr(expr,'_dual'):return expr._dual()
    if isinstance(expr,Atom):return expr if getattr(expr,'is_hermitian',False) else DaggerAtom(expr)
    if isinstance(expr,OperatorExpr):
        result=S.Integer(0)
        for word,coefficient in expr._terms.items():
            term=_conjugate_scalar(coefficient)
            for factor in reversed(word):term=term*Dagger(factor)
            result+=term
        return result
    if isinstance(expr,S.Expr) and expr.is_commutative is False:
        # Reference Dagger: arg.adjoint(), which reverses products and keeps
        # an unevaluated Dagger(A) for a noncommutative symbol.
        return S.adjoint(expr)
    return S.conjugate(expr)


# Python atoms of the native nodes standing for them inside scalar expressions.
_NATIVE_ATOMS={}


def _register_native(atom):
    native=atom._native_atom()
    try:_NATIVE_ATOMS.setdefault(native,atom)
    except TypeError:pass
    return native


def _python_atom(native):
    try:return _NATIVE_ATOMS.get(native,native)
    except TypeError:return native


class InnerProduct(S.Function):
    """The unevaluated inner product ``<b|k>`` of a bra and a ket (a scalar)."""
    # A host expression class with its own printer, like AccumulationBounds.
    _symox_user_hooks=True
    is_commutative=True
    def __new__(cls,*args,**options):
        if len(args)!=2:raise TypeError('InnerProduct expects a bra and a ket')
        bra,ket=args
        if isinstance(bra,Atom) or isinstance(ket,Atom):
            if not isinstance(bra,BraBase):raise TypeError('BraBase subclass expected, got: %r'%(bra,))
            if not isinstance(ket,KetBase):raise TypeError('KetBase subclass expected, got: %r'%(ket,))
            args=(_register_native(bra),_register_native(ket))
        return super().__new__(cls,*args,evaluate=False)
    @property
    def bra(self):return _python_atom(self.args[0])
    @property
    def ket(self):return _python_atom(self.args[1])
    def doit(self,**hints):
        bra,ket=self.bra,self.ket
        if isinstance(bra,BraBase) and isinstance(ket,KetBase):return bra._inner_with(ket)
        return self
    def _adjoint_product(self):
        bra,ket=self.bra,self.ket
        if isinstance(bra,StateBase) and isinstance(ket,StateBase):return InnerProduct(Dagger(ket),Dagger(bra))
        return S.conjugate(self)
    def _sympystr(self,printer,*args):
        bra,ket=self.bra,self.ket
        if isinstance(bra,StateBase) and isinstance(ket,StateBase):
            return bra._atom_repr()[:-1]+'|'+ket._atom_repr()[1:]
        return 'InnerProduct(%s, %s)'%tuple(printer._print(a) for a in self.args)


class OuterProduct(Operator):
    def __init__(self,ket,bra):
        if not isinstance(ket,KetBase) or not isinstance(bra,BraBase):raise TypeError('outer product requires a ket and a bra')
        self.ket,self.bra=ket,bra;Atom.__init__(self,ket,bra);self.hilbert_space=ket.hilbert_space
    def _apply_to(self,state,**options):
        if isinstance(state,OuterProduct):return self.bra._inner_with(state.ket)*OuterProduct(self.ket,state.bra)
        if not isinstance(state,KetBase):return None
        return self.bra._inner_with(state)*self.ket
    def _atom_repr(self):return str(self.ket)+str(self.bra)
    def _multiply_atom(self,other):
        # |k><b| * |c>  ->  <b|c>*|k>, as the reference Mul flattens it.
        if isinstance(other,KetBase):return self.bra._inner_with(other),(self.ket,)
        return None
    def _represent(self,**options):return represent(self.ket,**options)*represent(self.bra,**options)


_ORDERING_OF_CLASSES=['Zero','One','Half','Infinity','NaN','NegativeOne','NegativeInfinity','Integer','Rational','Float','Exp1','Pi','ImaginaryUnit','Symbol','Wild','Pow','Mul','Add']


def _class_name(value):
    if isinstance(value,DaggerAtom):return 'Dagger'
    if isinstance(value,PowerAtom):return 'Pow'
    if isinstance(value,OperatorExpr) and not isinstance(value,Atom):
        if len(value._terms)>1:return 'Add'
        word,coefficient=next(iter(value._terms.items()))
        return 'Mul'
    if isinstance(value,Atom):return getattr(type(value),'_sympy_class_name',type(value).__name__)
    return type(value).__name__


def _content(value):
    if isinstance(value,OperatorExpr) and not isinstance(value,Atom):
        if len(value._terms)>1:return tuple(OperatorExpr._make({w:c}) for w,c in value._terms.items())
        word,coefficient=next(iter(value._terms.items()))
        return ((coefficient,) if coefficient!=1 else ())+word
    if isinstance(value,Atom):return tuple(value.labels)
    return ()


def _compare(a,b):
    """SymPy's ``Basic.compare`` on operator expressions and their labels."""
    if a is b:return 0
    native=[isinstance(v,S.Expr) and not isinstance(v,OperatorExpr) for v in (a,b)]
    if all(native):
        try:return int(S.sympify(a).compare(S.sympify(b)))
        except Exception:return (str(a)>str(b))-(str(a)<str(b))
    n1,n2=_class_name(a),_class_name(b)
    if n1!=n2:
        unknown=len(_ORDERING_OF_CLASSES)+1
        i1=_ORDERING_OF_CLASSES.index(n1) if n1 in _ORDERING_OF_CLASSES else unknown
        i2=_ORDERING_OF_CLASSES.index(n2) if n2 in _ORDERING_OF_CLASSES else unknown
        if i1==unknown and i2==unknown:return (n1>n2)-(n1<n2)
        return (i1>i2)-(i1<i2)
    st,ot=_content(a),_content(b)
    if len(st)!=len(ot):return (len(st)>len(ot))-(len(st)<len(ot))
    for l,r in zip(st,ot):
        c=_compare(l,r)
        if c:return c
    return 0


def _split_commutative(value):
    """``args_cnc`` of an operator expression: its scalar factor and its word."""
    if isinstance(value,OperatorExpr) and not isinstance(value,Atom) and len(value._terms)==1:
        word,coefficient=next(iter(value._terms.items()))
        if coefficient!=1:return coefficient,OperatorExpr._make({word:S.Integer(1)})
    return S.Integer(1),value


def _is_commutative(value):
    if isinstance(value,OperatorExpr):return False
    if isinstance(value,(int,float,complex)):
        # SymPy does not sympify commutator arguments.
        raise AttributeError("'%s' object has no attribute 'is_commutative'"%type(value).__name__)
    return getattr(value,'is_commutative',True) is not False


def _word_factors(value):
    if isinstance(value,OperatorExpr) and not isinstance(value,Atom) and len(value._terms)==1:
        word,coefficient=next(iter(value._terms.items()))
        if coefficient==1 and len(word)>1:return word
    return None


class Commutator(Atom):
    """The commutator ``[A,B] = A*B - B*A`` of two operators, unevaluated."""
    _sympy_class_name='Commutator'
    def __new__(cls,a,b):
        value=cls.eval(a,b)
        if value is not None:return value
        return super().__new__(cls)
    def __init__(self,a,b):
        if getattr(self,'_initialized',False):return
        self.a,self.b=a,b;super().__init__(a,b);self._initialized=True
    @classmethod
    def eval(cls,a,b):
        if not (a and b):return S.Integer(0)
        if a==b:return S.Integer(0)
        if _is_commutative(a) or _is_commutative(b):return S.Integer(0)
        ca,nca=_split_commutative(a);cb,ncb=_split_commutative(b)
        if ca!=1 or cb!=1:return ca*cb*cls(nca,ncb)
        if _compare(a,b)==1:return -1*cls(b,a)
        return None
    def _atom_repr(self):return '[%s,%s]'%(self.a,self.b)
    def _dual(self):return Commutator(Dagger(self.b),Dagger(self.a))
    def doit(self,**hints):
        A,B=self.a,self.b
        if isinstance(A,Atom) and isinstance(B,Atom):
            value=A._commutator(B) if hasattr(A,'_commutator') else None
            if value is None and hasattr(B,'_commutator'):
                value=B._commutator(A)
                if value is not None:value=-value
            if value is not None:return value.doit(**hints) if hasattr(value,'doit') else value
        A=A.doit(**hints) if hasattr(A,'doit') else A
        B=B.doit(**hints) if hasattr(B,'doit') else B
        return A*B-B*A
    def _expand_commutator(self):
        cls=type(self)
        def sub(x,y):
            value=cls(x,y)
            return value._expand_commutator() if isinstance(value,Commutator) else value
        A,B=self.a,self.b
        if isinstance(A,OperatorExpr) and not isinstance(A,Atom) and len(A._terms)>1:
            return sum((sub(OperatorExpr._make({w:c}),B) for w,c in A._terms.items()),S.Integer(0))
        if isinstance(B,OperatorExpr) and not isinstance(B,Atom) and len(B._terms)>1:
            return sum((sub(A,OperatorExpr._make({w:c})) for w,c in B._terms.items()),S.Integer(0))
        word=_word_factors(A)
        if word is not None:
            a=word[0];b=OperatorExpr._make({word[1:]:S.Integer(1)})
            return a*sub(b,B)+sub(a,B)*b
        word=_word_factors(B)
        if word is not None:
            b=word[0];c=OperatorExpr._make({word[1:]:S.Integer(1)})
            return sub(A,b)*c+b*sub(A,c)
        return self
    def expand(self,**hints):
        if hints.get('commutator'):return self._expand_commutator()
        return self


class AntiCommutator(Commutator):
    """The anticommutator ``{A,B} = A*B + B*A`` of two operators, unevaluated."""
    _sympy_class_name='AntiCommutator'
    @classmethod
    def eval(cls,a,b):
        if not (a and b):return S.Integer(0)
        if a==b:return 2*a**2
        if _is_commutative(a) or _is_commutative(b):return 2*a*b
        ca,nca=_split_commutative(a);cb,ncb=_split_commutative(b)
        if ca!=1 or cb!=1:return ca*cb*cls(nca,ncb)
        if _compare(a,b)==1:return cls(b,a)
        return None
    def _dual(self):return AntiCommutator(Dagger(self.a),Dagger(self.b))
    def doit(self,**hints):
        A,B=self.a,self.b
        if isinstance(A,Atom) and isinstance(B,Atom):
            value=A._anticommutator(B) if hasattr(A,'_anticommutator') else None
            if value is None and hasattr(B,'_anticommutator'):value=B._anticommutator(A)
            if value is not None:return value.doit(**hints) if hasattr(value,'doit') else value
        A=A.doit(**hints) if hasattr(A,'doit') else A
        B=B.doit(**hints) if hasattr(B,'doit') else B
        return A*B+B*A
    def _atom_repr(self):return '{%s,%s}'%(self.a,self.b)
    def expand(self,**hints):return self


def _tensor_factor_repr(factor):
    # Reference str printer parenthesizes Add, Mul and Pow tensor factors.
    if isinstance(factor,OperatorExpr):
        compound=len(factor._terms)>1 or any(len(w)>1 or c!=1 for w,c in factor._terms.items())
    else:compound=isinstance(factor,S.Expr) and (factor.is_Add or factor.is_Mul or factor.is_Pow)
    return '(%s)'%factor if compound else str(factor)


class TensorProductAtom(Atom):
    def __init__(self,*factors):self.factors=tuple(factors);super().__init__(*factors)
    def _atom_repr(self):return 'x'.join(map(_tensor_factor_repr,self.factors))
    def _apply_to(self,state,**options):
        if not isinstance(state,TensorProductAtom) or len(self.factors)!=len(state.factors):return None
        return TensorProduct(*(qapply(a*b,**options) for a,b in zip(self.factors,state.factors)))
    def _represent(self,**options):return S.kronecker_product(*(represent(v,**options) for v in self.factors))
    def rewrite(self,*args,**hints):
        # Reference TensorProduct._eval_rewrite: rewrite the factors, then
        # distribute the tensor product over the resulting sums.
        factors=[v.rewrite(*args,**hints) if hasattr(v,'rewrite') else v for v in self.factors]
        result=TensorProduct(*factors)
        return result.expand(tensorproduct=True) if hasattr(result,'expand') else result
    def expand(self,**hints):
        # Expr.expand (deep) then TensorProduct._eval_expand_tensorproduct:
        # distribute over the first sum among the factors, recursively.
        native_hints={k:v for k,v in hints.items() if k not in ('tensorproduct','commutator')}
        factors=[v.expand(**(hints if isinstance(v,OperatorExpr) else native_hints)) if hasattr(v,'expand') else v for v in self.factors]
        if not hints.get('tensorproduct'):
            return self if factors==list(self.factors) else TensorProduct(*factors)
        for i,v in enumerate(factors):
            if isinstance(v,S.Expr) and v.is_Add:
                result=S.Integer(0)
                for term in v.args:
                    result+=TensorProduct(*(factors[:i]+[term]+factors[i+1:])).expand(**hints)
                return result
        return TensorProduct(*factors)
    def _multiply_atom(self,other):
        # Products of equal-length tensor products combine factorwise, as the
        # reference does automatically since tensor_product_simp was deprecated.
        if not isinstance(other,TensorProductAtom) or len(self.factors)!=len(other.factors):return None
        products=[a*b for a,b in zip(self.factors,other.factors)]
        if any(isinstance(v,OperatorExpr) and len(v._terms)!=1 for v in products):return None
        combined=TensorProduct(*products)
        terms=_terms(combined)
        if len(terms)!=1:return None
        word,coefficient=next(iter(terms.items()))
        return coefficient,word


def _tensor_factor_terms(value):
    """Commuting coefficient/factor pairs of one tensor-product argument."""
    if isinstance(value,OperatorExpr):
        return [(c,OperatorExpr._make({word:S.Integer(1)})) for word,c in value._terms.items()]
    value=S.sympify(value)
    if isinstance(value,S.Expr) and value.is_commutative is False:
        # Reference TensorProduct.flatten: args_cnc moves commutative parts out.
        commutative,noncommutative=value.args_cnc()
        return [(S.Mul(*commutative),S.Mul(*noncommutative) if noncommutative else S.Integer(1))]
    return [(value,S.Integer(1))]


def TensorProduct(*args):
    if args and all(isinstance(v,S.Matrix) for v in args):return S.kronecker_product(*args)
    if not args:return S.Integer(1)
    result=S.Integer(0)
    for combination in product(*(_tensor_factor_terms(v) for v in args)):
        coefficient=S.Mul(*(item[0] for item in combination));factors=[item[1] for item in combination]
        result+=coefficient*(factors[0] if len(factors)==1 else TensorProductAtom(*factors))
    return result


def tensor_product_simp(e,**hints):
    if not isinstance(e,OperatorExpr):return e
    result=S.Integer(0)
    for word,c in e._terms.items():
        factors=[]
        for item in word:
            if factors and isinstance(factors[-1],TensorProductAtom) and isinstance(item,TensorProductAtom) and len(factors[-1].factors)==len(item.factors):
                previous=factors.pop();factors.append(TensorProduct(*(a*b for a,b in zip(previous.factors,item.factors))))
            else:factors.append(item)
        term=c
        for factor in factors:term*=factor
        result+=term
    return result


def qapply(e,**options):
    if not isinstance(e,OperatorExpr):return e
    result=S.Integer(0)
    if any(isinstance(atom,TensorProductAtom) and any(isinstance(f,OperatorExpr) for f in atom.factors) for word in e._terms for atom in word):
        # Reference qapply applies itself inside each tensor-product factor.
        rebuilt=S.Integer(0)
        for word,coefficient in e._terms.items():
            term=coefficient
            for atom in word:
                term=term*(TensorProduct(*(qapply(f,**options) for f in atom.factors)) if isinstance(atom,TensorProductAtom) else atom)
            rebuilt+=term
        if rebuilt!=e:return qapply(rebuilt,**options)
    for word,coefficient in e._terms.items():
        value=coefficient;remaining=list(word)
        for _ in range(10000):
            changed=False
            for index in range(len(remaining)-2,-1,-1):
                left,right=remaining[index:index+2]
                action=left._apply_to(right,**options) if hasattr(left,'_apply_to') else None
                if action is None:continue
                prefix=value
                for atom in remaining[:index]:prefix=prefix*atom
                term=prefix*action
                for atom in remaining[index+2:]:term=term*atom
                result+=qapply(term,**options) if isinstance(term,OperatorExpr) and term!=e else term
                changed=True;break
            if changed:break
            for atom in remaining:value=value*atom
            result+=value;break
        else:raise RuntimeError('operator application exceeded its reduction limit')
    return result


class Wavefunction(KetBase):
    def __init__(self,expression,*variables):
        self.expr=S.sympify(expression);self._limits=[]
        for value in variables:
            if isinstance(value,(tuple,list)):
                if len(value)!=3:raise ValueError('wavefunction domain requires (variable, lower, upper)')
                self._limits.append(tuple(map(S.sympify,value)))
            else:self._limits.append((S.sympify(value),-S.oo,S.oo))
        if not self._limits:raise ValueError('wavefunction needs at least one coordinate')
        self.variables=tuple(limit[0] for limit in self._limits);super().__init__(self.expr,*self.variables);self.hilbert_space=L2(S.ProductSet(*(S.Interval(a,b) for x,a,b in self._limits)))
        # Reference arguments keep a bare coordinate for an unbounded domain.
        self._coordinate_args=tuple(var if a==-S.oo and b==S.oo else S.Tuple(var,a,b) for var,a,b in self._limits)
    def _key(self):return type(self),self.expr,tuple(self._limits)
    def _atom_repr(self):return 'Wavefunction(%s, %s)'%(self.expr,', '.join(map(str,self._coordinate_args)))
    @property
    def limits(self):return {var:(a,b) for var,a,b in self._limits}
    def _atom_diff(self,var):
        # Reference Wavefunction._eval_derivative differentiates the expression.
        return Wavefunction(S.diff(self.expr,var),*self._limits)
    def _subs_parameters(self,mapping,**kwargs):return Wavefunction(self.expr.subs(mapping,**kwargs),*self._limits)
    def __call__(self,*args,**options):
        # Reference Wavefunction.__call__: zero outside the domain, keyword
        # values substitute the expression's parameters by name.
        var=self.variables
        if len(args)!=len(var):raise NotImplementedError("Incorrect number of arguments to function!")
        ct=0
        for v in var:
            lower,upper=self.limits[v]
            if isinstance(args[ct],S.Expr) and not (lower in args[ct].free_symbols or upper in args[ct].free_symbols):
                continue
            if (args[ct]<lower)==True or (args[ct]>upper)==True:
                return S.Integer(0)
            ct+=1
        expr=self.expr
        for symbol in list(expr.free_symbols):
            if str(symbol) in options.keys():
                expr=expr.subs(symbol,options[str(symbol)])
        return expr.subs(list(zip(var,args)))
    def prob(self):
        # Reference prob(): the density as a Wavefunction of the bare coordinates.
        return Wavefunction(self.expr*S.conjugate(self.expr),*self.variables)
    @property
    def norm(self):
        # Reference norm: iterated integrals of expr*conjugate(expr).
        value=self.expr*S.conjugate(self.expr)
        for var,a,b in self._limits:value=S.integrate(value,(var,a,b))
        return S.sqrt(value)
    @property
    def is_normalized(self):return S.simplify(self.norm-1)==0
    def normalize(self):return Wavefunction(self.expr/self.norm,*self._limits)


class DifferentialOperator(Operator):
    @property
    def variables(self):return tuple(self.function.args)
    def __init__(self,expression,function):self.expr,self.function=S.sympify(expression),S.sympify(function);Atom.__init__(self,self.expr,self.function);self.hilbert_space=L2(S.Reals)
    def _apply_to(self,state,**options):
        if not isinstance(state,Wavefunction):return None
        if tuple(self.function.args)!=state.variables:raise ValueError('operator and wavefunction coordinates differ')
        replacements={self.function:state.expr};stack=[self.expr]
        while stack:
            node=stack.pop()
            if node.class_name=='Derivative' and node.args[0]==self.function:
                value=state.expr
                for var,count in node.variable_count:value=S.diff(value,var,int(count))
                replacements[node]=value
            stack.extend(v for v in node.args if isinstance(v,S.Expr))
        expression=self.expr.xreplace(replacements).doit()
        return Wavefunction(expression,*state._limits)


class _RepresentationUnsupported(NotImplementedError):
    """A representation the reference builds from objects symox lacks; no fallback."""


def _represent_dispatch(atom,basis=None,**options):
    # Reference QExpr._represent: _represent_<BasisClass>(basis) or the default basis.
    from ..paulialgebra import Pauli
    from ..matrices import msigma
    if isinstance(atom,Pauli):return msigma(atom.i)
    if isinstance(atom,DaggerAtom):
        value=represent(atom.argument,**dict(options,basis=basis) if basis is not None else options);return value.conjugate().T
    if basis is not None:
        name=basis.__name__ if isinstance(basis,type) else type(basis).__name__
        method=getattr(atom,'_represent_'+name,None)
        if method is not None:return method(basis,**options)
    if hasattr(atom,'_represent'):
        return atom._represent(basis=basis,**options) if basis is not None else atom._represent(**options)
    if basis is None and hasattr(atom,'_represent_default_basis'):return atom._represent_default_basis(**options)
    raise NotImplementedError('This object does not have a default basis' if basis is None else 'cannot represent %s in the basis %s'%(atom,basis))


_REPRESENT_KEYS=('replace_none','index','unities')


def _represent_atom(atom,**options):
    if isinstance(atom,(StateBase,Operator)) and not isinstance(atom,OuterProduct):
        options['replace_none']=False
        temporary=get_basis(atom,**options)
        if temporary is not None:options['basis']=temporary
    native={k:v for k,v in options.items() if k not in _REPRESENT_KEYS}
    try:return _represent_dispatch(atom,**native)
    except _RepresentationUnsupported:raise
    except NotImplementedError as error:
        options['replace_none']=True
        try:
            if isinstance(atom,StateBase):return rep_innerproduct(atom,**options)
            if isinstance(atom,Operator):return rep_expectation(atom,**options)
        except NotImplementedError:raise error
        except IndexError:pass
        if 'basis' in native or not isinstance(atom,StateBase) and not isinstance(atom,Operator):raise error
    # Without a basis or enumerable states, keep the abstract matrix symbol.
    n=options.get('ndim',S.Symbol('N',integer=True,positive=True))
    if isinstance(atom,BraBase):return S.MatrixSymbol(str(atom),1,n)
    if isinstance(atom,KetBase):return S.MatrixSymbol(str(atom),n,1)
    return S.MatrixSymbol(str(atom),n,n)


def _flatten_scalar(value):
    if isinstance(value,S.Matrix) and value.shape==(1,1):return value[0]
    return value


def _represent_word(word,**options):
    # Reference represent(Mul): right to left, tracking the enumeration index
    # and the unities to integrate over for continuous bases.
    if len(word)==1:return _represent_atom(word[0],**options)
    options['index']=options.get('index',0)+1
    options.setdefault('unities',[])
    result=_represent_atom(word[-1],**dict(options));last=word[-1]
    for atom in reversed(word[:-1]):
        if isinstance(last,Operator) and not isinstance(last,OuterProduct):
            options['index']+=1;options['unities'].append(options['index'])
        elif isinstance(last,BraBase) and isinstance(atom,KetBase):options['index']+=1
        elif isinstance(last,KetBase) and isinstance(atom,Operator):options['unities'].append(options['index'])
        elif isinstance(last,KetBase) and isinstance(atom,BraBase):options['unities'].append(options['index'])
        result=_represent_atom(atom,**dict(options))*result;last=atom
    result=_flatten_scalar(result)
    continuous=any(isinstance(a,StateBase) and hasattr(a,'_enumerate_state') for a in word) or 'basis' in options
    return _integrate_word(word,result,**options) if continuous else result


def represent(expr,**options):
    fmt=options.pop('format','sympy')
    if not isinstance(expr,OperatorExpr):return expr
    results=[];shape=None
    for word,c in expr._terms.items():
        if not word:results.append(c);continue
        value=_represent_word(word,**options)
        value=c*value;shape=getattr(value,'shape',shape);results.append(value)
    if shape is None:return S.Add(*results)
    output=S.zeros(*shape) if all(isinstance(v,(S.Expr,S.Matrix)) for v in results) else None
    for value in results:
        if isinstance(value,S.Expr) and not isinstance(value,S.Matrix):
            if shape[0]!=shape[1]:raise ValueError('cannot add scalar identity to a nonsquare state representation')
            value=value*S.eye(shape[0])
        output=value if output is None else output+value
    if fmt in ('sympy','symox'):return output
    if fmt in ('numpy','scipy.sparse'):
        import numpy as np
        actual=np.array([[complex(output[i,j].evalf()) for j in range(output.cols)] for i in range(output.rows)])
        if fmt=='scipy.sparse':
            from scipy.sparse import csr_matrix
            return csr_matrix(actual)
        return actual
    raise ValueError('unsupported representation format')


def rep_innerproduct(expr,**options):
    if isinstance(expr,InnerProduct):return represent(expr.bra,**options)*represent(expr.ket,**options)
    if not isinstance(expr,(KetBase,BraBase)):raise TypeError('expr passed is not a Bra or Ket')
    basis=get_basis(expr,**options)
    if not isinstance(basis,StateBase):raise NotImplementedError("Can't form this representation!")
    index=options.get('index',1)
    basis_kets=enumerate_states(basis,index,2)
    if isinstance(expr,BraBase):
        bra=expr;ket=basis_kets[1] if Dagger(basis_kets[0])==expr else basis_kets[0]
    else:
        bra=Dagger(basis_kets[1]) if basis_kets[0]==expr else Dagger(basis_kets[0]);ket=expr
    return bra._inner_with(ket)

def rep_expectation(expr,**options):
    if not isinstance(expr,Operator):raise TypeError('The passed expression is not an operator')
    index=options.get('index',1)
    basis_state=get_basis(expr,**options)
    if basis_state is None or not isinstance(basis_state,StateBase):raise NotImplementedError('Could not get basis kets for this operator')
    basis_kets=enumerate_states(basis_state,index,2)
    return qapply(Dagger(basis_kets[1])*expr*basis_kets[0])

def _integrate_word(word,result,**options):
    # Reference integrate_result over the enumerated unity coordinates.
    if not isinstance(result,S.Expr) or isinstance(result,S.Matrix):return result
    options['replace_none']=True
    if 'basis' not in options:options['basis']=get_basis(word[-1],**options)
    elif not isinstance(options['basis'],StateBase):options['basis']=get_basis(word[-1],**options)
    basis=options.pop('basis',None)
    if basis is None:return result
    unities=options.pop('unities',[])
    if not unities:return result
    from .operatorset import state_to_operators
    for ket in enumerate_states(basis,unities):
        coordinate=ket.label[0]
        if coordinate in result.free_symbols:
            space=state_to_operators(basis).hilbert_space.interval
            result=S.integrate(result,(coordinate,space.inf,space.sup))
    return result

def integrate_result(orig_expr,result,**options):
    if not isinstance(result,S.Expr):return result
    if 'unities' not in options and 'limits' in options:
        for limit in options['limits']:result=S.integrate(result,limit)
        return result
    word=next(iter(orig_expr._terms)) if isinstance(orig_expr,OperatorExpr) and len(orig_expr._terms)==1 else None
    if not word:return result
    return _integrate_word(word,result,**options)

def get_basis(expr,*,basis=None,replace_none=True,**options):
    # Reference get_basis: map operators and state classes to basis states.
    from .operatorset import operators_to_state
    if basis is None and not replace_none:return None
    if basis is None:
        if isinstance(expr,KetBase):return _make_default(type(expr))
        if isinstance(expr,BraBase):return _make_default(expr.dual_class())
        if isinstance(expr,Operator):return operators_to_state(expr)
        return None
    if isinstance(basis,Operator) or isinstance(basis,type) and issubclass(basis,Operator):
        state=operators_to_state(basis)
        if state is None:return None
        return state if isinstance(state,StateBase) else _make_default(state)
    if isinstance(basis,StateBase):return basis
    if isinstance(basis,type) and issubclass(basis,StateBase):return _make_default(basis)
    return None

def _make_default(cls):
    try:return cls()
    except TypeError:return cls

def enumerate_states(*args,**options):
    """Reference enumerate_states(state, start, num_states) or (state, index_list)."""
    state=args[0]
    if len(args) not in (2,3):raise NotImplementedError('Wrong number of arguments!')
    if not isinstance(state,StateBase):raise TypeError('First argument is not a state!')
    if len(args)==3:num_states=args[2];options['start_index']=args[1]
    else:num_states=len(args[1]);options['index_list']=args[1]
    enumerate_state=getattr(state,'_enumerate_state',None)
    if enumerate_state is None:return []
    try:return enumerate_state(num_states,**options)
    except NotImplementedError:return []


def _enumerate_continuous_1D(state,num_states,**options):
    index_list=list(options.pop('index_list',[]))
    if not index_list:
        start=options.pop('start_index',1);index_list=list(range(start,start+num_states))
    return [type(state)(str(state.labels[0])+'_'+str(index)) for index in index_list]


def _postprocess_state_mul(expr):return qapply(expr)
def _postprocess_state_pow(expr):return expr.doit() if hasattr(expr,'doit') else expr


__all__=['AntiCommutator','qapply','Commutator','Dagger','HilbertSpaceError','HilbertSpace','TensorProductHilbertSpace','TensorPowerHilbertSpace','DirectSumHilbertSpace','ComplexSpace','L2','FockSpace','InnerProduct','Operator','HermitianOperator','UnitaryOperator','IdentityOperator','OuterProduct','DifferentialOperator','represent','rep_innerproduct','rep_expectation','integrate_result','get_basis','enumerate_states','KetBase','BraBase','StateBase','State','Ket','Bra','TimeDepState','TimeDepBra','TimeDepKet','OrthogonalKet','OrthogonalBra','OrthogonalState','Wavefunction','TensorProduct','tensor_product_simp','hbar','HBar','_postprocess_state_mul','_postprocess_state_pow']


# The reference defines the state classes in sympy.physics.quantum.state.
for _cls in (KetBase,BraBase,StateBase,State,Ket,Bra,TimeDepState,TimeDepBra,TimeDepKet,OrthogonalKet,OrthogonalBra,OrthogonalState,Wavefunction):
    _cls.__module__='symox.physics.quantum.state'
del _cls
