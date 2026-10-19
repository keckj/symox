"""Ordered quantum operators, states, duality and finite representations."""
from itertools import product
import symox as S
from .._algebra import Atom,OperatorExpr,PowerAtom,_terms

hbar=S.Symbol('hbar',positive=True)


class HBar:
    def __new__(cls):return hbar


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


class FockSpace(HilbertSpace):pass


class TensorProductHilbertSpace(HilbertSpace):
    def __init__(self,*spaces):
        if not spaces or any(not isinstance(v,HilbertSpace) for v in spaces):raise HilbertSpaceError('tensor factors must be Hilbert spaces')
        self.args=spaces;self.dimension=S.Mul(*(s.dimension for s in spaces))
    def __repr__(self):return '*'.join(map(str,self.args))
class DirectSumHilbertSpace(TensorProductHilbertSpace):
    def __init__(self,*spaces):super().__init__(*spaces);self.dimension=S.Add(*(s.dimension for s in spaces))
class TensorPowerHilbertSpace(HilbertSpace):
    def __init__(self,space,power):self.args=(space,S.sympify(power));self.dimension=space.dimension**S.sympify(power)


class Operator(Atom):
    is_hermitian=None
    is_unitary=None
    def __init__(self,*labels,**kwargs):
        if not labels:labels=('A',)
        super().__init__(*(S.Symbol(v) if isinstance(v,str) else S.sympify(v) for v in labels))
        self.hilbert_space=HilbertSpace()
    def _atom_repr(self):return str(self.labels[0]) if len(self.labels)==1 else 'Operator('+', '.join(map(str,self.labels))+')'
    def adjoint(self):return Dagger(self)
    @property
    def inverse(self):return self**-1
    def _apply_to(self,state,**options):return None


class HermitianOperator(Operator):is_hermitian=True
class UnitaryOperator(Operator):is_unitary=True


class IdentityOperator(Operator):
    is_identity=True
    is_hermitian=True
    is_unitary=True
    def __init__(self,N=S.oo):self.N=S.sympify(N);super().__init__('I')
    def _key(self):return type(self),self.N
    def _apply_to(self,state,**options):return state
    def _represent(self,**options):return S.eye(int(self.N if self.N!=S.oo else options.get('ndim',2)))


class StateBase(Atom):
    is_ket=False
    is_bra=False
    is_orthogonal=False
    def __init__(self,*labels,**kwargs):
        super().__init__(*(S.Symbol(v) if isinstance(v,str) else S.sympify(v) for v in labels));self.hilbert_space=HilbertSpace()
    @property
    def label(self):return self.labels
    def _inner_with(self,ket):
        if self.is_orthogonal and ket.is_orthogonal:return S.Mul(*(S.KroneckerDelta(a,b) for a,b in zip(self.labels,ket.labels))) if len(self.labels)==len(ket.labels) else S.Integer(0)
        return S._core._apply('InnerProduct',self._native_atom(),ket._native_atom(),evaluate=False)
    def _multiply_atom(self,other):
        if self.is_bra and isinstance(other,StateBase) and other.is_ket:return self._inner_with(other),()
        if self.is_ket and isinstance(other,StateBase) and other.is_bra:return S.Integer(1),(OuterProduct(self,other),)
        return None
    @property
    def dual(self):return Dagger(self)


class State(StateBase):pass
class KetBase(StateBase):
    is_ket=True
    def _atom_repr(self):return '|'+','.join(map(str,self.labels))+'>'
class BraBase(StateBase):
    is_bra=True
    def _atom_repr(self):return '<'+','.join(map(str,self.labels))+'|'
class Ket(KetBase):pass
class Bra(BraBase):pass
class TimeDepState(State):
    @property
    def time(self):return self.labels[-1]
class TimeDepKet(Ket):
    @property
    def time(self):return self.labels[-1]
class TimeDepBra(Bra):
    @property
    def time(self):return self.labels[-1]
class OrthogonalState(State):is_orthogonal=True
class OrthogonalKet(Ket):is_orthogonal=True
class OrthogonalBra(Bra):is_orthogonal=True


_DUALS={Ket:Bra,Bra:Ket,KetBase:BraBase,BraBase:KetBase,TimeDepKet:TimeDepBra,TimeDepBra:TimeDepKet,OrthogonalKet:OrthogonalBra,OrthogonalBra:OrthogonalKet,State:State}


class DaggerAtom(Atom):
    def __init__(self,argument):self.argument=argument;super().__init__(argument)
    def _atom_repr(self):return 'Dagger(%s)'%self.argument


def Dagger(expr):
    if isinstance(expr,S.Matrix):return expr.conjugate().T
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
            term=S.conjugate(coefficient)
            for factor in reversed(word):term=term*Dagger(factor)
            result+=term
        return result
    return S.conjugate(expr)


class InnerProduct:
    def __init__(self,bra,ket):
        if not isinstance(bra,BraBase) or not isinstance(ket,KetBase):raise TypeError('inner product requires a bra and a ket')
        self.bra,self.ket=bra,ket
    @property
    def args(self):return self.bra,self.ket
    def doit(self,**hints):return self.bra._inner_with(self.ket)
    def _symox_(self):return self.doit()
    def __repr__(self):return str(self.doit())


class OuterProduct(Operator):
    def __init__(self,ket,bra):
        if not isinstance(ket,KetBase) or not isinstance(bra,BraBase):raise TypeError('outer product requires a ket and a bra')
        self.ket,self.bra=ket,bra;Atom.__init__(self,ket,bra);self.hilbert_space=ket.hilbert_space
    def _apply_to(self,state,**options):return self.bra._inner_with(state)*self.ket
    def _atom_repr(self):return str(self.ket)+str(self.bra)
    def _represent(self,**options):return represent(self.ket,**options)*represent(self.bra,**options)


class Commutator(Atom):
    def __new__(cls,a,b):
        if not isinstance(a,OperatorExpr) or not isinstance(b,OperatorExpr) or a==b:return S.Integer(0)
        return super().__new__(cls)
    def __init__(self,a,b):self.a,self.b=a,b;super().__init__(a,b)
    def _atom_repr(self):return '[%s,%s]'%(self.a,self.b)
    def doit(self,**hints):
        if hasattr(self.a,'_commutator'):
            value=self.a._commutator(self.b)
            if value is not None:return value
        if hasattr(self.b,'_commutator'):
            value=self.b._commutator(self.a)
            if value is not None:return -value
        return self.a*self.b-self.b*self.a
    def expand(self,**hints):return self.doit() if hints.get('commutator') else self


class AntiCommutator(Commutator):
    def __new__(cls,a,b):
        if not isinstance(a,OperatorExpr) or not isinstance(b,OperatorExpr):return 2*a*b
        if a==b:return 2*a**2
        return Atom.__new__(cls)
    def doit(self,**hints):return self.a*self.b+self.b*self.a
    def _atom_repr(self):return '{%s,%s}'%(self.a,self.b)


class TensorProductAtom(Atom):
    def __init__(self,*factors):self.factors=tuple(factors);super().__init__(*factors)
    def _atom_repr(self):return 'x'.join(map(str,self.factors))
    def _apply_to(self,state,**options):
        if not isinstance(state,TensorProductAtom) or len(self.factors)!=len(state.factors):return None
        return TensorProduct(*(qapply(a*b,**options) for a,b in zip(self.factors,state.factors)))
    def _represent(self,**options):return S.kronecker_product(*(represent(v,**options) for v in self.factors))


def TensorProduct(*args):
    if args and all(isinstance(v,S.Matrix) for v in args):return S.kronecker_product(*args)
    if not args:return S.Integer(1)
    result=S.Integer(0)
    for combination in product(*(list(_terms(v).items()) for v in args)):
        coefficient=S.Mul(*(item[1] for item in combination));factors=[]
        for word,c in combination:
            factors.append(OperatorExpr._make({word:S.Integer(1)}))
        result+=coefficient*TensorProductAtom(*factors)
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
    def _atom_repr(self):return 'Wavefunction(%s, %s)'%(self.expr,', '.join(map(str,self.variables)))
    def __call__(self,*values):
        if len(values)!=len(self.variables):raise ValueError('wrong number of wavefunction coordinates')
        for value,(var,a,b) in zip(values,self._limits):
            condition=S.And(S.sympify(value)>=a,S.sympify(value)<=b)
            if condition==S.false:return S.Integer(0)
        return self.expr.subs(dict(zip(self.variables,values)),simultaneous=True)
    def prob(self):return S.conjugate(self.expr)*self.expr
    @property
    def norm(self):return S.sqrt(S.integrate(self.prob(),*self._limits))
    @property
    def is_normalized(self):return S.simplify(self.norm-1)==0
    def normalize(self):return Wavefunction(self.expr/self.norm,*self._limits)


class DifferentialOperator(Operator):
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


def _represent_atom(atom,**options):
    from ..paulialgebra import Pauli
    from ..matrices import msigma
    if isinstance(atom,Pauli):return msigma(atom.i)
    if isinstance(atom,DaggerAtom):
        value=represent(atom.argument,**options);return value.conjugate().T
    if hasattr(atom,'_represent'):return atom._represent(**options)
    n=options.get('ndim',S.Symbol('N',integer=True,positive=True))
    if isinstance(atom,BraBase):return S.MatrixSymbol(str(atom),1,n)
    if isinstance(atom,KetBase):return S.MatrixSymbol(str(atom),n,1)
    return S.MatrixSymbol(str(atom),n,n)


def represent(expr,**options):
    fmt=options.pop('format','sympy')
    if not isinstance(expr,OperatorExpr):return expr
    results=[];shape=None
    for word,c in expr._terms.items():
        if not word:results.append(c);continue
        matrices=[_represent_atom(atom,**options) for atom in word];value=matrices[0]
        for matrix in matrices[1:]:value=value*matrix
        value=c*value;shape=getattr(value,'shape',shape);results.append(value)
    if shape is None:return S.Add(*results)
    output=S.zeros(*shape) if all(isinstance(v,(S.Expr,S.Matrix)) for v in results) else None
    for value in results:
        if isinstance(value,S.Expr):
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
    return represent(expr,**options)

def rep_expectation(expr,**options):return represent(qapply(expr),**options)

def integrate_result(orig_expr,result,**options):
    for limit in options.get('limits',()):result=S.integrate(result,limit)
    return result

def get_basis(expr,*,basis=None,replace_none=True,**options):
    if basis is not None:return basis
    if isinstance(expr,StateBase):return expr
    return getattr(expr,'basis',None)

def enumerate_states(*args,**options):
    indices=options.get('index_list',range(options.get('start',0),options.get('stop',options.get('ndim',4))))
    if not args:raise ValueError('a state prototype or class is required')
    prototype=args[0];cls=prototype if isinstance(prototype,type) else type(prototype)
    return [cls(index) for index in indices]

def _postprocess_state_mul(expr):return qapply(expr)
def _postprocess_state_pow(expr):return expr.doit() if hasattr(expr,'doit') else expr


__all__=['AntiCommutator','qapply','Commutator','Dagger','HilbertSpaceError','HilbertSpace','TensorProductHilbertSpace','TensorPowerHilbertSpace','DirectSumHilbertSpace','ComplexSpace','L2','FockSpace','InnerProduct','Operator','HermitianOperator','UnitaryOperator','IdentityOperator','OuterProduct','DifferentialOperator','represent','rep_innerproduct','rep_expectation','integrate_result','get_basis','enumerate_states','KetBase','BraBase','StateBase','State','Ket','Bra','TimeDepState','TimeDepBra','TimeDepKet','OrthogonalKet','OrthogonalBra','OrthogonalState','Wavefunction','TensorProduct','tensor_product_simp','hbar','HBar','_postprocess_state_mul','_postprocess_state_pow']
