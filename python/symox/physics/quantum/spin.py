"""Angular-momentum eigenstates, ladder actions, rotations and two-spin coupling."""
import symox as S
from ._core import (Operator, OrthogonalKet, OrthogonalBra, Atom, OperatorExpr,
    TensorProductAtom, TensorProduct, hbar, represent)
from ..wigner import clebsch_gordan, wigner_d, wigner_d_small


def m_values(j):
    j=S.sympify(j)
    if not (2*j).is_Integer or j<0:raise ValueError('spin must be an explicit nonnegative integer or half-integer')
    dimension=int(2*j+1);return dimension,[j-i for i in range(dimension)]


class _SpinState:
    axis='z'
    def _spin_init(self,j,m):
        self.j,self.m=S.sympify(j),S.sympify(m)
        if self.j.is_number and self.m not in m_values(self.j)[1]:raise ValueError('magnetic quantum number lies outside the spin multiplet')
    def _represent(self,**options):
        _,values=m_values(self.j);column=S.Matrix([int(m==self.m) for m in values])
        if self.axis!='z':
            # x/y eigenstates are rotations of the z eigenbasis.
            rotation=wigner_d(self.j,0,-S.pi/2,0) if self.axis=='x' else wigner_d(self.j,-S.pi/2,-S.pi/2,0)
            column=rotation*column
        return column.T.conjugate() if self.is_bra else column
class JzKet(OrthogonalKet,_SpinState):
    def __init__(self,j,m):self._spin_init(j,m);super().__init__(self.j,self.m)
    def _dual(self):return _BRA[self.axis](self.j,self.m)
    def _represent(self,**options):return _SpinState._represent(self,**options)
class JzBra(OrthogonalBra,_SpinState):
    def __init__(self,j,m):self._spin_init(j,m);super().__init__(self.j,self.m)
    def _dual(self):return _KET[self.axis](self.j,self.m)
    def _represent(self,**options):return _SpinState._represent(self,**options)
    def _inner_with(self,ket):
        if isinstance(ket,_SpinState):
            if self.j!=ket.j:return S.Integer(0)
            return (self._represent()*ket._represent())[0]
        return super()._inner_with(ket)
class JxKet(JzKet):axis='x'
class JxBra(JzBra):axis='x'
class JyKet(JzKet):axis='y'
class JyBra(JzBra):axis='y'
_KET={'x':JxKet,'y':JyKet,'z':JzKet};_BRA={'x':JxBra,'y':JyBra,'z':JzBra}


class _SpinOperator(Operator):
    axis='z';is_hermitian=True
    def _apply_to(self,state,**options):
        if not isinstance(state,_SpinState) or not state.is_ket:return None
        if self.axis=='2':return hbar**2*state.j*(state.j+1)*state
        if self.axis==state.axis:return hbar*state.m*state
        if state.axis=='z' and self.axis in ('+','-'):
            step=1 if self.axis=='+' else -1
            if state.m==step*state.j:return S.Integer(0)
            return hbar*S.sqrt(state.j*(state.j+1)-state.m*(state.m+step))*JzKet(state.j,state.m+step)
        if state.axis=='z' and self.axis in ('x','y'):
            up=Jplus._apply_to(state);down=Jminus._apply_to(state)
            return (up+down)/2 if self.axis=='x' else (up-down)/(2*S.I)
        # Change the finite eigenbasis explicitly for other axes.
        vector=self._represent(j=state.j)*state._represent()
        return sum((vector[i]*JzKet(state.j,m) for i,m in enumerate(m_values(state.j)[1])),S.Integer(0))
    def _dual(self):return Jminus if self.axis=='+' else Jplus if self.axis=='-' else self
    def _commutator(self,other):
        if not isinstance(other,_SpinOperator):return None
        if self.axis=='2' or other.axis=='2' or self.axis==other.axis:return S.Integer(0)
        cyclic={('x','y'):Jz,('y','z'):Jx,('z','x'):Jy}
        if (self.axis,other.axis) in cyclic:return S.I*hbar*cyclic[self.axis,other.axis]
        if (other.axis,self.axis) in cyclic:return -S.I*hbar*cyclic[other.axis,self.axis]
        if self.axis=='+' and other.axis=='-':return 2*hbar*Jz
        if self.axis=='-' and other.axis=='+':return -2*hbar*Jz
        return None
    def _represent(self,**options):
        j=S.sympify(options.get('j',S.Rational(1,2)));dimension,values=m_values(j)
        if self.axis=='2':return hbar**2*j*(j+1)*S.eye(dimension)
        if self.axis=='z':return hbar*S.diag(*values)
        up=S.zeros(dimension)
        for column,m in enumerate(values):
            if column>0:up[column-1,column]=hbar*S.sqrt(j*(j+1)-m*(m+1))
        if self.axis=='+':return up
        if self.axis=='-':return up.T
        return (up+up.T)/2 if self.axis=='x' else (up-up.T)/(2*S.I)
class JzOp(_SpinOperator):axis='z'
class J2Op(_SpinOperator):axis='2'
class JxOp(_SpinOperator):axis='x'
class JyOp(_SpinOperator):axis='y'
class JplusOp(_SpinOperator):axis='+';is_hermitian=False
class JminusOp(_SpinOperator):axis='-';is_hermitian=False
Jz=JzOp('Jz');J2=J2Op('J2');Jx=JxOp('Jx');Jy=JyOp('Jy');Jplus=JplusOp('J+');Jminus=JminusOp('J-')


class WignerD:
    def __init__(self,j,m,mp,alpha,beta,gamma):self.args=tuple(map(S.sympify,(j,m,mp,alpha,beta,gamma)))
    def doit(self,**hints):
        j,m,mp,a,b,g=self.args;_,values=m_values(j)
        if m not in values or mp not in values:return S.Integer(0)
        return wigner_d(j,-a,-b,-g)[values.index(m),values.index(mp)]
    def _symox_(self):return self.doit()
    def subs(self,*args,**kwargs):return WignerD(*(x.subs(*args,**kwargs) for x in self.args))
class Rotation(Operator):
    is_unitary=True
    def __init__(self,alpha,beta,gamma):self.alpha,self.beta,self.gamma=map(S.sympify,(alpha,beta,gamma));super().__init__(self.alpha,self.beta,self.gamma)
    @staticmethod
    def D(j,m,mp,alpha,beta,gamma):return WignerD(j,m,mp,alpha,beta,gamma)
    @staticmethod
    def d(j,m,mp,beta):return WignerD(j,m,mp,0,beta,0)
    def _represent(self,**options):return wigner_d(options.get('j',S.Rational(1,2)),-self.alpha,-self.beta,-self.gamma)
    def _apply_to(self,state,**options):
        if not isinstance(state,JzKet):return None
        values=m_values(state.j)[1];matrix=self._represent(j=state.j);column=values.index(state.m)
        return sum((matrix[i,column]*JzKet(state.j,m) for i,m in enumerate(values)),S.Integer(0))
    def _dual(self):return Rotation(-self.gamma,-self.beta,-self.alpha)


class JzKetCoupled(JzKet):
    def __init__(self,j,m,jn,jcoupling=None):
        self.jn=tuple(map(S.sympify,jn));self.jcoupling=jcoupling
        if len(self.jn)!=2:raise NotImplementedError('coupled states currently support two angular momenta')
        self._spin_init(j,m);Atom.__init__(self,self.j,self.m,S.Tuple(*self.jn))
    def _dual(self):return _COUPLED_BRA[self.axis](self.j,self.m,self.jn,self.jcoupling)
    def _subs_parameters(self,mapping,**kwargs):return type(self)(self.j.subs(mapping,**kwargs),self.m.subs(mapping,**kwargs),[j.subs(mapping,**kwargs) for j in self.jn],self.jcoupling)
class JzBraCoupled(JzBra):
    def __init__(self,j,m,jn,jcoupling=None):
        state=JzKetCoupled(j,m,jn,jcoupling);self.jn,self.jcoupling=state.jn,state.jcoupling
        self._spin_init(j,m);Atom.__init__(self,self.j,self.m,S.Tuple(*self.jn))
    def _dual(self):return _COUPLED_KET[self.axis](self.j,self.m,self.jn,self.jcoupling)
class JxKetCoupled(JzKetCoupled):axis='x'
class JxBraCoupled(JzBraCoupled):axis='x'
class JyKetCoupled(JzKetCoupled):axis='y'
class JyBraCoupled(JzBraCoupled):axis='y'
_COUPLED_KET={'x':JxKetCoupled,'y':JyKetCoupled,'z':JzKetCoupled};_COUPLED_BRA={'x':JxBraCoupled,'y':JyBraCoupled,'z':JzBraCoupled}


def uncouple(expr,jn=None,jcoupling_list=None):
    if isinstance(expr,JzKetCoupled):
        j1,j2=expr.jn;result=S.Integer(0)
        for m1 in m_values(j1)[1]:
            m2=expr.m-m1
            if m2 not in m_values(j2)[1]:continue
            coefficient=clebsch_gordan(j1,j2,expr.j,m1,m2,expr.m)
            if coefficient:result+=coefficient*TensorProduct(_KET[expr.axis](j1,m1),_KET[expr.axis](j2,m2))
        return result
    if isinstance(expr,OperatorExpr):return sum((c*_uncouple_word(w) for w,c in expr._terms.items()),S.Integer(0))
    return expr

def _uncouple_word(word):
    result=S.Integer(1)
    for atom in word:result=result*uncouple(atom) if isinstance(atom,JzKetCoupled) else result*atom
    return result

def couple(expr,jcoupling_list=None):
    if isinstance(expr,TensorProductAtom):
        if len(expr.factors)!=2 or any(not isinstance(k,JzKet) for k in expr.factors):raise NotImplementedError('coupling requires a tensor product of two spin kets')
        a,b=expr.factors
        if a.axis!=b.axis:raise ValueError('coupled spins must use a common quantization axis')
        m=a.m+b.m;minimum=abs(a.j-b.j);maximum=a.j+b.j;result=S.Integer(0)
        for n in range(int(maximum-minimum)+1):
            j=minimum+n
            if m not in m_values(j)[1]:continue
            coefficient=clebsch_gordan(a.j,b.j,j,a.m,b.m,m)
            if coefficient:result+=coefficient*_COUPLED_KET[a.axis](j,m,(a.j,b.j))
        return result
    if isinstance(expr,OperatorExpr):
        result=S.Integer(0)
        for word,c in expr._terms.items():
            term=c
            for atom in word:term=term*couple(atom) if isinstance(atom,TensorProductAtom) else term*atom
            result+=term
        return result
    return expr

__all__=['m_values','Jplus','Jminus','Jx','Jy','Jz','J2','Rotation','WignerD','JxKet','JxBra','JyKet','JyBra','JzKet','JzBra','JzOp','J2Op','JxKetCoupled','JxBraCoupled','JyKetCoupled','JyBraCoupled','JzKetCoupled','JzBraCoupled','couple','uncouple']
