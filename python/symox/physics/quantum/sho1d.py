"""Harmonic-oscillator ladder operators and number eigenstates."""
import symox as S
from ._core import Operator,HermitianOperator,OrthogonalKet,OrthogonalBra,StateBase,FockSpace,hbar
from .._algebra import Atom


class SHOOp(Operator):
    def __init__(self,label='a'):super().__init__(label);self.hilbert_space=FockSpace()
    def _represent(self,**options):
        n=int(options.get('ndim',4))
        if isinstance(self,NumberOp):return S.diag(*range(n))
        if isinstance(self,Hamiltonian):return hbar*S.Symbol('omega')*S.diag(*(S.Rational(2*i+1,2) for i in range(n)))
        matrix=S.zeros(n)
        for i in range(1,n):matrix[i-1,i]=S.sqrt(i)
        return matrix.T if isinstance(self,RaisingOp) else matrix
class LoweringOp(SHOOp):
    def _dual(self):return RaisingOp(*self.labels)
    def _apply_to(self,state,**options):
        if not isinstance(state,SHOKet):return None
        return S.Integer(0) if state.n==0 else S.sqrt(state.n)*SHOKet(state.n-1)
    def _commutator(self,other):return S.Integer(1) if isinstance(other,RaisingOp) else None
class RaisingOp(SHOOp):
    def _dual(self):return LoweringOp(*self.labels)
    def _apply_to(self,state,**options):return S.sqrt(state.n+1)*SHOKet(state.n+1) if isinstance(state,SHOKet) else None
class NumberOp(SHOOp):
    is_hermitian=True
    def _apply_to(self,state,**options):return state.n*state if isinstance(state,SHOKet) else None
class Hamiltonian(SHOOp):
    is_hermitian=True
    def _apply_to(self,state,**options):return hbar*S.Symbol('omega')*(state.n+S.Rational(1,2))*state if isinstance(state,SHOKet) else None


class SHOState(StateBase):
    def __init__(self,n):
        self.n=S.sympify(n)
        if self.n.is_integer is False or self.n.is_nonnegative is False:raise ValueError('oscillator level must be a nonnegative integer')
        StateBase.__init__(self,self.n);self.hilbert_space=FockSpace()
class SHOKet(OrthogonalKet,SHOState):
    def __init__(self,n):
        self.n=S.sympify(n)
        if self.n.is_integer is False or self.n.is_nonnegative is False:raise ValueError('oscillator level must be a nonnegative integer')
        super().__init__(self.n);self.hilbert_space=FockSpace()
    def _dual(self):return SHOBra(self.n)
    def _represent(self,**options):
        n=int(options.get('ndim',4))
        if not self.n.is_Integer or self.n>=n:raise ValueError('level lies outside the finite oscillator representation')
        return S.Matrix([int(i==self.n) for i in range(n)])
class SHOBra(OrthogonalBra,SHOState):
    def __init__(self,n):self.n=S.sympify(n);super().__init__(self.n);self.hilbert_space=FockSpace()
    def _dual(self):return SHOKet(self.n)
    def _represent(self,**options):return SHOKet(self.n)._represent(**options).T


a=LoweringOp('a');ad=RaisingOp('a');N=NumberOp('N');H=Hamiltonian('H')
__all__=['H','Hamiltonian','LoweringOp','N','NumberOp','RaisingOp','SHOBra','SHOKet','SHOOp','SHOState','a','ad']
