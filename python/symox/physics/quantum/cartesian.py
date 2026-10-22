"""Canonical position/momentum operators and continuous basis overlaps."""
import symox as S
from ._core import HermitianOperator,Ket,Bra,L2,hbar,Wavefunction,_enumerate_continuous_1D,Dagger


class XOp(HermitianOperator):
    axis=0
    def __init__(self,label='X'):super().__init__(label);self.hilbert_space=L2(S.Reals)
    def _commutator(self,other):return S.I*hbar if isinstance(other,PxOp) and self.axis==other.axis else None
    def _apply_to(self,state,**options):
        if isinstance(state,PositionKet3D):return state.labels[self.axis]*state
        if isinstance(state,XKet) and self.axis==0:return state.labels[0]*state
        if isinstance(state,Wavefunction) and self.axis<len(state.variables):return Wavefunction(state.variables[self.axis]*state.expr,*state._limits)
        return None
class YOp(XOp):axis=1
class ZOp(XOp):axis=2
class PxOp(HermitianOperator):
    axis=0
    def __init__(self,label='Px'):super().__init__(label);self.hilbert_space=L2(S.Reals)
    def _apply_to(self,state,**options):
        if isinstance(state,PxKet):return state.labels[0]*state
        if isinstance(state,Wavefunction):return Wavefunction(-S.I*hbar*S.diff(state.expr,state.variables[0]),*state._limits)
        return None


class XKet(Ket):
    def __init__(self,x):super().__init__(x);self.hilbert_space=L2(S.Reals)
    def _dual(self):return XBra(*self.labels)
    def _enumerate_state(self,num_states,**options):return _enumerate_continuous_1D(self,num_states,**options)
class XBra(Bra):
    def _dual(self):return XKet(*self.labels)
    def _enumerate_state(self,num_states,**options):return [Dagger(v) for v in self._dual()._enumerate_state(num_states,**options)]
    def _inner_with(self,state):
        if isinstance(state,XKet):return S.DiracDelta(self.labels[0]-state.labels[0])
        if isinstance(state,PxKet):return S.exp(S.I*state.labels[0]*self.labels[0]/hbar)/S.sqrt(2*S.pi*hbar)
        return super()._inner_with(state)
class PxKet(Ket):
    def __init__(self,p):super().__init__(p);self.hilbert_space=L2(S.Reals)
    def _dual(self):return PxBra(*self.labels)
    def _enumerate_state(self,num_states,**options):return _enumerate_continuous_1D(self,num_states,**options)
class PxBra(Bra):
    def _dual(self):return PxKet(*self.labels)
    def _enumerate_state(self,num_states,**options):return [Dagger(v) for v in self._dual()._enumerate_state(num_states,**options)]
    def _inner_with(self,state):
        if isinstance(state,PxKet):return S.DiracDelta(self.labels[0]-state.labels[0])
        if isinstance(state,XKet):return S.exp(-S.I*self.labels[0]*state.labels[0]/hbar)/S.sqrt(2*S.pi*hbar)
        return super()._inner_with(state)
class PositionState3D:
    def __init__(self,x,y,z):self.position=tuple(map(S.sympify,(x,y,z)));super().__init__(*self.position)
class PositionKet3D(PositionState3D,Ket):
    def _dual(self):return PositionBra3D(*self.labels)
class PositionBra3D(PositionState3D,Bra):
    def _dual(self):return PositionKet3D(*self.labels)
    def _inner_with(self,state):return S.Mul(*(S.DiracDelta(a-b) for a,b in zip(self.labels,state.labels))) if isinstance(state,PositionKet3D) else super()._inner_with(state)


X=XOp('X');Y=YOp('Y');Z=ZOp('Z');Px=PxOp('Px')
__all__=['XOp','YOp','ZOp','PxOp','X','Y','Z','Px','XKet','XBra','PxKet','PxBra','PositionState3D','PositionKet3D','PositionBra3D']
