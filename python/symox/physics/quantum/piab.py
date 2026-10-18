"""One-dimensional infinite square well eigenstates and Hamiltonian."""
import symox as S
from ._core import HermitianOperator,OrthogonalKet,OrthogonalBra,hbar,L2
m=S.Symbol('m');L=S.Symbol('L')
class PIABKet(OrthogonalKet):
    def __init__(self,n):
        self.n=S.sympify(n)
        if self.n.is_integer is False or self.n.is_nonpositive is True:raise ValueError('box quantum number must be a positive integer')
        super().__init__(self.n);self.hilbert_space=L2(S.Interval(0,L))
    def _dual(self):return PIABBra(self.n)
    def _represent(self,**options):
        x=options.get('x',S.Symbol('x'));return (S.sqrt(2/L)*S.sin(self.n*S.pi*x/L)).subs(options.get('subs',{}))
class PIABBra(OrthogonalBra):
    def __init__(self,n):self.n=PIABKet(n).n;super().__init__(self.n);self.hilbert_space=L2(S.Interval(0,L))
    def _dual(self):return PIABKet(self.n)
    def _represent(self,**options):return S.conjugate(PIABKet(self.n)._represent(**options))
class PIABHamiltonian(HermitianOperator):
    def _apply_to(self,state,**options):return state.n**2*S.pi**2*hbar**2/(2*m*L**2)*state if isinstance(state,PIABKet) else None
__all__=['PIABHamiltonian','PIABKet','PIABBra']
