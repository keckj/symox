"""Ordered Pauli algebra with exact Levi-Civita reduction."""
import symox as S
from ._algebra import Atom,OperatorExpr


class Pauli(Atom):
    is_hermitian = True
    def __init__(self,i,label='sigma'):
        value = S.sympify(i)
        if not value.is_Integer or int(value) not in (1,2,3):raise IndexError('Pauli index must be 1, 2 or 3')
        self.i,self.label = int(value),str(label);super().__init__(self.i,self.label)
    def _atom_repr(self):return self.label+str(self.i)
    def _multiply_atom(self,other):
        if not isinstance(other,Pauli) or self.label!=other.label:return None
        if self.i==other.i:return S.Integer(1),()
        sign = 1 if (self.i,other.i) in ((1,2),(2,3),(3,1)) else -1
        return sign*S.I,(Pauli(6-self.i-other.i,self.label),)
    def __pow__(self,power):
        power = S.sympify(power)
        return self if power.is_Integer and int(power)%2 else S.Integer(1) if power.is_Integer else super().__pow__(power)


def evaluate_pauli_product(arg):return arg.simplify() if isinstance(arg,OperatorExpr) else S.sympify(arg)


__all__ = ['evaluate_pauli_product']
