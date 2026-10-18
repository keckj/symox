"""Exact Grover oracle, diffusion and amplitude-amplification iterations."""
import math
import symox as S
from .gate import Gate,UGate
from .qubit import Qubit
from ._core import qapply
from .._algebra import Atom


class OracleGate(Gate):
    def __init__(self,nqubits,f):
        self.nqubits=int(nqubits);self.f=f
        if self.nqubits<1 or not callable(f):raise ValueError('oracle needs a positive qubit count and a predicate')
        super().__init__(*range(self.nqubits))
    def _key(self):return type(self),self.nqubits,self.f
    def _apply_to(self,state,**options):
        if not isinstance(state,Qubit):return None
        if state.nqubits!=self.nqubits:raise ValueError('oracle and state qubit counts differ')
        return -state if bool(self.f(state)) else state
class WGate(UGate):
    def __init__(self,nqubits):
        self.nqubits=int(nqubits);size=2**self.nqubits
        super().__init__(tuple(range(self.nqubits-1,-1,-1)),2*S.ones(size)/size-S.eye(size))


def superposition_basis(nqubits):
    nqubits=int(nqubits)
    if nqubits<1:raise ValueError('qubit count must be positive')
    return sum((Qubit(format(i,'0%db'%nqubits)) for i in range(2**nqubits)),S.Integer(0))/S.sqrt(2**nqubits)

def grover_iteration(q,oracle):return qapply(WGate(oracle.nqubits)*oracle*q)

def apply_grover(oracle,nqubits,iterations=None):
    oracle=oracle if isinstance(oracle,OracleGate) else OracleGate(nqubits,oracle)
    count=int(math.pi*math.sqrt(2**nqubits)/4) if iterations is None else int(iterations)
    state=superposition_basis(nqubits)
    for _ in range(count):state=grover_iteration(state,oracle)
    return state


__all__=['OracleGate','WGate','superposition_basis','grover_iteration','apply_grover']
