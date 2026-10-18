"""Finite qubit gates with explicit target and control bit semantics."""
import numpy as np
import symox as E
from ._core import UnitaryOperator,qapply,represent,ComplexSpace
from .qubit import Qubit,qubit_to_matrix,matrix_to_qubit
from .._algebra import Atom,OperatorExpr


class Gate(UnitaryOperator):
    order=None
    def __init__(self,*targets):
        self.targets=tuple(int(v) for v in targets)
        if not self.targets or any(v<0 for v in self.targets) or len(set(self.targets))!=len(self.targets):raise ValueError('gate targets must be distinct nonnegative indices')
        Atom.__init__(self,*self.targets);self.hilbert_space=ComplexSpace(2**(max(self.targets)+1))
    def _atom_repr(self):return type(self).__name__+'('+','.join(map(str,self.targets))+')'
    def _local_matrix(self):raise NotImplementedError('this abstract gate has no specified matrix')
    def _apply_to(self,state,**options):
        if not isinstance(state,Qubit) or not state.is_ket:return None
        if max(self.targets)>=state.nqubits:raise ValueError('gate target lies outside the input state')
        local=self._local_matrix();source=0
        for target in self.targets:source=(source<<1)|int(state[target])
        result=E.Integer(0)
        for output in range(local.rows):
            coefficient=local[output,source]
            if coefficient==0:continue
            bits=list(state.qubit_values)
            for i,target in enumerate(self.targets):bits[-target-1]=(output>>(len(self.targets)-i-1))&1
            result+=coefficient*Qubit(*bits)
        return result
    def _represent(self,**options):
        n=int(options.get('nqubits',options.get('ndim',max(self.targets)+1)))
        if n<max(self.targets)+1:raise ValueError('representation has too few qubits')
        columns=[qubit_to_matrix(self._apply_to(Qubit(format(i,'0%db'%n)))) for i in range(2**n)]
        return E.Matrix.hstack(*columns)
    def _multiply_atom(self,other):
        if self.order==2 and self==other:return E.Integer(1),()
        return None
    def _dual(self):
        if self.order in (1,2):return self
        if isinstance(self,CGate):return CGate(self.controls,self.gate._dual())
        return UGate(self.targets,self._local_matrix().conjugate().T)
    def __pow__(self,power):
        power=E.sympify(power)
        if self.order and power.is_Integer:
            residue=int(power)%self.order
            return E.Integer(1) if residue==0 else super().__pow__(residue)
        return super().__pow__(power)


class OneQubitGate(Gate):
    def __init__(self,target):super().__init__(target)
class TwoQubitGate(Gate):
    def __init__(self,first,second):super().__init__(first,second)
class IdentityGate(OneQubitGate):
    order=1
    def _local_matrix(self):return E.eye(2)
class HadamardGate(OneQubitGate):
    order=2
    def _local_matrix(self):return E.Matrix([[1,1],[1,-1]])/E.sqrt(2)
class XGate(OneQubitGate):
    order=2
    def _local_matrix(self):return E.Matrix([[0,1],[1,0]])
class YGate(OneQubitGate):
    order=2
    def _local_matrix(self):return E.Matrix([[0,-E.I],[E.I,0]])
class ZGate(OneQubitGate):
    order=2
    def _local_matrix(self):return E.diag(1,-1)
class PhaseGate(OneQubitGate):
    order=4
    def _local_matrix(self):return E.diag(1,E.I)
class TGate(OneQubitGate):
    order=8
    def _local_matrix(self):return E.diag(1,E.exp(E.I*E.pi/4))
class SwapGate(TwoQubitGate):
    order=2
    def _local_matrix(self):return E.Matrix([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])
class CNotGate(TwoQubitGate):
    order=2
    def __init__(self,control,target):super().__init__(control,target);self.controls=(int(control),)
    def _local_matrix(self):return E.Matrix([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])


class UGate(Gate):
    def __init__(self,targets,U):
        targets=(targets,) if isinstance(targets,(int,E.Expr)) else tuple(targets);self.matrix=E.Matrix(U)
        if self.matrix.shape!=(2**len(targets),2**len(targets)):raise ValueError('gate matrix dimension differs from its targets')
        super().__init__(*targets)
    def _key(self):return type(self),self.targets,tuple(self.matrix)
    def _local_matrix(self):return self.matrix


class CGate(Gate):
    def __init__(self,controls,gate):
        controls=(controls,) if isinstance(controls,(int,E.Expr)) else tuple(controls)
        if not isinstance(gate,Gate):raise TypeError('controlled gate requires a target gate')
        self.controls=tuple(map(int,controls));self.gate=gate
        if set(self.controls)&set(gate.targets):raise ValueError('controls cannot overlap targets')
        super().__init__(*(self.controls+gate.targets))
    def _key(self):return type(self),self.controls,self.gate
    def _apply_to(self,state,**options):
        if not isinstance(state,Qubit):return None
        return self.gate._apply_to(state,**options) if all(state[c]==1 for c in self.controls) else state
    def _local_matrix(self):raise ValueError('controlled gate uses an explicit full-system representation')
class CGateS(CGate):pass


def CPHASE(control,target):return CGate((control,),PhaseGate(target))

def normalized(state):
    matrix=qubit_to_matrix(state);norm=E.sqrt((matrix.conjugate().T*matrix)[0]);return state/norm

def gate_simp(expr):return expr.simplify() if isinstance(expr,OperatorExpr) else expr

def gate_sort(expr):return gate_simp(expr)

def random_circuit(ngates,nqubits,gate_space=(XGate,YGate,ZGate,HadamardGate,PhaseGate,TGate),seed=None):
    rng=np.random.default_rng(seed);result=E.Integer(1)
    for i in range(int(ngates)):result=gate_space[int(rng.integers(len(gate_space)))](int(rng.integers(nqubits)))*result
    return result


CNOT=CNotGate;SWAP=SwapGate;H=HadamardGate;X=XGate;Y=YGate;Z=ZGate;T=TGate;S=PhaseGate;Phase=PhaseGate
__all__=['Gate','CGate','UGate','OneQubitGate','TwoQubitGate','IdentityGate','HadamardGate','XGate','YGate','ZGate','TGate','PhaseGate','SwapGate','CNotGate','CNOT','SWAP','H','X','Y','Z','T','S','Phase','normalized','gate_sort','gate_simp','random_circuit','CPHASE','CGateS']
