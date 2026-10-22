"""Computational basis states and Born-rule projective measurements."""
import numpy as np
import symox as S
from ._core import Ket,Bra,ComplexSpace,StateBase,represent,Dagger
from .._algebra import Atom,OperatorExpr,_terms


class Qubit(Ket):
    is_orthogonal=True
    def __init__(self,*bits):
        if len(bits)==1 and isinstance(bits[0],str):bits=tuple(bits[0])
        elif len(bits)==1 and isinstance(bits[0],(list,tuple)):bits=tuple(bits[0])
        self.qubit_values=tuple(S.sympify(v) for v in bits)
        if not bits or any(v not in (0,1) for v in self.qubit_values):raise ValueError('qubit labels must be 0 or 1')
        self.nqubits=len(bits);Atom.__init__(self,*self.qubit_values);self.hilbert_space=ComplexSpace(2**self.nqubits)
    def _atom_repr(self):return '|'+''.join(map(str,self.qubit_values))+'>'
    def __getitem__(self,index):return self.qubit_values[-int(index)-1]
    def flip(self,*bits):
        values=list(self.qubit_values)
        for bit in bits:values[-int(bit)-1]=1-values[-int(bit)-1]
        return type(self)._from_bits(values)
    @classmethod
    def _from_bits(cls,bits):return cls(*bits)
    def _dual(self):return QubitBra(*self.qubit_values)
    def _represent(self,**options):
        value=int(''.join(map(str,self.qubit_values)),2);return S.Matrix([int(i==value) for i in range(2**self.nqubits)])


class QubitBra(Qubit,Bra):
    is_ket=False
    is_bra=True
    def _atom_repr(self):return '<'+''.join(map(str,self.qubit_values))+'|'
    def _dual(self):return Qubit(*self.qubit_values)
    def _represent(self,**options):return Qubit(*self.qubit_values)._represent(**options).T


def _int_qubit_bits(args,nqubits=None):
    """Bit values of the reference IntQubit argument conventions."""
    if len(args)==1 and isinstance(args[0],Qubit):
        if nqubits is not None:raise ValueError('nqubits cannot accompany a qubit argument')
        return tuple(args[0].qubit_values)
    if not all(isinstance(a,int) or (isinstance(a,S.Expr) and a.is_Integer) for a in args):
        raise ValueError('values must be integers, got (%s)'%(tuple(type(a) for a in args),))
    args=tuple(int(a) for a in args)
    def with_nqubits(number,count):
        if count<abs(number).bit_length():raise ValueError('cannot represent %s with %s bits'%(number,count))
        return tuple((number>>i)&1 for i in reversed(range(count)))
    if nqubits is not None:
        if not (isinstance(nqubits,int) or (isinstance(nqubits,S.Expr) and nqubits.is_Integer)):
            raise ValueError('nqubits must be an integer, got (%s)'%type(nqubits))
        if len(args)!=1:raise ValueError('too many positional arguments (%s). should be (number, nqubits=n)'%(args,))
        return with_nqubits(args[0],int(nqubits))
    # One integer uses its minimal binary expansion; a second integer above
    # one is a bit count; anything else is an explicit list of bit values.
    if len(args)==1 and args[0]>1:return with_nqubits(args[0],abs(args[0]).bit_length())
    if len(args)==2 and args[1]>1:return with_nqubits(args[0],args[1])
    return args


class IntQubit(Qubit):
    def __init__(self,*args,nqubits=None):
        bits=_int_qubit_bits(args,nqubits)
        if not bits:raise ValueError('IntQubit needs at least one value')
        super().__init__(*bits)
    # IntQubit(*bits) always reads explicit bits back, as in the reference.
    def as_int(self):
        number=0
        for bit in self.qubit_values:number=2*number+int(bit)
        return S.Integer(number)
    def _atom_repr(self):return '|%s>'%self.as_int()
    def _dual(self):return IntQubitBra(*self.qubit_values)
class IntQubitBra(IntQubit,QubitBra):
    is_ket=False
    is_bra=True
    def _atom_repr(self):return '<%s|'%self.as_int()
    def _dual(self):return IntQubit(*self.qubit_values)


def qubit_to_matrix(qubit,format='sympy'):
    if isinstance(qubit,S.Matrix):return qubit
    return represent(qubit,format=format)


def matrix_to_qubit(matrix):
    matrix=S.Matrix(matrix)
    if matrix.cols!=1 or matrix.rows<2 or matrix.rows&(matrix.rows-1):raise ValueError('qubit state must be a column of power-of-two length')
    n=matrix.rows.bit_length()-1
    return sum((matrix[i]*Qubit(format(i,'0%db'%n)) for i in range(matrix.rows) if matrix[i]!=0),S.Integer(0))


def matrix_to_density(mat):
    from .density import Density
    if isinstance(mat,S.Matrix) and mat.cols!=1:raise ValueError('matrix_to_density requires a state vector')
    return Density((matrix_to_qubit(S.Matrix(mat)),1))


def _state_probabilities(state):
    vector=qubit_to_matrix(state);weights=[S.simplify(c*S.conjugate(c)) for c in vector];total=S.simplify(sum(weights))
    if total==0:raise ValueError('zero state cannot be measured')
    return vector,weights,total


def measure_all(qubit,format='sympy',normalize=True):
    vector,weights,total=_state_probabilities(qubit);n=vector.rows.bit_length()-1
    return [(Qubit(format_bits(i,n)),S.simplify(weight/total) if normalize else weight) for i,weight in enumerate(weights) if weight!=0]


def format_bits(value,n):return format(value,'0%db'%n)


def measure_partial(qubit,bits,format='sympy',normalize=True):
    vector,weights,total=_state_probabilities(qubit);n=vector.rows.bit_length()-1;groups={}
    bits=tuple(int(b) for b in bits)
    if any(b<0 or b>=n for b in bits):raise ValueError('measurement qubit index is outside the state')
    for index,coefficient in enumerate(vector):
        if coefficient==0:continue
        outcome=tuple((index>>bit)&1 for bit in bits);groups.setdefault(outcome,[]).append(index)
    results=[]
    for indices in groups.values():
        mass=S.simplify(sum(weights[i] for i in indices));collapsed=S.Matrix([vector[i]/S.sqrt(mass) if i in indices else 0 for i in range(vector.rows)])
        results.append((matrix_to_qubit(collapsed),S.simplify(mass/total) if normalize else mass))
    return results


def _choose(results,seed=None):
    rng=np.random.default_rng(seed);probabilities=np.array([float(p.evalf()) for state,p in results]);probabilities/=probabilities.sum()
    return results[int(rng.choice(len(results),p=probabilities))][0]

def measure_all_oneshot(qubit,format='sympy',seed=None):return _choose(measure_all(qubit,format),seed)
def measure_partial_oneshot(qubit,bits,format='sympy',seed=None):return _choose(measure_partial(qubit,bits,format),seed)


__all__=['Qubit','QubitBra','IntQubit','IntQubitBra','qubit_to_matrix','matrix_to_qubit','matrix_to_density','measure_all','measure_partial','measure_partial_oneshot','measure_all_oneshot']
