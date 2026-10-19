"""Exact modular arithmetic and a classical simulation of Shor's order-finding step."""
from fractions import Fraction
from math import gcd,isqrt
import symox as S
from .gate import Gate
from .qubit import Qubit
from .._algebra import Atom


class OrderFindingException(Exception):pass
class CMod(Gate):
    """Reversible modular-exponentiation oracle |k,y> -> |k,y XOR a**k mod N>."""
    def __init__(self,t,a,N):
        self.t,self.a,self.N=map(int,(t,a,N))
        if self.t<1 or not 1<self.N<=2**self.t or gcd(self.a,self.N)!=1:raise ValueError('invalid modular oracle register or non-coprime base')
        self.targets=tuple(range(2*self.t));Atom.__init__(self,self.t,self.a,self.N)
    def _key(self):return type(self),self.t,self.a,self.N
    def _apply_to(self,state,**options):
        if not isinstance(state,Qubit):return None
        if state.nqubits!=2*self.t:raise ValueError('modular oracle requires two registers of t qubits')
        source=int(''.join(map(str,state.qubit_values)),2);mask=(1<<self.t)-1
        k=source>>self.t;y=(source&mask)^pow(self.a,k,self.N)
        return Qubit(format((k<<self.t)|y,'0%db'%(2*self.t)))
    def _dual(self):return self
    def _local_matrix(self):
        size=1<<(2*self.t)
        if size>4096:raise ValueError('explicit modular-oracle matrix is limited to 12 qubits')
        matrix=S.zeros(size)
        for i in range(size):
            state=self._apply_to(Qubit(format(i,'0%db'%(2*self.t))));j=int(''.join(map(str,state.qubit_values)),2);matrix[j,i]=1
        return matrix
    def _represent(self,**options):return self._local_matrix()
    def _subs_parameters(self,mapping,**kwargs):return self

def ratioize(list,N):
    """Evaluate a continued fraction, keeping convergent denominators at most N."""
    p0,p1,q0,q1=0,1,1,0;result=Fraction(0)
    for a in list:
        p0,p1=p1,int(a)*p1+p0;q0,q1=q1,int(a)*q1+q0
        if q1>int(N):break
        result=Fraction(p1,q1)
    return S.Rational(result.numerator,result.denominator)

def getr(x,y,N):
    x,y=int(x),int(y)
    if y==0:raise ZeroDivisionError('zero phase denominator')
    values=[]
    while y:quotient,remainder=divmod(x,y);values.append(quotient);x,y=y,remainder
    approximation=ratioize(values,N)
    return approximation.q

def period_find(a,N):
    a,N=int(a),int(N)
    if N<2 or gcd(a,N)!=1:raise OrderFindingException('modular order requires a base coprime to N')
    # Exact classical simulation: no claim of quantum speedup or hardware execution.
    value=1
    for r in range(1,N+1):
        value=(value*a)%N
        if value==1:return S.Integer(r)
    raise OrderFindingException('modular order was not found')

def shor(N):
    N=int(N)
    if N<4 or S.isprime(N):raise ValueError('factorization requires a composite integer at least 4')
    if N%2==0:return S.Integer(2),S.Integer(N//2)
    root=isqrt(N)
    if root*root==N:return S.Integer(root),S.Integer(root)
    for a in range(2,N):
        factor=gcd(a,N)
        if 1<factor<N:return S.Integer(factor),S.Integer(N//factor)
        order=int(period_find(a,N))
        if order%2:continue
        factor=gcd(pow(a,order//2,N)-1,N)
        if 1<factor<N:return S.Integer(factor),S.Integer(N//factor)
    raise OrderFindingException('no nontrivial factor was found')
__all__=['CMod','OrderFindingException','getr','period_find','ratioize','shor']
