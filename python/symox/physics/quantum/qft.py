"""Finite quantum Fourier transforms and controlled phase increments."""
import symox as S
from .gate import UGate,OneQubitGate


class QFT(UGate):
    sign=1
    def __init__(self,start,finish):
        self.start,self.finish=int(start),int(finish)
        if self.start<0 or self.finish<=self.start:raise ValueError('Fourier bit interval must be nonempty')
        size=2**(self.finish-self.start)
        matrix=S.Matrix(size,size,lambda i,j:S.exp(self.sign*2*S.pi*S.I*i*j/size)/S.sqrt(size))
        super().__init__(tuple(range(self.finish-1,self.start-1,-1)),matrix)
    def _dual(self):return IQFT(self.start,self.finish) if self.sign==1 else QFT(self.start,self.finish)
class IQFT(QFT):sign=-1
class RkGate(OneQubitGate):
    def __init__(self,target,k):self.k=S.sympify(k);super().__init__(target)
    def _key(self):return type(self),self.targets,self.k
    def _local_matrix(self):return S.diag(1,S.exp(2*S.pi*S.I/2**self.k))
Rk=RkGate
__all__=['QFT','IQFT','RkGate','Rk']
