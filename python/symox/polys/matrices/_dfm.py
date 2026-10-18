"""Dense domain matrices backed by Symox's native domain engine."""
from .ddm import DDM
class DFM(DDM):
    is_DFM=True;is_DDM=False
    def to_dfm(self):return self.copy()
    def to_dfm_or_ddm(self):return self.copy()
    def __repr__(self):return 'DFM(%r, %r, %s)'%(list(self),self.shape,self.domain)
__all__=['DFM']
