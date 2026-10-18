"""Mutable dense row storage over native exact coefficient domains."""
from ._storage import MatrixStorage
from .domainmatrix import DomainMatrix
class DDM(MatrixStorage,list):
    fmt='dense';is_DDM=True
    def __init__(self,rowslist,shape,domain):
        matrix=DomainMatrix(rowslist,shape,domain);list.__init__(self,matrix.to_list())
        self.shape=tuple(shape);self.rows,self.cols=self.shape;self.domain=domain
    def to_list(self):return [list(row) for row in self]
    def to_dod(self):return {i:{j:v for j,v in enumerate(row) if v} for i,row in enumerate(self) if any(row)}
    def __repr__(self):return 'DDM(%r, %r, %s)'%(list(self),self.shape,self.domain)
__all__=['DDM']
