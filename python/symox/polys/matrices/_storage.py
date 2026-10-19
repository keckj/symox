"""Dense/sparse storage protocols sharing the exact native DomainMatrix engine."""
from .domainmatrix import DomainMatrix
from .exceptions import DMShapeError,DMNonInvertibleMatrixError


class MatrixStorage:
    is_DDM=False;is_DFM=False
    def _dm(self):return DomainMatrix(self.to_dod(),self.shape,self.domain,fmt=self.fmt)
    @classmethod
    def _from_dm(cls,matrix):return cls(matrix.to_list() if cls.fmt=='dense' else matrix.to_dod(),matrix.shape,matrix.domain)
    def copy(self):return type(self)._from_dm(self._dm())
    def to_list(self):return self._dm().to_list()
    def to_list_flat(self):return [v for row in self.to_list() for v in row]
    def flat(self):return self.to_list_flat()
    def flatiter(self):return iter(self.to_list_flat())
    def to_dok(self):return {(i,j):v for i,row in self.to_dod().items() for j,v in row.items()}
    def iter_items(self):return iter(self.to_dok().items())
    def iter_values(self):return iter(self.to_dok().values())
    def nnz(self):return len(self.to_dok())
    @property
    def is_zero_matrix(self):return not self.to_dok()
    @property
    def is_diagonal(self):return all(i==j for i,j in self.to_dok())
    @property
    def is_upper(self):return all(i<=j for i,j in self.to_dok())
    @property
    def is_lower(self):return all(i>=j for i,j in self.to_dok())
    def getitem(self,i,j):return self.to_dod().get(i,{}).get(j,self.domain.zero)
    def setitem(self,i,j,value):
        if not 0<=i<self.rows or not 0<=j<self.cols:raise IndexError('matrix index out of range')
        value=self.domain.convert(value)
        if self.fmt=='dense':self[i][j]=value
        elif value:self.setdefault(i,{})[j]=value
        elif i in self:
            self[i].pop(j,None)
            if not self[i]:del self[i]
    def convert_to(self,domain):return type(self)._from_dm(self._dm().convert_to(domain))
    def add(self,other):return type(self)._from_dm(self._dm().add(other._dm()))
    def sub(self,other):return type(self)._from_dm(self._dm().sub(other._dm()))
    def neg(self):return type(self)._from_dm(self._dm().neg())
    def mul(self,scalar):return type(self)._from_dm(self._dm().mul(scalar))
    def rmul(self,scalar):return self.applyfunc(lambda v:scalar*v,self.domain)
    def matmul(self,other):return type(self)._from_dm(self._dm().matmul(other._dm()))
    def mul_elementwise(self,other):
        if self.shape!=other.shape:raise DMShapeError('matrix shapes differ')
        return type(self).from_dod({i:{j:self.getitem(i,j)*other.getitem(i,j) for j in range(self.cols)} for i in range(self.rows)},self.shape,self.domain)
    def transpose(self):return type(self)._from_dm(self._dm().transpose())
    def hstack(self,*others):return type(self)._from_dm(self._dm().hstack(*(o._dm() for o in others)))
    def vstack(self,*others):return type(self)._from_dm(self._dm().vstack(*(o._dm() for o in others)))
    def extract(self,rows,cols):return type(self)._from_dm(self._dm().extract(rows,cols))
    def extract_slice(self,rowslice,colslice):return self.extract(list(range(self.rows))[rowslice],list(range(self.cols))[colslice])
    def applyfunc(self,func,domain):
        return type(self).from_list([[func(v) for v in row] for row in self.to_list()],self.shape,domain)
    def rref(self):
        result,pivots=self._dm().rref();return type(self)._from_dm(result),list(pivots)
    def rref_den(self):
        result,den,pivots=self._dm().rref_den();return type(self)._from_dm(result),den,list(pivots)
    def inv(self):return type(self)._from_dm(self._dm().inv())
    def det(self):return self._dm().det()
    def charpoly(self):return self._dm().charpoly()
    def nullspace(self):
        reduced,den,pivots=self._dm().rref_den()
        result=reduced.nullspace_from_rref(pivots)
        return type(self)._from_dm(result),[j for j in range(self.cols) if j not in pivots]
    def nullspace_from_rref(self,pivots=None):
        if pivots is None:pivots=[next((j for j,v in enumerate(row) if v),self.cols) for row in self.to_list()];pivots=[j for j in pivots if j<self.cols]
        nonpivots=[j for j in range(self.cols) if j not in pivots]
        return type(self)._from_dm(self._dm().nullspace_from_rref(pivots)),nonpivots
    def particular(self):
        result=type(self).zeros((1,self.cols-1),self.domain)
        for i,row in enumerate(self.to_list()):
            pivot=next((j for j,v in enumerate(row[:-1]) if v),None)
            if pivot is not None:result.setitem(0,pivot,row[-1]/row[pivot])
            elif row[-1]:raise DMNonInvertibleMatrixError('inconsistent augmented matrix')
        return result
    def lu(self):
        from .dense import ddm_ilu_split
        lower=[[self.domain.one if i==j else self.domain.zero for j in range(self.rows)] for i in range(self.rows)];upper=self.to_list()
        swaps=ddm_ilu_split(lower,upper,self.domain)
        return type(self).from_list(lower,(self.rows,self.rows),self.domain),type(self).from_list(upper,self.shape,self.domain),swaps
    def lu_solve(self,rhs):
        from .dense import ddm_ilu_solve
        l,u,swaps=self.lu();solution=[[self.domain.zero]*rhs.cols for _ in range(self.cols)]
        ddm_ilu_solve(solution,l.to_list(),u.to_list(),swaps,rhs.to_list())
        return type(self).from_list(solution,(self.cols,rhs.cols),self.domain)
    def to_ddm(self):
        from .ddm import DDM
        return DDM(self.to_list(),self.shape,self.domain)
    def to_sdm(self):
        from .sdm import SDM
        return SDM(self.to_dod(),self.shape,self.domain)
    def to_dfm(self):
        from ._dfm import DFM
        return DFM(self.to_list(),self.shape,self.domain)
    def to_dfm_or_ddm(self):return self.to_dfm()
    def diagonal(self):return [self.getitem(i,i) for i in range(min(self.shape))]
    def to_flat_nz(self):
        entries=self.to_dok();keys=tuple(entries);return list(entries.values()),(keys,self.shape)
    @classmethod
    def from_flat_nz(cls,elements,data,domain):
        keys,shape=data;return cls.from_dok(dict(zip(keys,elements)),shape,domain)
    @classmethod
    def from_dod(cls,data,shape,domain):return cls._from_dm(DomainMatrix(data,shape,domain))
    @classmethod
    def from_dok(cls,data,shape,domain):
        rows={}
        for (i,j),v in data.items():rows.setdefault(i,{})[j]=v
        return cls.from_dod(rows,shape,domain)
    @classmethod
    def from_list(cls,rows,shape,domain):return cls._from_dm(DomainMatrix(rows,shape,domain))
    @classmethod
    def from_list_flat(cls,values,shape,domain):
        m,n=shape;return cls.from_list([list(values[i*n:(i+1)*n]) for i in range(m)],shape,domain)
    @classmethod
    def zeros(cls,shape,domain):return cls._from_dm(DomainMatrix.zeros(shape,domain))
    @classmethod
    def ones(cls,shape,domain):return cls._from_dm(DomainMatrix.ones(shape,domain))
    @classmethod
    def eye(cls,shape,domain):
        if isinstance(shape,int):shape=(shape,shape)
        return cls._from_dm(DomainMatrix.eye(shape,domain))
    @classmethod
    def diag(cls,diagonal,domain,shape=None):return cls._from_dm(DomainMatrix.diag(diagonal,domain,shape))
    @classmethod
    def from_ddm(cls,other):return cls.from_list(other.to_list(),other.shape,other.domain)
    def scc(self):
        from ...utilities.iterables import strongly_connected_components
        if self.rows!=self.cols:raise DMShapeError('SCC requires a square adjacency matrix')
        edges=[(i,j) for i,row in self.to_dod().items() for j in row]
        return strongly_connected_components((list(range(self.rows)),edges))
    def __add__(self,other):return self.add(other)
    def __sub__(self,other):return self.sub(other)
    def __neg__(self):return self.neg()
    def __matmul__(self,other):return self.matmul(other)
