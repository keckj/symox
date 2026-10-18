"""Sparse dictionary matrix storage and zero-preserving row algorithms."""
from ._storage import MatrixStorage
from .domainmatrix import DomainMatrix
class SDM(MatrixStorage,dict):
    fmt='sparse'
    def __init__(self,elemsdict,shape,domain):
        matrix=DomainMatrix(elemsdict,shape,domain);dict.__init__(self,matrix.to_dod())
        self.shape=tuple(shape);self.rows,self.cols=self.shape;self.domain=domain
    @classmethod
    def new(cls,elemsdict,shape,domain):return cls(elemsdict,shape,domain)
    def to_dod(self):return {i:dict(row) for i,row in self.items()}
    def __repr__(self):return 'SDM(%r, %r, %s)'%(dict(self),self.shape,self.domain)

def unop_dict(A,f):return {i:{j:f(v) for j,v in row.items() if f(v)} for i,row in A.items() if any(f(v) for v in row.values())}
def binop_dict(A,B,fab,fa,fb):
    result={}
    for i in set(A)|set(B):
        row={};a,b=A.get(i,{}),B.get(i,{})
        for j in set(a)|set(b):
            value=fab(a[j],b[j]) if j in a and j in b else fa(a[j]) if j in a else fb(b[j])
            if value:row[j]=value
        if row:result[i]=row
    return result

def sdm_transpose(M):
    result={}
    for i,row in M.items():
        for j,v in row.items():result.setdefault(j,{})[i]=v
    return result

def sdm_dotvec(A,B,K):return sum((v*B[j] for j,v in A.items() if j in B),K.zero)
def sdm_matvecmul(A,B,K):return {i:value for i,row in A.items() if (value:=sdm_dotvec(row,B,K))}
def sdm_matmul(A,B,K,m,o):
    result={}
    for i,row in A.items():
        output={}
        for k,a in row.items():
            for j,b in B.get(k,{}).items():output[j]=output.get(j,K.zero)+a*b
        output={j:v for j,v in output.items() if v}
        if output:result[i]=output
    return result

def sdm_matmul_exraw(A,B,K,m,o):
    # EXRAW can contain infinities: an absent zero multiplied by infinity still
    # contributes NaN, so retain full multiplication semantics for this domain.
    denseA=DomainMatrix(A,(m,max((max(r,default=-1) for r in A.values()),default=-1)+1),K)
    n=max(denseA.shape[1],max(B,default=-1)+1);result={}
    for i in range(m):
        for j in range(o):
            value=sum((A.get(i,{}).get(k,K.zero)*B.get(k,{}).get(j,K.zero) for k in range(n)),K.zero)
            if value:result.setdefault(i,{})[j]=value
    return result

def sdm_irref(A):
    rows=[dict(row) for _,row in sorted(A.items()) if row];pivots=[];result={}
    while rows:
        row=min(rows,key=lambda r:min(r));rows.remove(row);pivot=min(row);divisor=row[pivot]
        row={j:v/divisor for j,v in row.items()}
        for other in rows+list(result.values()):
            if pivot not in other:continue
            factor=other[pivot]
            for j,v in row.items():
                value=other.get(j,0)-factor*v
                if value:other[j]=value
                else:other.pop(j,None)
        rows=[r for r in rows if r];result[len(result)]=row;pivots.append(pivot)
    nonzero={}
    for i,row in result.items():
        for j in row:
            if j not in pivots:nonzero.setdefault(j,set()).add(i)
    return result,pivots,nonzero

def sdm_rref_den(A,K):
    m=max(A,default=-1)+1;n=max((max(row,default=-1) for row in A.values()),default=-1)+1
    reduced,den,pivots=DomainMatrix(A,(m,n),K).rref_den();return reduced.to_dod(),den,list(pivots)

def sdm_nullspace_from_rref(A,one,ncols,pivots,nonzero_cols):
    nonpivots=[j for j in range(ncols) if j not in pivots];basis=[]
    for free in nonpivots:
        row={free:one}
        for i,pivot in enumerate(pivots):
            if free in A.get(i,{}):row[pivot]=-A[i][free]
        basis.append(row)
    return basis,nonpivots

def sdm_particular_from_rref(A,ncols,pivots):
    return {pivot:A[i][ncols-1] for i,pivot in enumerate(pivots) if ncols-1 in A.get(i,{}) and pivot<ncols-1}
def sdm_berk(M,n,K):return {i:v for i,v in enumerate(DomainMatrix(M,(n,n),K).charpoly()) if v}
__all__=['SDM','unop_dict','binop_dict','sdm_transpose','sdm_dotvec','sdm_matvecmul',
 'sdm_matmul','sdm_matmul_exraw','sdm_irref','sdm_rref_den','sdm_nullspace_from_rref',
 'sdm_particular_from_rref','sdm_berk']
