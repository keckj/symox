"""Mutable dense row algorithms over generic exact ring and field elements."""
from typing import TypeVar
from .exceptions import DMNonInvertibleMatrixError,DMShapeError
R=TypeVar('R');T=TypeVar('T')

def ddm_transpose(matrix):return [list(row) for row in zip(*matrix)]
def ddm_iadd(a,b):
    for i,row in enumerate(a):
        for j in range(len(row)):row[j]+=b[i][j]
def ddm_isub(a,b):
    for i,row in enumerate(a):
        for j in range(len(row)):row[j]-=b[i][j]
def ddm_ineg(a):
    for row in a:
        for j in range(len(row)):row[j]=-row[j]
def ddm_imul(a,b):
    for row in a:
        for j in range(len(row)):row[j]=row[j]*b
def ddm_irmul(a,b):
    for row in a:
        for j in range(len(row)):row[j]=b*row[j]
def ddm_imatmul(a,b,c):
    for i,row in enumerate(b):
        for k,value in enumerate(row):
            if value:
                for j in range(len(c[k])):a[i][j]+=value*c[k][j]
def ddm_irref(a,_partial_pivot=False):
    pivots=[];row=0;m=len(a);n=len(a[0]) if m else 0
    for col in range(n):
        pivot=next((i for i in range(row,m) if a[i][col]),None)
        if pivot is None:continue
        if _partial_pivot:pivot=max((i for i in range(row,m) if a[i][col]),key=lambda i:abs(a[i][col]))
        a[row],a[pivot]=a[pivot],a[row];value=a[row][col];a[row]=[v/value for v in a[row]]
        for i in range(m):
            if i==row:continue
            value=a[i][col]
            if value:a[i]=[v-value*w for v,w in zip(a[i],a[row])]
        pivots.append(col);row+=1
        if row==m:break
    return pivots

def ddm_irref_den(a,K):
    from .domainmatrix import DomainMatrix
    matrix=DomainMatrix(a,(len(a),len(a[0]) if a else 0),K)
    reduced,den,pivots=matrix.rref_den();a[:]=reduced.to_list();return den,list(pivots)

def ddm_idet(a,K):
    n=len(a)
    if any(len(row)!=n for row in a):raise DMShapeError('determinant requires a square matrix')
    sign=K.one;previous=K.one
    for k in range(n-1):
        pivot=next((i for i in range(k,n) if a[i][k]),None)
        if pivot is None:return K.zero
        if pivot!=k:a[k],a[pivot]=a[pivot],a[k];sign=-sign
        value=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):a[i][j]=K.exquo(value*a[i][j]-a[i][k]*a[k][j],previous)
            a[i][k]=K.zero
        previous=value
    return sign*(a[-1][-1] if n else K.one)

def ddm_iinv(ainv,a,K):
    n=len(a)
    if any(len(row)!=n for row in a):raise DMShapeError('inverse requires a square matrix')
    augmented=[list(row)+[K.one if i==j else K.zero for j in range(n)] for i,row in enumerate(a)]
    pivots=ddm_irref(augmented)
    if pivots[:n]!=list(range(n)):raise DMNonInvertibleMatrixError('matrix is singular')
    ainv[:]=[row[n:] for row in augmented]

def ddm_ilu(a):
    m=len(a);n=len(a[0]) if m else 0;swaps=[]
    for k in range(min(m,n)):
        if not a[k][k]:
            pivot=next((i for i in range(k+1,m) if a[i][k]),None)
            if pivot is None:continue
            a[k],a[pivot]=a[pivot],a[k];swaps.append((k,pivot))
        for i in range(k+1,m):
            a[i][k]=a[i][k]/a[k][k]
            for j in range(k+1,n):a[i][j]-=a[i][k]*a[k][j]
    return swaps

def ddm_ilu_split(L,U,K):
    swaps=ddm_ilu(U)
    for i,row in enumerate(U):
        for j in range(min(i,len(row))):L[i][j]=row[j];row[j]=K.zero
    return swaps

def ddm_ilu_solve(x,L,U,swaps,b):
    rhs=[list(row) for row in b]
    for i,j in swaps:rhs[i],rhs[j]=rhs[j],rhs[i]
    m=len(U);n=len(U[0]) if m else 0;columns=len(rhs[0]) if rhs else 0
    for i in range(m):
        for k in range(i):
            for j in range(columns):rhs[i][j]-=L[i][k]*rhs[k][j]
    for i in range(n,m):
        if any(rhs[i]):raise DMNonInvertibleMatrixError('inconsistent system')
    if m<n:raise DMShapeError('LU solving requires at least as many rows as columns')
    for i in reversed(range(n)):
        if not U[i][i]:raise DMNonInvertibleMatrixError('matrix is singular')
        for j in range(columns):x[i][j]=(rhs[i][j]-sum(U[i][k]*x[k][j] for k in range(i+1,n)))/U[i][i]

def ddm_berk(M,K):
    from .domainmatrix import DomainMatrix
    coefficients=DomainMatrix(M,(len(M),len(M)),K).charpoly()
    return [[v] for v in coefficients]

__all__=['R','T','ddm_transpose','ddm_iadd','ddm_isub','ddm_ineg','ddm_imul','ddm_irmul',
 'ddm_imatmul','ddm_irref','ddm_irref_den','ddm_idet','ddm_iinv','ddm_ilu','ddm_ilu_split','ddm_ilu_solve','ddm_berk']
