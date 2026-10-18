"""Pauli, Dirac, Fourier and parallel-axis matrices."""
import operator
import symox as S


def msigma(i):
    i = operator.index(i)
    if i==1:return S.Matrix([[0,1],[1,0]])
    if i==2:return S.Matrix([[0,-S.I],[S.I,0]])
    if i==3:return S.diag(1,-1)
    raise IndexError('Pauli matrix index must be 1, 2 or 3')


def mgamma(mu,lower=False):
    mu = operator.index(mu)
    if mu==0:result = S.diag(1,1,-1,-1)
    elif mu in (1,2,3):
        block = msigma(mu);result = S.zeros(2).row_join(block).col_join((-block).row_join(S.zeros(2)))
    elif mu==5:result = S.zeros(2).row_join(S.eye(2)).col_join(S.eye(2).row_join(S.zeros(2)))
    else:raise IndexError('Dirac index must be 0, 1, 2, 3 or 5')
    return -result if lower and mu!=0 else result


def mdft(n):
    n = operator.index(n)
    if n<1:raise ValueError('Fourier matrix dimension must be positive')
    return S.Matrix(n,n,lambda i,j:S.exp(-2*S.pi*S.I*i*j/n)/S.sqrt(n))


def pat_matrix(m,dx,dy,dz):
    r = S.Matrix([dx,dy,dz]);return S.sympify(m)*((r.T*r)[0]*S.eye(3)-r*r.T)


minkowski_tensor = S.diag(1,-1,-1,-1)
__all__ = ['msigma','mgamma','mdft','pat_matrix']
