"""Gaussian and circular beta ensembles, with normalized spectral laws."""
import operator
import numpy as np
import symox as S
from ._base import Distribution,REGISTRY,scalar,positive,check,HALF
from ._matrix import MATRIX_REGISTRY,MatrixSupport
from ._joint import matrix


def _haar(rng,n):
    z = rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
    q,r = np.linalg.qr(z);diagonal = np.diag(r)
    return q*(diagonal/np.abs(diagonal))


class EnsembleDistribution(Distribution):
    kind = 'matrix'
    def __init__(self,owner):
        super().__init__(owner.n,owner.beta);self.owner,self.family = owner,owner.family
    @property
    def set(self):return self.owner.support
    def pdf(self,value):
        if self.owner.circular:raise NotImplementedError('circular ensembles have Haar measure, not an ambient Lebesgue density; use joint_eigen_distribution for eigenphases')
        n,beta = self.owner.n,self.owner.beta;value = matrix(value)
        if beta == 4 and not isinstance(value,S.MatrixBase) and value.shape == (n,n):
            # An abstract n x n matrix stands for the quaternion matrix, as in
            # SymPy: exp(-n*Tr(H**2)) in its n + 4*n*(n-1)/2 real coordinates,
            # with the evaluated normalising integral.
            normalizer = (S.pi/n)**S.Rational(n,2)*(S.pi/(2*n))**(n*(n-1))
            return S.exp(-n*S.Trace(value**2))/normalizer
        check(value.shape == self.owner.shape,'matrix dimensions differ from ensemble')
        if beta == 1:normalizer = (4*S.pi/n)**(S.Rational(n,2))*(2*S.pi/n)**(S.Rational(n*(n-1),4))
        elif beta == 2:normalizer = 2**S.Rational(n,2)*(S.pi/n)**S.Rational(n*n,2)
        elif beta == 4:
            # Density in the n diagonal and 4*n*(n-1)/2 real quaternion
            # coordinates. The complex representation duplicates the trace.
            normalizer = (S.pi/n)**S.Rational(n,2)*(S.pi/(2*n))**(n*(n-1))
            return S.exp(-n*S.trace(value*value)/2)/normalizer
        else:raise NotImplementedError('an ambient matrix density requires beta = 1, 2 or 4')
        return S.exp(-beta*n*S.trace(value*value)/4)/normalizer
    def sample(self,rng,size=()):return self.owner.sample(rng,size)


class _Ensemble:
    _stats_matrix = True
    _symox_matrixexpr = True
    is_random = True
    circular = False
    def __init__(self,sym,dim,beta=None):
        dim = scalar(dim);check(dim.is_integer is not False and dim.is_positive is not False,'ensemble dimension must be a positive integer')
        if not dim.is_Integer:raise NotImplementedError('matrix ensembles require an explicit dimension')
        self.n = int(dim);self.beta = S.Symbol('beta',positive=True) if beta is None else scalar(beta)
        self.name = str(sym);self.family = type(self).__name__
        self.shape = (2*self.n,2*self.n) if self.beta == 4 else (self.n,self.n)
        self.rows,self.cols = self.shape;self.symbol = S.MatrixSymbol(self.name,*self.shape)
        self.group,self.process = self,None;self.support = MatrixSupport(*self.shape)
        self.distribution = EnsembleDistribution(self);MATRIX_REGISTRY[self.symbol] = self;REGISTRY[self.symbol] = self
    def _symox_(self):return self.symbol
    def __repr__(self):return self.name
    def __hash__(self):return hash(self.symbol)
    @property
    def free_symbols(self):return {self}
    @property
    def pspace(self):
        from ._base import PSpace,Domain
        result = PSpace((self,));result.domain = Domain((self,),self.support);return result
    def mean(self):return S.zeros(*self.shape)
    def variance(self):
        if self.circular:return S.zeros(*self.shape)
        if self.beta == 1:return (S.Integer(self.n+1)/self.n)*S.eye(self.rows)
        if self.beta in (2,4):return 2/(self.beta*self.n)*S.eye(self.rows)
        raise NotImplementedError('matrix variance requires a classical beta ensemble')
    def joint_eigen_distribution(self):
        n,beta = self.n,self.beta
        variables = S.symbols('lambda1:'+str(n+1),real=True)
        if self.circular:
            normalizer = (2*S.pi)**n*S.gamma(1+beta*n/2)/S.gamma(1+beta/2)**n
            vandermonde = S.Mul(*(S.Abs(S.exp(S.I*a)-S.exp(S.I*b))**beta for i,a in enumerate(variables) for b in variables[i+1:]))
            support = S.And(*(S.And(x>=0,x<=2*S.pi) for x in variables))
            value = S.Piecewise((vandermonde/normalizer,support),(0,True))
        else:
            normalizer = (2/(beta*n))**(S.Rational(n,2)+beta*n*(n-1)/4)*(2*S.pi)**S.Rational(n,2)*S.Mul(*(S.gamma(1+beta*j/2)/S.gamma(1+beta/2) for j in range(1,n+1)))
            value = S.exp(-beta*n*S.Add(*(x*x for x in variables))/4)*S.Mul(*(S.Abs(a-b)**beta for i,a in enumerate(variables) for b in variables[i+1:]))/normalizer
        return S.Lambda(variables,S.simplify(value))
    def level_spacing_distribution(self):
        beta = self.beta;s = S.Dummy('spacing',nonnegative=True)
        if beta not in (1,2,4):raise ValueError('level spacing requires a classical beta = 1, 2 or 4 ensemble')
        b = (S.gamma((beta+2)/2)/S.gamma((beta+1)/2))**2
        a = 2*b**((beta+1)/2)/S.gamma((beta+1)/2)
        return S.Lambda(s,S.Piecewise((S.simplify(a*s**beta*S.exp(-b*s*s)),s>=0),(0,True)))
    def _draw(self,rng):
        n,beta = self.n,self.beta
        if beta not in (1,2,4):raise NotImplementedError('matrix sampling requires a concrete classical beta = 1, 2 or 4 ensemble')
        if self.circular:
            u = _haar(rng,2*n if beta == 4 else n)
            if beta == 2:return u
            if beta == 1:return u@u.T
            J = np.block([[np.zeros((n,n)),np.eye(n)],[-np.eye(n),np.zeros((n,n))]])
            return (-J@u.T@J)@u
        if beta == 1:
            a = rng.normal(size=(n,n));return (a+a.T)/np.sqrt(2*n)
        if beta == 2:
            a = rng.normal(size=(n,n))+1j*rng.normal(size=(n,n));return (a+a.conj().T)/(2*np.sqrt(n))
        A = np.zeros((n,n),dtype=complex);B = np.zeros_like(A)
        A[np.diag_indices(n)] = rng.normal(scale=1/np.sqrt(2*n),size=n)
        for i in range(n):
            for j in range(i):
                A[i,j] = (rng.normal()+1j*rng.normal())/(2*np.sqrt(n));A[j,i] = A[i,j].conjugate()
                B[i,j] = (rng.normal()+1j*rng.normal())/(2*np.sqrt(n));B[j,i] = -B[i,j]
        return np.block([[A,B],[-B.conjugate(),A.conjugate()]])
    def sample(self,rng,size=()):
        if not size:return self._draw(rng)
        return np.array([self._draw(rng) for _ in range(int(np.prod(size)))]).reshape(tuple(size)+self.shape)


class GaussianEnsemble(_Ensemble):pass
class GaussianOrthogonalEnsemble(GaussianEnsemble):
    def __init__(self,sym,dim):super().__init__(sym,dim,1)
class GaussianUnitaryEnsemble(GaussianEnsemble):
    def __init__(self,sym,dim):super().__init__(sym,dim,2)
class GaussianSymplecticEnsemble(GaussianEnsemble):
    def __init__(self,sym,dim):super().__init__(sym,dim,4)
class CircularEnsemble(_Ensemble):circular = True
class CircularOrthogonalEnsemble(CircularEnsemble):
    def __init__(self,sym,dim):super().__init__(sym,dim,1)
class CircularUnitaryEnsemble(CircularEnsemble):
    def __init__(self,sym,dim):super().__init__(sym,dim,2)
class CircularSymplecticEnsemble(CircularEnsemble):
    def __init__(self,sym,dim):super().__init__(sym,dim,4)


class JointEigenDistribution:
    """Exact pushforward of a random explicit matrix under its eigenvalue map."""
    def __init__(self,mat):
        from ._queries import random_symbols
        self.matrix = matrix(mat)
        check(self.matrix.rows == self.matrix.cols,'eigenvalue distributions require a square matrix')
        if not random_symbols(self.matrix):raise ValueError('matrix contains no random entries')
        self.args = tuple(value for value,multiplicity in self.matrix.eigenvals().items() for _ in range(int(multiplicity)))
        self.dimension = len(self.args)
    def pdf(self,*values):
        from ._queries import E
        if len(values)!=self.dimension:raise ValueError('wrong number of eigenvalues')
        return E(S.Mul(*(S.DiracDelta(scalar(x)-eigenvalue) for x,eigenvalue in zip(values,self.args))))
    def __call__(self,*values):return self.pdf(*values)
    def sample(self,rng,size=()):
        from ._queries import sample
        if not size:
            from ._queries import random_symbols,_groups
            mapping = {}
            for group,rvs in _groups(random_symbols(self.matrix)).items():
                if hasattr(group,'sample_mapping'):mapping.update(group.sample_mapping(rng,rvs))
                else:
                    for rv in rvs:mapping[rv._atom] = scalar(rv.distribution.sample(rng,()))
            actual = self.matrix.subs(mapping)
            return np.linalg.eigvals(np.array([[complex(actual[i,j]) for j in range(actual.cols)] for i in range(actual.rows)]))
        return np.array([self.sample(rng) for _ in range(int(np.prod(size)))]).reshape(tuple(size)+(self.dimension,))


def joint_eigen_distribution(mat):
    if isinstance(mat,_Ensemble):return mat.joint_eigen_distribution()
    result = JointEigenDistribution(mat);variables = S.symbols('lambda1:'+str(result.dimension+1),real=True)
    return S.Lambda(variables,result.pdf(*variables))


def level_spacing_distribution(mat):
    if not isinstance(mat,_Ensemble):raise ValueError('level spacing requires a random matrix ensemble')
    return mat.level_spacing_distribution()


__all__ = ['CircularEnsemble','CircularUnitaryEnsemble','CircularOrthogonalEnsemble','CircularSymplecticEnsemble','GaussianEnsemble','GaussianUnitaryEnsemble','GaussianOrthogonalEnsemble','GaussianSymplecticEnsemble','joint_eigen_distribution','JointEigenDistribution','level_spacing_distribution']
