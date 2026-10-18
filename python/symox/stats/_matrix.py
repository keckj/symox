"""Positive-definite matrix laws and matrix-valued Gaussian models."""
import operator
import numpy as np
import symox as S
from ._base import Distribution,RandomSymbol,Domain,PSpace,REGISTRY,HALF,scalar,positive,check,variable
from ._joint import matrix,positive_matrix,_normal_raw
from ._univariate import continuous_model,_shift

MATRIX_REGISTRY = {}


def _numeric_matrix(value):return np.array([[float(value[i,j]) for j in range(value.cols)] for i in range(value.rows)])


def _wishart_draw(rng,df,scale):
    dimension = scale.shape[0]
    triangular = np.zeros_like(scale,dtype=float)
    for i in range(dimension):
        triangular[i,i] = np.sqrt(rng.chisquare(df-i))
        triangular[i,:i] = rng.normal(size=i)
    factor = np.linalg.cholesky(scale)@triangular
    return factor@factor.T


class MatrixSupport:
    def __init__(self,rows,cols,positive=False):self.shape,self.positive = (rows,cols),positive
    def __contains__(self,value):
        try:
            value = matrix(value)
            if value.shape != self.shape:return False
            if self.positive:positive_matrix(value)
            return True
        except (ValueError,TypeError):return False
    def __repr__(self):return 'PositiveDefiniteMatrices%s' % (self.shape,) if self.positive else 'RealMatrices%s' % (self.shape,)


class MatrixDistribution(Distribution):
    kind = 'matrix'
    def __init__(self,owner):
        super().__init__();self.owner = owner;self.family = owner.family;self.args = owner.params
    @property
    def set(self):return self.owner.support
    def pdf(self,value):return self.owner.pdf(value)
    def sample(self,rng,size=()):return self.owner.sample(rng,size)


class MatrixRandomVariable:
    _stats_matrix = True
    _symox_matrixexpr = True
    is_random = True
    def __init__(self,name,family,params,shape):
        self.name,self.family,self.params,self.shape = str(name),family,tuple(params),shape
        self.rows,self.cols = shape
        self.symbol = S.MatrixSymbol(self.name,*shape)
        self.group,self.process = self,None
        self.distribution = MatrixDistribution(self)
        self.support = MatrixSupport(*shape,positive=family in ('Wishart','MatrixGamma'))
        self._entries = {}
        MATRIX_REGISTRY[self.symbol] = self
        REGISTRY[self.symbol] = self
    def __repr__(self):return self.name
    def _symox_(self):return self.symbol
    def __hash__(self):return hash(self.symbol)
    @property
    def free_symbols(self):return {self}
    @property
    def pspace(self):
        space = PSpace((self,));space.domain = Domain((self,),self.support)
        return space
    @property
    def args(self):return self.symbol,self.distribution
    def __add__(self,other):return self.symbol+other
    def __radd__(self,other):return other+self.symbol
    def __sub__(self,other):return self.symbol-other
    def __rsub__(self,other):return other-self.symbol
    def __mul__(self,other):return self.symbol*other
    def __rmul__(self,other):return other*self.symbol
    @property
    def T(self):return self.symbol.T
    def __getitem__(self,index):
        if not isinstance(index,tuple) or len(index) != 2:raise IndexError('matrix indices must be a row/column pair')
        i,j = map(operator.index,index)
        if not 0<=i<self.rows or not 0<=j<self.cols:raise IndexError('matrix index out of range')
        if (i,j) in self._entries:return self._entries[i,j]
        if self.family == 'MatrixNormal':
            mu,U,V = self.params;model = continuous_model('Normal',(mu[i,j],S.sqrt(U[i,i]*V[j,j])))
        elif self.family == 'MatrixStudentT':
            nu,mu,Omega,Sigma = self.params;scale = S.sqrt(Sigma[i,i]*Omega[j,j]/nu);base = continuous_model('StudentT',(nu,))
            model = Distribution(nu,mu[i,j],scale);model.family = 'MatrixTEntry'
            model._pdf = lambda x:base.pdf((x-mu[i,j])/scale)/scale
            model._moment = lambda n:self.joint_moment({self._entries[i,j]:int(n)})
        else:
            n,scale = self._wishart_parameters();variance = scale[i,i]*scale[j,j]-scale[i,j]**2
            if i == j:model = continuous_model('Gamma',(n/2,2*scale[i,i]))
            else:
                model = Distribution(n,scale[i,i],scale[j,j],scale[i,j]);model.family = 'WishartEntry'
                w = S.sqrt(scale[i,i]*scale[j,j]);nu = n/2;c = scale[i,j]
                model._pdf = lambda x:(S.Abs(x)/(2*w))**(nu-HALF)*S.exp(c*x/variance)*S.besselk(nu-HALF,S.Abs(x)*w/variance)/(S.sqrt(S.pi*variance)*S.gamma(nu))
                model._mgf = lambda t:(1-2*c*t-variance*t*t)**(-nu)
                model._moment = lambda q:self.joint_moment({self._entries[i,j]:int(q)})
        rv = variable('%s[%s,%s]' % (self.name,i,j),model,group=self)
        rv.matrix_index = i,j;rv.index = i*self.cols+j;self._entries[i,j] = rv
        return rv
    def _wishart_parameters(self):
        if self.family == 'Wishart':return self.params
        alpha,beta,scale = self.params
        return 2*alpha,scale*beta/2
    def mean(self):
        if self.family == 'MatrixNormal':return self.params[0]
        if self.family == 'MatrixStudentT':
            nu,mu,omega,sigma = self.params
            check(nu>1,'matrix-t mean requires degrees of freedom above one')
            return mu
        n,scale = self._wishart_parameters()
        return n*scale
    def entry_covariance(self,first,second):
        i,j = first;k,l = second
        if self.family == 'MatrixNormal':M,U,V = self.params;return U[i,k]*V[j,l]
        if self.family == 'MatrixStudentT':
            nu,M,Omega,Sigma = self.params
            check(nu>2,'matrix-t covariance requires degrees of freedom above two')
            return Sigma[i,k]*Omega[j,l]/(nu-2)
        n,scale = self._wishart_parameters()
        return n*(scale[i,k]*scale[j,l]+scale[i,l]*scale[j,k])
    def variance(self):
        if self.cols == 1:
            return S.Matrix([[self.entry_covariance((i,0),(j,0)) for j in range(self.rows)] for i in range(self.rows)])
        return S.Matrix([[S.Add(*(self.entry_covariance((i,j),(i,k)) for i in range(self.rows))) for k in range(self.cols)] for j in range(self.cols)])
    def joint_moment(self,factors):
        entries = list(factors);powers = [factors[rv] for rv in entries];degree = sum(powers)
        mean = self.mean()
        if self.family == 'MatrixNormal':
            locations = [mean[rv.matrix_index] for rv in entries]
            cov = S.Matrix([[self.entry_covariance(a.matrix_index,b.matrix_index) for b in entries] for a in entries])
            return _normal_raw(powers,locations,cov)
        if degree == 1:return mean[entries[0].matrix_index]
        if degree == 2:
            a = entries[0];b = entries[-1]
            return mean[a.matrix_index]*mean[b.matrix_index]+self.entry_covariance(a.matrix_index,b.matrix_index)
        if len(entries) == 1 and self.family in ('Wishart','MatrixGamma'):
            rv = entries[0];i,j = rv.matrix_index;n,scale = self._wishart_parameters();t = S.Dummy('wishart_t')
            mgf = (1-2*scale[i,j]*t-(scale[i,i]*scale[j,j]-scale[i,j]**2)*t*t)**(-n/2)
            return S.simplify(mgf.diff(t,degree).subs(t,0))
        if len(entries) == 1 and self.family == 'MatrixStudentT':
            nu,M,Omega,Sigma = self.params;i,j = entries[0].matrix_index
            scale = S.sqrt(Sigma[i,i]*Omega[j,j]/nu);base = continuous_model('StudentT',(nu,))
            return _shift(S.Integer(degree),M[i,j],lambda q:scale**q*base.moment(q))
        raise NotImplementedError('higher mixed moments for this matrix law are not implemented')
    def pdf(self,value):
        x = matrix(value)
        if x.shape != self.shape:raise ValueError('density argument has the wrong matrix shape')
        f,p = self.family,self.params
        if f in ('Wishart','MatrixGamma'):
            n,scale = self._wishart_parameters();dimension = self.rows
            return S.exp(-(scale.inv()*x).trace()/2)*x.det()**((n-dimension-1)/2)/(2**(dimension*n/2)*scale.det()**(n/2)*S.multigamma(n/2,dimension))
        if f == 'MatrixNormal':
            mu,U,V = p;delta = x-mu
            return S.exp(-(V.inv()*delta.T*U.inv()*delta).trace()/2)/((2*S.pi)**S.Rational(self.rows*self.cols,2)*U.det()**S.Rational(self.cols,2)*V.det()**S.Rational(self.rows,2))
        nu,mu,Omega,Sigma = p;delta = x-mu;n,p = self.shape
        normalizer = S.multigamma((nu+n+p-1)/2,p)/(S.multigamma((nu+p-1)/2,p)*S.pi**S.Rational(n*p,2)*Omega.det()**S.Rational(n,2)*Sigma.det()**S.Rational(p,2))
        return normalizer*(S.eye(n)+Sigma.inv()*delta*Omega.inv()*delta.T).det()**(-(nu+n+p-1)/2)
    def sample(self,rng,size=()):
        shape = tuple(size);count = int(np.prod(shape)) if shape else 1;draws = []
        for i in range(count):
            if self.family in ('Wishart','MatrixGamma'):
                n,scale = self._wishart_parameters();draw = _wishart_draw(rng,float(n),_numeric_matrix(scale))
            elif self.family == 'MatrixNormal':
                M,U,V = map(_numeric_matrix,self.params)
                draw = M+np.linalg.cholesky(U)@rng.normal(size=self.shape)@np.linalg.cholesky(V).T
            else:
                nu,M,Omega,Sigma = self.params;M,Omega,Sigma = map(_numeric_matrix,(M,Omega,Sigma))
                precision = _wishart_draw(rng,float(nu)+self.cols-1,np.linalg.inv(Omega))
                covariance = np.linalg.inv(precision)
                draw = M+np.linalg.cholesky(Sigma)@rng.normal(size=self.shape)@np.linalg.cholesky(covariance).T
            draws.append(draw)
        return np.asarray(draws).reshape(shape+self.shape)
    def sample_mapping(self,rng,variables):
        draw = self.sample(rng)
        return {rv._atom:scalar(float(draw[rv.matrix_index])) for rv in variables}


def Wishart(symbol,n,scale_matrix):
    n,scale = scalar(n),positive_matrix(scale_matrix)
    check(n>scale.rows-1,'nonsingular Wishart density requires n > dimension-1')
    return MatrixRandomVariable(symbol,'Wishart',(n,scale),(int(scale.rows),int(scale.cols)))


def MatrixGamma(symbol,alpha,beta,scale_matrix):
    alpha,beta,scale = scalar(alpha),scalar(beta),positive_matrix(scale_matrix)
    positive(alpha,beta);check(alpha>S.Rational(scale.rows-1,2),'matrix-gamma density requires alpha > (dimension-1)/2')
    return MatrixRandomVariable(symbol,'MatrixGamma',(alpha,beta,scale),(int(scale.rows),int(scale.cols)))


def MatrixNormal(symbol,location_matrix,scale_matrix_1,scale_matrix_2):
    M,U,V = matrix(location_matrix),positive_matrix(scale_matrix_1),positive_matrix(scale_matrix_2)
    if U.rows != M.rows or V.rows != M.cols:raise ValueError('matrix-normal row and column scale dimensions differ from location')
    return MatrixRandomVariable(symbol,'MatrixNormal',(M,U,V),(int(M.rows),int(M.cols)))


def MatrixStudentT(symbol,nu,location_matrix,scale_matrix_1,scale_matrix_2):
    nu,M,Omega,Sigma = scalar(nu),matrix(location_matrix),positive_matrix(scale_matrix_1),positive_matrix(scale_matrix_2)
    positive(nu)
    if Omega.rows != M.cols or Sigma.rows != M.rows:raise ValueError('matrix-t column and row scale dimensions differ from location')
    return MatrixRandomVariable(symbol,'MatrixStudentT',(nu,M,Omega,Sigma),(int(M.rows),int(M.cols)))


__all__ = ['MatrixGamma','Wishart','MatrixNormal','MatrixStudentT']
