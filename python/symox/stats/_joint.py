"""Joint laws, component dependence and exact multivariate moments."""
from functools import lru_cache
from itertools import product
import operator
import numpy as np
import symox as S
from ._base import Distribution,RandomSymbol,PSpace,Domain,HALF,scalar,positive,check,variable,integral
from ._univariate import continuous_model,finite_model,discrete_model,_normal_moment,_exists


def matrix(value):
    return value if getattr(value,'_symox_matrixexpr',False) else S.ImmutableMatrix(value)


def positive_matrix(value):
    result = matrix(value)
    if result.rows != result.cols:raise ValueError('scale/covariance matrix must be square')
    if isinstance(result,S.Matrix):
        for i in range(result.rows):
            for j in range(i):check(S.simplify(result[i,j]-result[j,i]).is_zero is not False,'scale/covariance matrix must be symmetric')
        for size in range(1,result.rows+1):
            minor = result[:size,:size].det()
            check(minor.is_nonpositive is not True,'scale/covariance matrix must be positive definite')
    return result


def _normal_raw(powers,mean,cov):
    @lru_cache(maxsize=None)
    def recur(counts):
        if not any(counts):return S.Integer(1)
        i = next(j for j,p in enumerate(counts) if p)
        remainder = list(counts);remainder[i] -= 1
        value = mean[i]*recur(tuple(remainder))
        for j,power in enumerate(remainder):
            if power:
                lower = remainder[:];lower[j] -= 1
                value += power*cov[i,j]*recur(tuple(lower))
        return value
    return S.simplify(recur(tuple(powers)))


def _compositions(n,k):
    if k == 1:yield (n,);return
    for first in range(n+1):
        for tail in _compositions(n-first,k-1):yield (first,)+tail


def _ewens_counts(n):
    def visit(index,remaining,values):
        if index == 1:
            yield tuple([remaining]+values)
            return
        for count in range(remaining//index+1):yield from visit(index-1,remaining-index*count,[count]+values)
    yield from visit(n,n,[])


class JointDistribution(Distribution):
    def __init__(self,owner):
        super().__init__();self.owner = owner;self.family = owner.family;self.kind = owner.kind;self.args = owner.params
    def pdf(self,*values):return self.owner.pdf(*values)
    pmf = pdf
    @property
    def set(self):return self.owner.support
    def sample(self,rng,size=()):return self.owner.sample(rng,size)


class JointVariable:
    _stats_joint = True
    is_random = True
    def __init__(self,name,family,params,dimension,kind='continuous'):
        self.name,self.family,self.params,self.dimension,self.kind = str(name),family,tuple(params),operator.index(dimension),kind
        self.symbol = S.Symbol(self.name)
        self.support = S.ProductSet(*([S.Reals]*self.dimension))
        self.distribution = JointDistribution(self)
        self.components = []
    def initialize(self,models):
        self.components = tuple(variable('%s[%s]' % (self.name,i),model,group=self) for i,model in enumerate(models))
        for i,rv in enumerate(self.components):rv.index = i
        return self
    def __repr__(self):return self.name
    def __getitem__(self,index):return self.components[operator.index(index)]
    def __iter__(self):return iter(self.components)
    def __len__(self):return self.dimension
    @property
    def free_symbols(self):return set(self.components)
    @property
    def pspace(self):
        result = PSpace(self.components);result.domain = Domain(self.components,self.support)
        return result
    def conditioned_replacements(self,bindings):
        observed = [i for i,v in enumerate(self.components) if v._atom in bindings]
        if not observed:return None
        selected = [i for i in range(self.dimension) if i not in observed]
        replacements = {self.components[i]._atom:bindings[self.components[i]._atom] for i in observed}
        if not selected:return replacements
        if self.family == 'MultivariateNormal':
            mu,cov = self.params
            Coo = S.Matrix([[cov[i,j] for j in observed] for i in observed])
            Cfo = S.Matrix([[cov[i,j] for j in observed] for i in selected])
            residual = S.Matrix([bindings[self.components[i]._atom]-mu[i] for i in observed])
            mean = S.Matrix([mu[i] for i in selected])+Cfo*Coo.inv()*residual
            conditional_cov = S.Matrix([[cov[i,j] for j in selected] for i in selected])-Cfo*Coo.inv()*Cfo.T
            model = MultivariateNormal(self.name+'|condition',mean,conditional_cov)
            replacements.update({self.components[i]._atom:scalar(model[j]) for j,i in enumerate(selected)})
            return replacements
        if self.family == 'NormalGamma':
            mu,lamda,alpha,beta = self.params
            if observed == [1]:
                tau = bindings[self.components[1]._atom];positive(tau)
                model = continuous_model('Normal',(mu,1/S.sqrt(lamda*tau)))
            else:
                value = bindings[self.components[0]._atom]
                model = continuous_model('Gamma',(alpha+HALF,1/(beta+lamda*(value-mu)**2/2)))
            replacements[self.components[selected[0]]._atom] = scalar(variable(self.name+'|condition',model))
            return replacements
        raise NotImplementedError('point conditioning for this joint family requires a conditional joint model')
    def finite_worlds(self,variables):
        if self.family not in ('Multinomial','MultivariateEwens'):return None
        n = int(self.params[0])
        counts = _compositions(n,self.dimension) if self.family == 'Multinomial' else _ewens_counts(n)
        return (({rv._atom:S.Integer(v) for rv,v in zip(self.components,row)},self.pdf(*row)) for row in counts)
    def pdf(self,*values):
        if len(values) == 1 and isinstance(values[0],(list,tuple,S.Matrix)):values = tuple(values[0])
        if len(values) != self.dimension:raise ValueError('wrong number of joint density coordinates')
        x = tuple(map(scalar,values));f,p = self.family,self.params;k = self.dimension
        if f in ('MultivariateNormal','MultivariateT','MultivariateLaplace'):
            mu,cov = p[:2];vector = S.Matrix(x)
            if f == 'MultivariateNormal':
                delta = vector-mu
                # The reference takes the exponential of the 1x1 matrix
                # -x.T*(Sigma**-1*x)/2, which yields the expanded quadratic.
                return S.exp(S.expand((S.Rational(-1,2)*delta.T*(cov.inv()*delta))[0]))/S.sqrt((2*S.pi)**k*cov.det())
            if f == 'MultivariateT':
                nu = p[2];delta = vector-mu
                return S.gamma((nu+k)/2)/(S.gamma(nu/2)*(nu*S.pi)**(S.Rational(k,2))*S.sqrt(cov.det()))*(1+(delta.T*cov.inv()*delta)[0]/nu)**(-(nu+k)/2)
            inv = cov.inv();q = (vector.T*inv*vector)[0];a = (mu.T*inv*mu)[0];v = 1-S.Rational(k,2)
            return 2*(q/(2+a))**(v/2)*S.besselk(v,S.sqrt((2+a)*q))*S.exp((vector.T*inv*mu)[0])/((2*S.pi)**S.Rational(k,2)*S.sqrt(cov.det()))
        if f == 'Dirichlet':
            alpha, = p
            return S.gamma(S.Add(*alpha))*S.Mul(*(v**(a-1)/S.gamma(a) for v,a in zip(x,alpha)))
        if f == 'Multinomial':
            n,probabilities = p
            body = S.factorial(n)*S.Mul(*(prob**v/S.factorial(v) for v,prob in zip(x,probabilities)))
            return S.Piecewise((body,S.And(S.Eq(S.Add(*x),n),*(v>=0 for v in x))),(0,True))
        if f == 'NegativeMultinomial':
            r,probabilities = p;p0 = 1-S.Add(*probabilities)
            return S.gamma(r+S.Add(*x))/S.gamma(r)*p0**r*S.Mul(*(prob**v/S.factorial(v) for v,prob in zip(x,probabilities)))
        if f == 'MultivariateEwens':
            n,theta = p
            body = S.factorial(n)/S.rf(theta,n)*S.Mul(*(theta**v/((i+1)**v*S.factorial(v)) for i,v in enumerate(x)))
            return S.Piecewise((body,S.Eq(S.Add(*((i+1)*v for i,v in enumerate(x))),n)),(0,True))
        if f == 'NormalGamma':
            mu,lamda,alpha,beta = p;value,tau = x
            return beta**alpha*S.sqrt(lamda)/(S.gamma(alpha)*S.sqrt(2*S.pi))*tau**(alpha-HALF)*S.exp(-1*beta*tau)*S.exp(-1*(lamda*tau*(value-mu)**2)/S.Integer(2))
        if f == 'GeneralizedMultivariateLogGamma':
            # The reference sums over a plain Symbol n (not a Dummy) with
            # gamma(n + 1) in the denominator.
            delta,v,scales,rates = p;n = S.Symbol('n',negative=False,integer=True)
            summand = (1-delta)**n/((S.gamma(v+n)**(k-1))*S.gamma(v)*S.gamma(n+1))*S.Mul(*(rate*scale**(-v-n) for rate,scale in zip(rates,scales)))*S.exp((v+n)*S.Add(*(rate*value for rate,value in zip(rates,x)))-S.Add(*(S.exp(rate*value)/scale for rate,value,scale in zip(rates,x,scales))))
            return summand.subs(n,0) if delta == 1 else delta**v*S.Sum(summand,(n,0,S.oo))
        if f == 'JointHandmade':return self._lambda(*x)
        raise ValueError('unknown joint probability family')
    def joint_moment(self,factors):
        powers = tuple(int(factors.get(rv,0)) for rv in self.components)
        total = sum(powers);f,p,k = self.family,self.params,self.dimension
        if total == 0:return S.Integer(1)
        if f == 'GeneralizedMultivariateLogGamma' and (len(factors) == 1 or p[0] == 1):
            return S.Mul(*(rv.distribution.moment(power) for rv,power in factors.items()))
        if f == 'MultivariateNormal':return _normal_raw(powers,p[0],p[1])
        if f == 'Dirichlet':
            alpha, = p
            return S.Mul(*(S.rf(a,n) for a,n in zip(alpha,powers)))/S.rf(S.Add(*alpha),total)
        if f in ('Multinomial','NegativeMultinomial'):
            parameter,probabilities = p;terms = []
            for orders in product(*(range(n+1) for n in powers)):
                weight = S.Mul(*(S.stirling(n,r,kind=2) for n,r in zip(powers,orders)))
                if weight == 0:continue
                if f == 'Multinomial':value = S.ff(parameter,sum(orders))*S.Mul(*(q**r for q,r in zip(probabilities,orders)))
                else:value = S.rf(parameter,sum(orders))*S.Mul(*((q/(1-S.Add(*probabilities)))**r for q,r in zip(probabilities,orders)))
                terms.append(weight*value)
            return S.simplify(S.Add(*terms))
        if f == 'MultivariateEwens':
            return S.simplify(S.Add(*(weight*S.Mul(*(values[rv._atom]**n for rv,n in factors.items())) for values,weight in self.finite_worlds(self.components))))
        if f == 'NormalGamma':
            mu,lamda,alpha,beta = p;a,b = powers
            terms = []
            for j in range(a//2+1):
                exponent = b-j
                gamma_moment = beta**(-exponent)*S.gamma(alpha+exponent)/S.gamma(alpha)
                coefficient = S.binomial(a,2*j)*mu**(a-2*j)*S.factorial(2*j)/(2**j*S.factorial(j)*lamda**j)
                terms.append(S.Piecewise((coefficient*gamma_moment,alpha+exponent>0),(S.oo,True)))
            result = S.simplify(S.Add(*terms))
            return S.Piecewise((result,alpha+b>S.Rational(a,2)),(S.nan if a%2 else S.oo,True))
        if f in ('MultivariateT','MultivariateLaplace'):
            mu,cov = p[:2];terms = []
            for selected in product(*(range(n+1) for n in powers)):
                degree = sum(selected)
                if degree%2:continue
                normal = _normal_raw(selected,[S.Integer(0)]*k,cov)
                coefficient = S.Mul(*(S.binomial(n,q)*mu[i]**(n-q) for i,(n,q) in enumerate(zip(powers,selected))))
                if f == 'MultivariateT':
                    nu = p[2];scale = (nu/2)**S.Rational(degree,2)*S.gamma((nu-degree)/2)/S.gamma(nu/2)
                else:scale = S.factorial(total-degree//2)
                terms.append(coefficient*normal*scale)
            result = S.simplify(S.Add(*terms))
            return S.Piecewise((result,p[2]>total),(S.nan if any(n%2 for n in powers) else S.oo,True)) if f == 'MultivariateT' else result
        return self.expectation(S.Mul(*(rv._atom**n for rv,n in factors.items())),list(factors))
    def expectation(self,expr,variables):
        coordinates = [S.Dummy('joint_x%s' % i,real=True) for i in range(self.dimension)]
        body = scalar(expr).subs({rv._atom:x for rv,x in zip(self.components,coordinates)},simultaneous=True)
        if self.family == 'Dirichlet':
            body = body.subs(coordinates[-1],1-S.Add(*coordinates[:-1]))
            body *= self.pdf(*coordinates).subs(coordinates[-1],1-S.Add(*coordinates[:-1]))
            for i in range(self.dimension-2,-1,-1):body = integral(body,coordinates[i],S.Integer(0),1-S.Add(*coordinates[:i]))
            return body
        body *= self.pdf(*coordinates)
        for rv,x in zip(reversed(self.components),reversed(coordinates)):
            if self.kind == 'discrete':body = S.Sum(body,(x,rv.distribution.lower,rv.distribution.upper))
            else:body = integral(body,x,rv.distribution.lower,rv.distribution.upper)
        return body
    def normal_linear(self,expr,variables):
        if self.family != 'MultivariateNormal':return None
        atoms = [rv._atom for rv in self.components]
        try:poly = S.Poly(expr,*atoms)
        except ValueError:return None
        if poly.total_degree() > 1:return None
        coefficients = S.Matrix([expr.diff(atom) for atom in atoms]);mu,cov = self.params
        mean = expr.subs(dict(zip(atoms,tuple(mu))),simultaneous=True)
        return mean,S.simplify((coefficients.T*cov*coefficients)[0])
    def independent_components(self,a,b):
        return self.family == 'MultivariateNormal' and self.params[1][a.index,b.index] == 0
    def sample(self,rng,size=()):
        f,p,k = self.family,self.params,self.dimension
        shape = tuple(size)
        if f in ('MultivariateNormal','MultivariateT','MultivariateLaplace'):
            mu = np.array([float(v) for v in p[0]])
            cov = np.array([[float(p[1][i,j]) for j in range(k)] for i in range(k)])
            if f == 'MultivariateNormal':return rng.multivariate_normal(mu,cov,size=shape)
            z = rng.multivariate_normal(np.zeros(k),cov,size=shape)
            if f == 'MultivariateT':return mu+z/np.sqrt(rng.chisquare(float(p[2]),size=shape)[...,None]/float(p[2]))
            mixing = rng.exponential(size=shape)[...,None]
            return mixing*mu+np.sqrt(mixing)*z
        if f == 'Dirichlet':return rng.dirichlet([float(a) for a in p[0]],size=shape)
        if f == 'Multinomial':return rng.multinomial(int(p[0]),[float(q) for q in p[1]],size=shape)
        if f == 'NegativeMultinomial':
            r,probabilities = p
            rate = rng.gamma(float(r),1/float(1-S.Add(*probabilities)),size=shape)
            return rng.poisson(rate[...,None]*np.array([float(q) for q in probabilities]))
        if f == 'NormalGamma':
            mu,lamda,alpha,beta = map(float,p)
            tau = rng.gamma(alpha,1/beta,size=shape)
            x = rng.normal(mu,1/np.sqrt(lamda*tau),size=shape)
            return np.stack((x,tau),axis=-1)
        if f == 'GeneralizedMultivariateLogGamma':
            delta,v,scales,rates = p
            latent = rng.negative_binomial(float(v),float(delta),size=shape)
            gammas = rng.gamma(float(v)+latent[...,None],1,size=shape+(k,))
            return np.log(gammas*np.array(list(map(float,scales))))/np.array(list(map(float,rates)))
        if f == 'MultivariateEwens':
            n,theta = int(p[0]),float(p[1]);draws = []
            for trial in range(int(np.prod(shape)) if shape else 1):
                tables = []
                for customer in range(n):
                    choice = rng.random()*(theta+customer)
                    if choice < theta:tables.append(1)
                    else:
                        choice -= theta
                        for i,count in enumerate(tables):
                            if choice < count:tables[i] += 1;break
                            choice -= count
                draws.append([tables.count(j) for j in range(1,n+1)])
            return np.array(draws).reshape(shape+(n,))
        raise NotImplementedError('sampling a handmade joint density requires a specialized joint sampler')
    def sample_mapping(self,rng,variables):
        draw = self.sample(rng,())
        return {rv._atom:scalar(np.asarray(draw[rv.index]).item()) for rv in variables}


def _vector(values):
    result = matrix(values)
    if result.cols != 1:raise ValueError('a location vector must be a column vector')
    if not getattr(result.rows,'is_Integer',isinstance(result.rows,int)):
        # SymPy holds a density of symbolic dimension as the matrix element
        # (exp(matrix)/...)[0, 0]; the native exp of an abstract matrix is
        # not implemented, so the dimension must be explicit.
        raise NotImplementedError('multivariate distributions of symbolic dimension (MatrixSymbol parameters) are not supported')
    return result


def MultivariateNormal(name,mu,sigma):
    mu,cov = _vector(mu),positive_matrix(sigma)
    if mu.rows != cov.rows:raise ValueError('location and covariance dimensions differ')
    k = int(mu.rows);joint = JointVariable(name,'MultivariateNormal',(mu,cov),k)
    return joint.initialize([continuous_model('Normal',(mu[i],S.sqrt(cov[i,i]))) for i in range(k)])


def MultivariateT(syms,mu,sigma,v):
    mu,cov,nu = _vector(mu),positive_matrix(sigma),scalar(v);positive(nu)
    if mu.rows != cov.rows:raise ValueError('location and shape dimensions differ')
    joint = JointVariable(syms,'MultivariateT',(mu,cov,nu),int(mu.rows));models = []
    for i in range(joint.dimension):
        base = continuous_model('StudentT',(nu,));location,scale = mu[i],S.sqrt(cov[i,i]);model = Distribution(location,scale,nu);model.family = 'ScaledStudentT'
        model._pdf = lambda x,b=base,a=location,c=scale:b.pdf((x-a)/c)/c
        model._cdf = lambda x,b=base,a=location,c=scale:b.cdf((x-a)/c)
        model._moment = lambda n,i=i:joint.joint_moment({joint.components[i]:int(n)})
        models.append(model)
    return joint.initialize(models)


def MultivariateLaplace(name,mu,sigma):
    mu,cov = _vector(mu),positive_matrix(sigma)
    if mu.rows != cov.rows:raise ValueError('location and shape dimensions differ')
    joint = JointVariable(name,'MultivariateLaplace',(mu,cov),int(mu.rows));models = []
    for i in range(joint.dimension):
        component = Distribution(mu[i],cov[i,i]);component.family = 'AsymmetricLaplaceMarginal'
        component._moment = lambda n,i=i:joint.joint_moment({joint.components[i]:int(n)})
        a,b = mu[i],cov[i,i]
        component._pdf = lambda x,a=a,b=b:S.exp(a*x/b-S.Abs(x)*S.sqrt(a*a+2*b)/b)/S.sqrt(a*a+2*b)
        models.append(component)
    return joint.initialize(models)


def Dirichlet(syms,*alpha):
    if len(alpha) == 1 and isinstance(alpha[0],(list,tuple,S.Matrix)):alpha = tuple(alpha[0])
    alpha = tuple(map(scalar,alpha))
    if len(alpha)<2:raise ValueError('Dirichlet requires at least two categories')
    positive(*alpha);total = S.Add(*alpha)
    joint = JointVariable(syms,'Dirichlet',(alpha,),len(alpha))
    joint.support = S.ProductSet(*([S.Interval(0,1)]*len(alpha)))
    return joint.initialize([continuous_model('Beta',(a,total-a)) for a in alpha])


MultivariateBeta = Dirichlet


def Multinomial(syms,n,*p):
    if len(p) == 1 and isinstance(p[0],(list,tuple)):p = tuple(p[0])
    n = scalar(n);check(n.is_integer is not False,'trial count must be integral');check(n>=0,'trial count must be nonnegative')
    probabilities = tuple(map(scalar,p))
    if len(probabilities)<2:raise ValueError('multinomial needs at least two categories')
    for q in probabilities:check(S.And(q>=0,q<=1),'probability must lie in [0,1]')
    check(S.simplify(S.Add(*probabilities)-1).is_zero is not False,'probabilities must sum to one')
    joint = JointVariable(syms,'Multinomial',(n,probabilities),len(probabilities),'discrete')
    joint.support = S.ProductSet(*([S.Range(0,n+1)]*len(probabilities)))
    return joint.initialize([finite_model('Binomial',(n,q,1,0)) for q in probabilities])


def NegativeMultinomial(syms,k0,*p):
    if len(p) == 1 and isinstance(p[0],(list,tuple)):p = tuple(p[0])
    r = scalar(k0);positive(r);probabilities = tuple(map(scalar,p))
    if not probabilities:raise ValueError('negative multinomial needs categories')
    for q in probabilities:check(S.And(q>=0,q<1),'probability must lie in [0,1)')
    p0 = 1-S.Add(*probabilities);positive(p0)
    joint = JointVariable(syms,'NegativeMultinomial',(r,probabilities),len(probabilities),'discrete')
    return joint.initialize([discrete_model('NegativeBinomial',(r,p0/(p0+q))) for q in probabilities])


def MultivariateEwens(syms,n,theta):
    n,theta = scalar(n),scalar(theta);positive(theta)
    if not n.is_Integer or n<1:raise ValueError('Ewens sample size must be a positive explicit integer')
    if n>30:raise ValueError('Ewens partition enumeration is limited to sample size 30')
    joint = JointVariable(syms,'MultivariateEwens',(n,theta),int(n),'discrete');models = []
    for i in range(int(n)):
        model = Distribution(n,theta);model.family = 'EwensMarginal';model.kind = 'finite';model._table = {}
        models.append(model)
    joint.initialize(models)
    for values,weight in joint.finite_worlds(joint.components):
        for rv in joint.components:
            value = values[rv._atom];table = rv.distribution._table;table[value] = table.get(value,S.Integer(0))+weight
    return joint


def NormalGamma(sym,mu,lamda,alpha,beta):
    mu,lamda,alpha,beta = map(scalar,(mu,lamda,alpha,beta));positive(lamda,alpha,beta)
    joint = JointVariable(sym,'NormalGamma',(mu,lamda,alpha,beta),2)
    student = continuous_model('StudentT',(2*alpha,));scale = S.sqrt(beta/(lamda*alpha))
    xmodel = Distribution(mu,lamda,alpha,beta);xmodel.family = 'NormalGammaLocation'
    xmodel._pdf = lambda x:student.pdf((x-mu)/scale)/scale
    xmodel._cdf = lambda x:student.cdf((x-mu)/scale)
    xmodel._moment = lambda n:joint.joint_moment({joint.components[0]:int(n)})
    return joint.initialize([xmodel,continuous_model('Gamma',(alpha,1/beta))])


def GeneralizedMultivariateLogGamma(syms,delta,v,lamda,mu):
    delta,v = scalar(delta),scalar(v);scales,rates = tuple(map(scalar,lamda)),tuple(map(scalar,mu))
    if len(scales)<2 or len(scales) != len(rates):raise ValueError('scale/rate vectors must have equal length, at least two')
    check(S.And(delta>0,delta<=1),'delta must lie in (0,1]');positive(v,*scales,*rates)
    joint = JointVariable(syms,'GeneralizedMultivariateLogGamma',(delta,v,scales,rates),len(scales));models = []
    for scale,rate in zip(scales,rates):
        model = Distribution(delta,v,scale,rate);model.family = 'LogGammaMarginal'
        model._pdf = lambda x,l=scale,r=rate:r*(delta/l)**v*S.exp(v*r*x-delta*S.exp(r*x)/l)/S.gamma(v)
        model._moment = lambda n,l=scale,r=rate:(S.digamma(v)+S.log(l/delta))/r if n == 1 else ((S.digamma(v)+S.log(l/delta))**2+S.polygamma(1,v))/r**2 if n == 2 else None
        models.append(model)
    return joint.initialize(models)


def GeneralizedMultivariateLogGammaOmega(syms,omega,v,lamda,mu):
    omega = positive_matrix(omega)
    if omega.rows != len(lamda) or omega.cols != len(mu):raise ValueError('omega and vector dimensions differ')
    for i in range(omega.rows):
        check(omega[i,i] == 1,'omega diagonal must equal one')
        for j in range(omega.cols):check(S.And(omega[i,j]>=0,omega[i,j]<=1),'omega entries must lie in [0,1]')
    delta = omega.det()**S.Rational(1,len(lamda)-1)
    return GeneralizedMultivariateLogGamma(syms,delta,v,lamda,mu)


def JointRV(symbol,pdf,_set=None):
    expression = scalar(pdf)
    if getattr(expression,'is_Lambda',False):variables = expression.variables;body = expression.expr
    else:
        from .._public_matrix import _indexed
        references = [(atom,entry.indices[0]) for atom,entry in _indexed.items() if str(entry.base) == str(symbol) and atom in expression.free_symbols and len(entry.indices) == 1]
        if not references:raise ValueError('joint density requires indexed coordinates or an explicit Lambda')
        variables = tuple(atom for atom,index in sorted(references,key=lambda row:int(row[1])));body = expression
    joint = JointVariable(symbol,'JointHandmade',(),len(variables));joint._lambda = S.Lambda(tuple(variables),body)
    if _set is not None:joint.support = _set
    models = []
    for i in range(joint.dimension):
        model = Distribution();model.family = 'HandmadeJointMarginal'
        def pdf_marginal(x,i=i):
            coords = [S.Dummy('marginal_%s' % j,real=True) for j in range(joint.dimension)];coords[i] = x
            value = joint.pdf(*coords)
            for j,var in enumerate(coords):
                if j != i:value = integral(value,var,-S.oo,S.oo)
            return value
        model._pdf = pdf_marginal;models.append(model)
    return joint.initialize(models)


def marginal_distribution(rv,*indices):
    if not getattr(rv,'_stats_joint',False):raise TypeError('a joint random variable is required')
    selected = [index.index if isinstance(index,RandomSymbol) else operator.index(index) for index in indices]
    if not selected or len(set(selected)) != len(selected) or any(not 0<=i<rv.dimension for i in selected):raise ValueError('select distinct valid component indices')
    args = [S.Dummy('marginal_x%s' % i,real=True) for i in selected]
    if len(selected) == 1:return S.Lambda(args[0],rv.components[selected[0]].distribution.pdf(args[0]))
    if len(selected) == rv.dimension:return S.Lambda(tuple(args),rv.pdf(*(args[selected.index(i)] for i in range(rv.dimension))))
    f,p = rv.family,rv.params
    if f in ('MultivariateNormal','MultivariateT','MultivariateLaplace'):
        mu,cov = p[:2];new_mu = [mu[i] for i in selected];new_cov = [[cov[i,j] for j in selected] for i in selected]
        factory = {'MultivariateNormal':MultivariateNormal,'MultivariateT':MultivariateT,'MultivariateLaplace':MultivariateLaplace}[f]
        result = factory(rv.name+'_marginal',new_mu,new_cov,*p[2:])
        return S.Lambda(tuple(args),result.pdf(*args))
    if f == 'Dirichlet':
        alpha, = p;rest = S.Add(*(a for i,a in enumerate(alpha) if i not in selected))
        body = S.gamma(S.Add(*alpha))/S.gamma(rest)*(1-S.Add(*args))**(rest-1)*S.Mul(*(x**(alpha[i]-1)/S.gamma(alpha[i]) for x,i in zip(args,selected)))
        return S.Lambda(tuple(args),body)
    if f == 'Multinomial':
        n,probabilities = p;remaining = n-S.Add(*args);p0 = 1-S.Add(*(probabilities[i] for i in selected))
        body = S.factorial(n)/S.factorial(remaining)*p0**remaining*S.Mul(*(probabilities[i]**x/S.factorial(x) for x,i in zip(args,selected)))
        return S.Lambda(tuple(args),S.Piecewise((body,remaining>=0),(0,True)))
    if f == 'NegativeMultinomial':
        r,probabilities = p;denominator = 1-S.Add(*(q for i,q in enumerate(probabilities) if i not in selected))
        model = NegativeMultinomial(rv.name+'_marginal',r,*[probabilities[i]/denominator for i in selected])
        return S.Lambda(tuple(args),model.pdf(*args))
    coordinates = [S.Dummy('integrated_%s' % i,real=True) for i in range(rv.dimension)]
    for i,x in zip(selected,args):coordinates[i] = x
    value = rv.pdf(*coordinates)
    for i,x in enumerate(coordinates):
        if i not in selected:value = integral(value,x,rv.components[i].distribution.lower,rv.components[i].distribution.upper)
    return S.Lambda(tuple(args),value)


__all__ = ['JointRV','Dirichlet','GeneralizedMultivariateLogGamma','GeneralizedMultivariateLogGammaOmega',
 'Multinomial','MultivariateBeta','MultivariateEwens','MultivariateT','NegativeMultinomial','NormalGamma',
 'MultivariateNormal','MultivariateLaplace','marginal_distribution']
