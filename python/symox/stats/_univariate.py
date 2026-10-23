"""Scalar finite, countable and continuous probability laws."""
import inspect
import math
import numpy as np
import symox as S
from ._base import Distribution,HALF as H,scalar,positive,check,variable,sum_expr


def _integer(value,minimum=0):
    check(value.is_integer is not False,'parameter must be integral')
    check(value >= minimum,'integer parameter is outside its range')


def _probability(p,closed=True):
    check(S.And(p >= 0,p <= 1) if closed else S.And(p > 0,p < 1),'probability is outside its range')


def _shift(n,location,moment):
    if not n.is_Integer or n < 0:return None
    return S.Add(*(S.binomial(n,k)*location**(n-k)*moment(S.Integer(k)) for k in range(int(n)+1)))


def _normal_moment(n,mu,sigma):
    if not n.is_Integer or n < 0:return None
    return S.Add(*(S.binomial(n,2*j)*mu**(n-2*j)*sigma**(2*j)*S.factorial(2*j)/(2**j*S.factorial(j)) for j in range(int(n)//2+1)))


def _exists(value,condition,n):
    failure = S.oo
    return S.Piecewise((value,condition),(failure,True))


def _bounded(m,a,b):
    m.lower,m.upper = a,b
    check(b > a,'upper bound must exceed lower bound')


def _polynomial_moment(n,pieces):
    if not n.is_Integer or n < 0:return None
    x = S.Dummy('moment_x')
    total = S.Integer(0)
    for polynomial,a,b in pieces(x):
        if a == b:continue
        for (degree,),coefficient in S.Poly(S.expand(polynomial),x).terms():
            power = n+degree+1
            total += coefficient*(b**power-a**power)/power
    return S.simplify(total)


CONTINUOUS = {
 'Arcsin':('a b',(0,1)), 'Benini':('alpha beta sigma',()),
 'Beta':('alpha beta',()), 'BetaNoncentral':('alpha beta lamda',()),
 'BetaPrime':('alpha beta',()), 'BoundedPareto':('alpha left right',()),
 'Cauchy':('x0 gamma',()), 'Chi':('k',()), 'ChiNoncentral':('k l',()),
 'ChiSquared':('k',()), 'Dagum':('p a b',()), 'Davis':('b n mu',()),
 'Erlang':('k l',()), 'ExGaussian':('mean std rate',()),
 'Exponential':('rate',()), 'ExponentialPower':('mu alpha beta',()),
 'FDistribution':('d1 d2',()), 'FisherZ':('d1 d2',()),
 'Frechet':('a s m',(1,0)), 'Gamma':('k theta',()),
 'GammaInverse':('a b',()), 'Gompertz':('b eta',()),
 'Gumbel':('beta mu minimum',(False,)), 'Kumaraswamy':('a b',()),
 'Laplace':('mu b',()), 'Levy':('mu c',()), 'Logistic':('mu s',()),
 'LogCauchy':('mu sigma',()), 'LogLogistic':('alpha beta',()),
 'LogitNormal':('mu s',()), 'LogNormal':('mean std',()),
 'Lomax':('alpha lamda',()), 'Moyal':('mu sigma',()), 'Maxwell':('a',()),
 'Nakagami':('mu omega',()), 'Normal':('mean std',()),
 'GaussianInverse':('mean shape',()), 'Pareto':('xm alpha',()),
 'PowerFunction':('alpha a b',()), 'QuadraticU':('a b',()),
 'RaisedCosine':('mu s',()), 'Rayleigh':('sigma',()), 'Reciprocal':('a b',()),
 'StudentT':('nu',()), 'ShiftedGompertz':('b eta',()),
 'Trapezoidal':('a b c d',()), 'Triangular':('a b c',()),
 'Uniform':('left right',()), 'UniformSum':('n',()),
 'VonMises':('mu k',()), 'Wald':('mean shape',()),
 'Weibull':('alpha beta',()), 'WignerSemicircle':('R',()),
}


def continuous_model(family,params):
    values = tuple(map(scalar,params))
    m = Distribution(*values);m.family = family;m._argnames = tuple(CONTINUOUS[family][0].split())
    if family == 'Erlang':
        k,rate = values;_integer(k,1);positive(rate)
        return continuous_model('Gamma',(k,1/rate))
    if family == 'Wald':return continuous_model('GaussianInverse',values)
    # Formulae use native scalars throughout: fractional exponents stay exact.
    if family == 'Normal':
        mu,sigma = values;positive(sigma)
        m._pdf = lambda x:S.exp(-(x-mu)**2/(2*sigma**2))/(S.sqrt(2*S.pi)*sigma)
        m._cdf = lambda x:(1+S.erf((x-mu)/(S.sqrt(2)*sigma)))/2
        m._moment = lambda n:_normal_moment(n,mu,sigma)
        m._mgf = lambda t:S.exp(mu*t+sigma**2*t**2/2)
        m._quantile = lambda p:mu+sigma*S.sqrt(2)*S.erfinv(2*p-1)
        m._sampler = lambda r,z:r.normal(float(mu),float(sigma),z)
        m._entropy = (1+S.log(2*S.pi*sigma**2))/2
    elif family == 'Uniform':
        a,b = values;_bounded(m,a,b)
        m._pdf = lambda x:S.Piecewise((1/(b-a),S.And(x>=a,x<=b)),(0,True))
        m._cdf = lambda x:(x-a)/(b-a)
        m._cdf_form = lambda x:S.Piecewise((0,x < a),((x-a)/(b-a),x <= b),(1,True))
        m._moment = lambda n:(b**(n+1)-a**(n+1))/((n+1)*(b-a))
        m._mgf = lambda t:(S.exp(b*t)-S.exp(a*t))/(t*(b-a))
        m._quantile = lambda p:a+(b-a)*p
        m._sampler = lambda r,z:r.uniform(float(a),float(b),z)
        m._entropy = S.log(b-a)
    elif family in ('Gamma','ChiSquared','Exponential'):
        if family == 'Gamma':k,theta = values;positive(k,theta)
        elif family == 'ChiSquared':k0, = values;_integer(k0,1);k,theta = k0/2,S.Integer(2)
        else:rate, = values;positive(rate);k,theta = S.Integer(1),1/rate
        m.lower = S.Integer(0)
        m._pdf = lambda x:x**(k-1)*S.exp(-x/theta)/(S.gamma(k)*theta**k)
        m._cdf = lambda x:S.lowergamma(k,x/theta)/S.gamma(k)
        if family == 'Gamma':
            m._cdf_form = lambda x:S.Piecewise((S.lowergamma(k,x/theta)/S.gamma(k),x > 0),(0,True))
            # cdf(X, meijerg=True): the reference integrates the density.
            m._cdf_integrated = lambda x:S.Piecewise((k*S.lowergamma(k,x/theta)/S.gamma(k+1),x >= 0),(0,True))
        elif family == 'ChiSquared':
            m._cdf_form = lambda x:S.Piecewise((S.lowergamma(k,x/2)/S.gamma(k),x >= 0),(0,True))
        else:
            m._cdf_form = lambda x:S.Piecewise((1-S.exp(-rate*x),x >= 0),(0,True))
        m._moment = lambda n:theta**n*S.rf(k,n)
        m._mgf = lambda t:(1-theta*t)**(-k)
        m._sampler = lambda r,z:r.gamma(float(k),float(theta),z)
        m._entropy = k+S.log(theta*S.gamma(k))+(1-k)*S.digamma(k)
        if family == 'Exponential':m._quantile = lambda p:-theta*S.log(1-p)
    elif family in ('Beta','BetaPrime','BetaNoncentral'):
        a,b = values[:2];positive(a,b);m.lower = S.Integer(0)
        if family != 'BetaPrime':m.upper = S.Integer(1)
        if family == 'Beta':
            m._pdf = lambda x:x**(a-1)*(1-x)**(b-1)/S.beta(a,b)
            m._cdf = lambda x:S.betainc_regularized(a,b,0,x)
            m._moment = lambda n:S.rf(a,n)/S.rf(a+b,n)
            m._mgf = lambda t:S.hyper((a,),(a+b,),t)
            m._sampler = lambda r,z:r.beta(float(a),float(b),z)
            m._entropy = S.log(S.beta(a,b))-(a-1)*S.digamma(a)-(b-1)*S.digamma(b)+(a+b-2)*S.digamma(a+b)
        elif family == 'BetaPrime':
            m._pdf = lambda x:x**(a-1)*(1+x)**(-a-b)/S.beta(a,b)
            m._cdf = lambda x:S.betainc_regularized(a,b,0,x/(1+x))
            m._moment = lambda n:_exists(S.gamma(a+n)*S.gamma(b-n)/(S.gamma(a)*S.gamma(b)),b>n,n)
            m._sampler = lambda r,z:r.gamma(float(a),1,z)/r.gamma(float(b),1,z)
        else:
            nc = values[2];check(nc>=0,'noncentrality must be nonnegative')
            if nc == 0:return continuous_model('Beta',(a,b))
            def mixture(x,moment=None):
                j = S.Dummy('j',integer=True,nonnegative=True)
                weight = S.exp(-nc/2)*(nc/2)**j/S.factorial(j)
                body = (x**(a+j-1)*(1-x)**(b-1)/S.beta(a+j,b)) if moment is None else S.rf(a+j,moment)/S.rf(a+b+j,moment)
                return S.Sum(weight*body,(j,0,S.oo))
            m._pdf = mixture;m._moment = lambda n:mixture(None,n)
            m._sampler = lambda r,z:r.beta(float(a)+r.poisson(float(nc)/2,z),float(b),z)
    elif family in ('Cauchy','LogCauchy'):
        mu,scale = values;positive(scale)
        if family == 'Cauchy':
            m._pdf = lambda x:1/(S.pi*scale*(1+((x-mu)/scale)**2))
            m._cdf = lambda x:H+S.atan((x-mu)/scale)/S.pi
            m._quantile = lambda p:mu+scale*S.tan(S.pi*(p-H))
            m._moment = lambda n:S.oo if n.is_even else S.nan
            m._sampler = lambda r,z:float(mu)+float(scale)*r.standard_cauchy(z)
            m._entropy = S.log(4*S.pi*scale)
        else:
            m.lower = S.Integer(0)
            m._pdf = lambda x:scale/(S.pi*x*((S.log(x)-mu)**2+scale**2))
            m._cdf = lambda x:H+S.atan((S.log(x)-mu)/scale)/S.pi
            m._quantile = lambda p:S.exp(mu+scale*S.tan(S.pi*(p-H)))
            m._moment = lambda n:S.oo
            m._sampler = lambda r,z:np.exp(float(mu)+float(scale)*r.standard_cauchy(z))
    elif family in ('Chi','Maxwell','Rayleigh','Nakagami','ChiNoncentral'):
        m.lower = S.Integer(0)
        if family == 'Chi':k, = values;_integer(k,1);scale = S.Integer(1)
        elif family == 'Maxwell':scale, = values;positive(scale);k = S.Integer(3)
        elif family == 'Rayleigh':scale, = values;positive(scale);k = S.Integer(2)
        elif family == 'Nakagami':
            shape,omega = values;check(shape>=H,'Nakagami shape must be at least 1/2');positive(omega)
            k,scale = 2*shape,S.sqrt(omega/(2*shape))
        else:
            k,nc = values;_integer(k,1);check(nc>=0,'noncentrality must be nonnegative');scale = S.Integer(1)
        if family == 'ChiNoncentral':
            if nc == 0:return continuous_model('Chi',(k,))
            m._pdf = lambda x:S.exp(-(x*x+nc*nc)/2)*x**k*nc/(nc*x)**(k/2)*S.besseli(k/2-1,nc*x)
            m._moment = lambda n:2**(n/2)*S.gamma((k+n)/2)/S.gamma(k/2)*S.hyper((-n/2,),(k/2,),-nc**2/2)
            m._sampler = lambda r,z:np.sqrt(r.noncentral_chisquare(float(k),float(nc*nc),z))
        else:
            m._pdf = lambda x:2**(1-k/2)*x**(k-1)*S.exp(-x*x/(2*scale*scale))/(S.gamma(k/2)*scale**k)
            m._cdf = lambda x:S.lowergamma(k/2,x*x/(2*scale*scale))/S.gamma(k/2)
            m._moment = lambda n:scale**n*2**(n/2)*S.gamma((k+n)/2)/S.gamma(k/2)
            m._sampler = lambda r,z:float(scale)*np.sqrt(r.chisquare(float(k),z))
            if family == 'Rayleigh':m._quantile = lambda p:scale*S.sqrt(-2*S.log(1-p))
            if family == 'Nakagami':
                m._pdf = lambda x:2*shape**shape/(S.gamma(shape)*omega**shape)*x**(2*shape-1)*S.exp(-shape/omega*x**2)
                m._cdf_form = lambda x:S.Piecewise((S.lowergamma(shape,(shape/omega)*x**2)/S.gamma(shape),x > 0),(0,True))
    elif family in ('LogNormal','LogitNormal'):
        mu,sigma = values;positive(sigma)
        m.lower = S.Integer(0)
        if family == 'LogNormal':
            m._pdf = lambda x:S.exp(-(S.log(x)-mu)**2/(2*sigma*sigma))/(x*S.sqrt(2*S.pi)*sigma)
            m._cdf = lambda x:(1+S.erf((S.log(x)-mu)/(S.sqrt(2)*sigma)))/2
            m._moment = lambda n:S.exp(n*mu+n*n*sigma*sigma/2)
            m._quantile = lambda p:S.exp(mu+sigma*S.sqrt(2)*S.erfinv(2*p-1))
            m._sampler = lambda r,z:r.lognormal(float(mu),float(sigma),z)
            m._entropy = mu+H+S.log(sigma*S.sqrt(2*S.pi))
        else:
            m.upper = S.Integer(1)
            m._pdf = lambda x:S.exp(-(S.log(x/(1-x))-mu)**2/(2*sigma*sigma))/(x*(1-x)*S.sqrt(2*S.pi)*sigma)
            m._cdf = lambda x:(1+S.erf((S.log(x/(1-x))-mu)/(S.sqrt(2)*sigma)))/2
            m._quantile = lambda p:1/(1+S.exp(-mu-sigma*S.sqrt(2)*S.erfinv(2*p-1)))
            m._sampler = lambda r,z:1/(1+np.exp(-r.normal(float(mu),float(sigma),z)))
    elif family in ('Pareto','BoundedPareto','Lomax','LogLogistic','Frechet','GammaInverse','FDistribution','Dagum'):
        if family == 'Pareto':
            xm,alpha = values;positive(xm,alpha);m.lower = xm
            m._pdf = lambda x:alpha*xm**alpha/x**(alpha+1)
            m._cdf = lambda x:1-(xm/x)**alpha
            m._moment = lambda n:_exists(alpha*xm**n/(alpha-n),alpha>n,n)
            m._quantile = lambda p:xm/(1-p)**(1/alpha)
            m._sampler = lambda r,z:float(xm)*(1+r.pareto(float(alpha),z))
        elif family == 'BoundedPareto':
            alpha,left,right = values;positive(alpha,left);_bounded(m,left,right)
            m._pdf = lambda x:alpha*left**alpha*x**(-alpha-1)/(1-(left/right)**alpha)
            m._cdf = lambda x:(1-(left/x)**alpha)/(1-(left/right)**alpha)
            m._moment = lambda n:S.Piecewise((alpha*left**alpha*S.log(right/left)/(1-(left/right)**alpha),S.Eq(n,alpha)),(alpha*left**alpha*(right**(n-alpha)-left**(n-alpha))/((n-alpha)*(1-(left/right)**alpha)),True))
            m._quantile = lambda p:left/(1-p*(1-(left/right)**alpha))**(1/alpha)
        elif family == 'Lomax':
            alpha,scale = values;positive(alpha,scale);m.lower = S.Integer(0)
            m._pdf = lambda x:alpha/scale*(1+x/scale)**(-alpha-1)
            m._cdf = lambda x:1-(1+x/scale)**(-alpha)
            m._cdf_form = lambda x:S.Piecewise((1-(1+x/scale)**(-alpha),x >= 0),(0,True))
            m._moment = lambda n:_exists(scale**n*S.gamma(n+1)*S.gamma(alpha-n)/S.gamma(alpha),alpha>n,n)
            m._quantile = lambda p:scale*((1-p)**(-1/alpha)-1)
            m._sampler = lambda r,z:float(scale)*r.pareto(float(alpha),z)
        elif family == 'LogLogistic':
            alpha,beta = values;positive(alpha,beta);m.lower = S.Integer(0)
            m._pdf = lambda x:beta/alpha*(x/alpha)**(beta-1)/(1+(x/alpha)**beta)**2
            m._cdf = lambda x:1/(1+(x/alpha)**(-beta))
            m._moment = lambda n:_exists(alpha**n*S.gamma(1+n/beta)*S.gamma(1-n/beta),beta>n,n)
            m._quantile = lambda p:alpha*(p/(1-p))**(1/beta)
        elif family == 'Frechet':
            alpha,scale,location = values;positive(alpha,scale);m.lower = location
            m._pdf = lambda x:alpha/scale*((x-location)/scale)**(-1-alpha)*S.exp(-((x-location)/scale)**(-alpha))
            m._cdf = lambda x:S.exp(-((x-location)/scale)**(-alpha))
            m._cdf_form = lambda x:S.Piecewise((S.exp(-((x-location)/scale)**(-alpha)),x >= location),(0,True))
            m._moment = lambda n:_shift(n,location,lambda j:S.Integer(1) if j == 0 else _exists(scale**j*S.gamma(1-j/alpha),alpha>j,j))
            m._quantile = lambda p:location+scale*(-S.log(p))**(-1/alpha)
        elif family == 'GammaInverse':
            alpha,beta = values;positive(alpha,beta);m.lower = S.Integer(0)
            m._pdf = lambda x:beta**alpha*x**(-alpha-1)*S.exp(-beta/x)/S.gamma(alpha)
            m._cdf = lambda x:S.uppergamma(alpha,beta/x)/S.gamma(alpha)
            m._cdf_form = lambda x:S.Piecewise((S.uppergamma(alpha,beta/x)/S.gamma(alpha),x > 0),(0,True))
            m._moment = lambda n:_exists(beta**n*S.gamma(alpha-n)/S.gamma(alpha),alpha>n,n)
            m._sampler = lambda r,z:float(beta)/r.gamma(float(alpha),1,z)
        elif family == 'FDistribution':
            d1,d2 = values;_integer(d1,1);_integer(d2,1);m.lower = S.Integer(0)
            m._pdf = lambda x:S.sqrt((d1*x)**d1*d2**d2/(d1*x+d2)**(d1+d2))/(x*S.beta(d1/2,d2/2))
            m._cdf = lambda x:S.betainc_regularized(d1/2,d2/2,0,d1*x/(d1*x+d2))
            m._moment = lambda n:_exists((d2/d1)**n*S.gamma(d1/2+n)*S.gamma(d2/2-n)/(S.gamma(d1/2)*S.gamma(d2/2)),d2>2*n,n)
            m._sampler = lambda r,z:r.f(float(d1),float(d2),z)
        else:
            p,a,b = values;positive(p,a,b);m.lower = S.Integer(0)
            m._pdf = lambda x:a*p*(x/b)**(a*p)/(x*(1+(x/b)**a)**(p+1))
            m._cdf = lambda x:(1+(x/b)**(-a))**(-p)
            m._cdf_form = lambda x:S.Piecewise(((1+(x/b)**(-a))**(-p),x >= 0),(0,True))
            m._moment = lambda n:_exists(b**n*S.gamma(1-n/a)*S.gamma(p+n/a)/S.gamma(p),a>n,n)
            m._quantile = lambda q:b*(q**(-1/p)-1)**(-1/a)
    elif family == 'StudentT':
        nu, = values;positive(nu)
        m._pdf = lambda x:(1+x*x/nu)**(-(nu+1)/2)/(S.sqrt(nu)*S.beta(H,nu/2))
        m._cdf = lambda x:H+x*S.gamma((nu+1)/2)*S.hyper((H,(nu+1)/2),(S.Rational(3,2),),-x*x/nu)/(S.sqrt(S.pi*nu)*S.gamma(nu/2))
        m._moment = lambda n:S.Piecewise((0,nu>n),(S.nan,True)) if n.is_odd else _exists(nu**(n/2)*S.gamma((n+1)/2)*S.gamma((nu-n)/2)/(S.sqrt(S.pi)*S.gamma(nu/2)),nu>n,n)
        m._sampler = lambda r,z:r.standard_t(float(nu),z)
    elif family in ('Laplace','Logistic','ExponentialPower'):
        mu,scale = values[:2];positive(scale)
        if family == 'Laplace':
            m._pdf = lambda x:S.exp(-S.Abs(x-mu)/scale)/(2*scale)
            m._cdf = lambda x:S.Piecewise((S.exp((x-mu)/scale)/2,x<mu),(1-S.exp(-(x-mu)/scale)/2,True))
            m._moment = lambda n:_shift(n,mu,lambda j:S.factorial(j)*scale**j if j.is_even else S.Integer(0))
            m._mgf = lambda t:S.exp(mu*t)/(1-scale*scale*t*t)
            m._quantile = lambda p:S.Piecewise((mu+scale*S.log(2*p),p<H),(mu-scale*S.log(2*(1-p)),True))
            m._sampler = lambda r,z:r.laplace(float(mu),float(scale),z)
            m._entropy = 1+S.log(2*scale)
        elif family == 'Logistic':
            m._pdf = lambda x:S.exp(-(x-mu)/scale)/(scale*(1+S.exp(-(x-mu)/scale))**2)
            m._cdf = lambda x:1/(1+S.exp(-(x-mu)/scale))
            m._moment = lambda n:_shift(n,mu,lambda j:S.Integer(1) if j == 0 else (2*(1-2**(1-j))*S.zeta(j)*S.factorial(j)*scale**j if j.is_even else S.Integer(0)))
            m._mgf = lambda t:S.exp(mu*t)*S.beta(1-scale*t,1+scale*t)
            m._quantile = lambda p:mu+scale*S.log(p/(1-p))
            m._sampler = lambda r,z:r.logistic(float(mu),float(scale),z)
            m._entropy = S.log(scale)+2
        else:
            beta = values[2];positive(beta)
            m._pdf = lambda x:beta*S.exp(-(S.Abs(x-mu)/scale)**beta)/(2*scale*S.gamma(1/beta))
            m._cdf = lambda x:H+S.sign(x-mu)*S.lowergamma(1/beta,(S.Abs(x-mu)/scale)**beta)/(2*S.gamma(1/beta))
            m._moment = lambda n:_shift(n,mu,lambda j:scale**j*S.gamma((j+1)/beta)/S.gamma(1/beta) if j.is_even else S.Integer(0))
            m._sampler = lambda r,z:float(mu)+float(scale)*r.choice([-1,1],size=z)*r.gamma(1/float(beta),1,z)**(1/float(beta))
    elif family == 'ExGaussian':
        mu,sigma,rate = values;positive(sigma,rate)
        m._pdf = lambda x:rate*S.exp(rate*(2*mu+rate*sigma*sigma-2*x)/2)*S.erfc((mu+rate*sigma*sigma-x)/(S.sqrt(2)*sigma))/2
        normalcdf = lambda x:(1+S.erf(x/S.sqrt(2)))/2
        m._cdf = lambda x:normalcdf((x-mu)/sigma)-S.exp(-rate*(x-mu)+rate*rate*sigma*sigma/2)*normalcdf((x-mu)/sigma-rate*sigma)
        def exgaussian_cdf(x):
            # The reference composes two Gaussian CDFs in the scaled variable.
            u,v = rate*(x-mu),rate*sigma
            gauss = lambda value,mean,std:S.erf(S.sqrt(2)*(-mean+value)/(2*std))/2+H
            return gauss(u,0,v)-S.exp(-u+v**2/2+S.log(gauss(u,v**2,v)))
        m._cdf_form = exgaussian_cdf
        m._moment = lambda n:S.Add(*(S.binomial(n,j)*_normal_moment(n-j,mu,sigma)*S.factorial(j)/rate**j for j in range(int(n)+1))) if n.is_Integer else None
        m._mgf = lambda t:S.exp(mu*t+sigma*sigma*t*t/2)/(1-t/rate)
        m._sampler = lambda r,z:r.normal(float(mu),float(sigma),z)+r.exponential(1/float(rate),z)
    elif family == 'GaussianInverse':
        mu,shape = values;positive(mu,shape);m.lower = S.Integer(0)
        m._pdf = lambda x:S.sqrt(shape/(2*S.pi*x**3))*S.exp(-shape*(x-mu)**2/(2*x*mu*mu))
        phi = lambda x:(1+S.erf(x/S.sqrt(2)))/2
        m._cdf = lambda x:phi(S.sqrt(shape/x)*(x/mu-1))+S.exp(2*shape/mu)*phi(-S.sqrt(shape/x)*(x/mu+1))
        m._moment = lambda n:mu if n == 1 else mu*mu+mu**3/shape if n == 2 else None
        m._mgf = lambda t:S.exp(shape/mu*(1-S.sqrt(1-2*mu*mu*t/shape)))
        m._sampler = lambda r,z:r.wald(float(mu),float(shape),z)
    elif family == 'Weibull':
        alpha,beta = values;positive(alpha,beta);m.lower = S.Integer(0)
        m._pdf = lambda x:beta/alpha*(x/alpha)**(beta-1)*S.exp(-(x/alpha)**beta)
        m._cdf = lambda x:1-S.exp(-(x/alpha)**beta)
        m._moment = lambda n:alpha**n*S.gamma(1+n/beta)
        m._quantile = lambda p:alpha*(-S.log(1-p))**(1/beta)
        m._sampler = lambda r,z:float(alpha)*r.weibull(float(beta),z)
    elif family == 'Arcsin':
        a,b = values;_bounded(m,a,b)
        m._pdf = lambda x:1/(S.pi*S.sqrt((x-a)*(b-x)))
        m._cdf = lambda x:2*S.asin(S.sqrt((x-a)/(b-a)))/S.pi
        m._cdf_form = lambda x:S.Piecewise((0,x < a),(2*S.asin(S.sqrt((x-a)/(b-a)))/S.pi,x <= b),(1,True))
        m._moment = lambda n:_shift(n,a,lambda j:(b-a)**j*S.rf(H,j)/S.factorial(j))
        m._quantile = lambda p:a+(b-a)*S.sin(S.pi*p/2)**2
        m._sampler = lambda r,z:float(a)+float(b-a)*r.beta(.5,.5,z)
    elif family == 'Kumaraswamy':
        a,b = values;positive(a,b);m.lower,m.upper = S.Integer(0),S.Integer(1)
        m._pdf = lambda x:a*b*x**(a-1)*(1-x**a)**(b-1)
        m._cdf = lambda x:1-(1-x**a)**b
        m._cdf_form = lambda x:S.Piecewise((0,x < 0),(1-(1-x**a)**b,x <= 1),(1,True))
        m._moment = lambda n:b*S.beta(1+n/a,b)
        m._quantile = lambda p:(1-(1-p)**(1/b))**(1/a)
    elif family == 'PowerFunction':
        alpha,a,b = values;positive(alpha);_bounded(m,a,b)
        m._pdf = lambda x:alpha*(x-a)**(alpha-1)/(b-a)**alpha
        m._cdf = lambda x:((x-a)/(b-a))**alpha
        m._moment = lambda n:_shift(n,a,lambda j:(b-a)**j*alpha/(alpha+j))
        m._quantile = lambda p:a+(b-a)*p**(1/alpha)
    elif family == 'Reciprocal':
        a,b = values;positive(a);_bounded(m,a,b)
        m._pdf = lambda x:1/(x*S.log(b/a))
        m._cdf = lambda x:S.log(x/a)/S.log(b/a)
        m._moment = lambda n:(b**n-a**n)/(n*S.log(b/a))
        m._quantile = lambda p:a*(b/a)**p
    elif family == 'Benini':
        a,b,scale = values;positive(a,b,scale);m.lower = scale
        m._pdf = lambda x:S.exp(-a*S.log(x/scale)-b*S.log(x/scale)**2)*(a/x+2*b*S.log(x/scale)/x)
        m._cdf = lambda x:1-S.exp(-a*S.log(x/scale)-b*S.log(x/scale)**2)
        m._cdf_form = lambda x:S.Piecewise((1-S.exp(-a*S.log(x/scale)-b*S.log(x/scale)**2),x >= scale),(0,True))
        m._moment = lambda n:scale**n*(1+n*S.sqrt(S.pi)/(2*S.sqrt(b))*S.exp((a-n)**2/(4*b))*S.erfc((a-n)/(2*S.sqrt(b))))
        m._quantile = lambda p:scale*S.exp((S.sqrt(a*a-4*b*S.log(1-p))-a)/(2*b))
    elif family == 'Davis':
        b,n,mu = values;positive(b,mu);check(n>1,'Davis shape must exceed one');m.lower = mu
        m._pdf = lambda x:b**n*(x-mu)**(-n-1)/((S.exp(b/(x-mu))-1)*S.gamma(n)*S.zeta(n))
        m._moment = lambda r:_shift(r,mu,lambda j:S.Integer(1) if j == 0 else _exists(b**j*S.gamma(n-j)*S.zeta(n-j)/(S.gamma(n)*S.zeta(n)),n>j+1,j))
    elif family in ('Gompertz','ShiftedGompertz'):
        b,eta = values;positive(b,eta);m.lower = S.Integer(0)
        if family == 'Gompertz':
            m._pdf = lambda x:b*eta*S.exp(b*x)*S.exp(eta)*S.exp(-eta*S.exp(b*x))
            m._cdf = lambda x:1-S.exp(eta-eta*S.exp(b*x))
            m._mgf = lambda t:S.exp(eta)*eta**(-t/b)*S.uppergamma(1+t/b,eta)
            m._moment = lambda n:S.exp(eta)*S.expint(1,eta)/b if n == 1 else None
            m._quantile = lambda p:S.log(1-S.log(1-p)/eta)/b
        else:
            m._pdf = lambda x:b*S.exp(-b*x)*S.exp(-eta*S.exp(-b*x))*(1+eta*(1-S.exp(-b*x)))
            m._cdf = lambda x:(1-S.exp(-b*x))*S.exp(-eta*S.exp(-b*x))
    elif family == 'Gumbel':
        beta,mu,minimum = values;positive(beta);sign = -1 if bool(minimum) else 1
        z = lambda x:sign*(x-mu)/beta
        m._pdf = lambda x:S.exp(-z(x)-S.exp(-z(x)))/beta
        m._cdf = lambda x:S.exp(-S.exp(-z(x))) if sign == 1 else 1-S.exp(-S.exp(-z(x)))
        mean = mu+sign*beta*S.EulerGamma
        m._moment = lambda n:mean if n == 1 else mean*mean+beta*beta*S.pi**2/6 if n == 2 else None
        m._mgf = lambda t:S.exp(mu*t)*S.gamma(1-sign*beta*t)
        m._quantile = lambda p:mu-sign*beta*S.log(-S.log(p if sign == 1 else 1-p))
        m._sampler = lambda r,z:float(mu)+sign*r.gumbel(0,float(beta),z)
    elif family == 'Levy':
        mu,c = values;positive(c);m.lower = mu
        m._pdf = lambda x:S.sqrt(c/(2*S.pi))*S.exp(-c/(2*(x-mu)))/(x-mu)**S.Rational(3,2)
        m._cdf = lambda x:S.erfc(S.sqrt(c/(2*(x-mu))))
        m._moment = lambda n:S.oo
        m._quantile = lambda p:mu+c/(2*S.erfcinv(p)**2)
        m._sampler = lambda r,z:float(mu)+float(c)/r.standard_normal(z)**2
    elif family == 'Moyal':
        mu,sigma = values;positive(sigma)
        m._pdf = lambda x:S.exp(-(S.exp(-(x-mu)/sigma)+(x-mu)/sigma)/2)/(S.sqrt(2*S.pi)*sigma)
        m._cdf = lambda x:S.erfc(S.exp(-(x-mu)/(2*sigma))/S.sqrt(2))
        mean = mu+sigma*(S.EulerGamma+S.log(2))
        m._moment = lambda n:mean if n == 1 else mean*mean+sigma*sigma*S.pi**2/2 if n == 2 else None
        m._mgf = lambda t:S.exp(mu*t)*2**(-sigma*t)*S.gamma(H-sigma*t)/S.sqrt(S.pi)
        m._quantile = lambda p:mu-2*sigma*S.log(S.sqrt(2)*S.erfcinv(p))
        m._sampler = lambda r,z:float(mu)-float(sigma)*np.log(r.chisquare(1,z))
    elif family == 'FisherZ':
        d1,d2 = values;positive(d1,d2)
        m._pdf = lambda x:2*d1**(d1/2)*d2**(d2/2)*S.exp(d1*x)/(S.beta(d1/2,d2/2)*(d1*S.exp(2*x)+d2)**((d1+d2)/2))
        m._cdf = lambda x:S.betainc_regularized(d1/2,d2/2,0,d1*S.exp(2*x)/(d2+d1*S.exp(2*x)))
        m._moment = lambda n:(S.digamma(d1/2)-S.digamma(d2/2)+S.log(d2/d1))/2 if n == 1 else None
        m._sampler = lambda r,z:np.log(r.f(float(d1),float(d2),z))/2
    elif family in ('Triangular','Trapezoidal','QuadraticU'):
        if family == 'Triangular':
            a,b,c = values;_bounded(m,a,b);check(S.And(c>=a,c<=b),'triangular mode must lie in the interval')
            def pieces(x):
                result = []
                if c != a:result.append((2*(x-a)/((b-a)*(c-a)),a,c))
                if c != b:result.append((2*(b-x)/((b-a)*(b-c)),c,b))
                return result
            m._pdf = lambda x:S.Piecewise(*((p,S.And(x>=l,x<u)) for p,l,u in pieces(x)),(2/(b-a),S.Eq(x,c)),(0,True))
            m._cdf = lambda x:S.Piecewise(((x-a)**2/((b-a)*(c-a)),x<c),(1-(b-x)**2/((b-a)*(b-c)),True)) if c != a and c != b else ((x-a)**2/(b-a)**2 if c == b else 1-(b-x)**2/(b-a)**2)
            m._quantile = lambda p:S.Piecewise((a+S.sqrt(p*(b-a)*(c-a)),p<(c-a)/(b-a)),(b-S.sqrt((1-p)*(b-a)*(b-c)),True))
            m._sampler = lambda r,z:r.triangular(float(a),float(c),float(b),z)
        elif family == 'Trapezoidal':
            a,b,c,d = values;_bounded(m,a,d);check(S.And(a<=b,b<c,c<=d),'trapezoid parameters must be ordered')
            height = 2/(d+c-a-b)
            def pieces(x):
                result = []
                if b != a:result.append((height*(x-a)/(b-a),a,b))
                result.append((height,b,c))
                if d != c:result.append((height*(d-x)/(d-c),c,d))
                return result
            def trapezoid_pdf(x):
                # Reference pieces: half-open [a, b), [b, c), closed [c, d].
                result = []
                if b != a:result.append((2*(x-a)/((b-a)*(d+c-a-b)),S.And(a <= x,x < b)))
                result.append((2/(d+c-a-b),S.And(b <= x,x < c)))
                result.append((2*(d-x)/((d-c)*(d+c-a-b)),S.And(c <= x,x <= d)) if d != c else (height,S.Eq(x,c)))
                return S.Piecewise(*result,(0,True))
            m._pdf = trapezoid_pdf
            rising = lambda x:S.Integer(0) if b == a else height*(x-a)**2/(2*(b-a))
            falling = lambda x:S.Integer(1) if d == c else 1-height*(d-x)**2/(2*(d-c))
            m._cdf = lambda x:S.Piecewise((rising(x),x<b),(height*((b-a)/2+x-b),x<=c),(falling(x),True))
            left_mass,right_start = height*(b-a)/2,1-height*(d-c)/2
            m._quantile = lambda p:S.Piecewise((a+S.sqrt(2*p*(b-a)/height),p<left_mass),(b+p/height-(b-a)/2,p<=right_start),(d-S.sqrt(2*(1-p)*(d-c)/height),True))
        else:
            a,b = values;_bounded(m,a,b)
            def pieces(x):return [(12*(x-(a+b)/2)**2/(b-a)**3,a,b)]
            m._pdf = lambda x:S.Piecewise((pieces(x)[0][0],S.And(x>=a,x<=b)),(0,True))
            m._cdf = lambda x:H+4*(x-(a+b)/2)**3/(b-a)**3
        m._moment = lambda n:_polynomial_moment(n,pieces)
    elif family == 'RaisedCosine':
        mu,scale = values;positive(scale);m.lower,m.upper = mu-scale,mu+scale
        m._pdf = lambda x:S.Piecewise(((1+S.cos(S.pi*(x-mu)/scale))/(2*scale),S.And(x>=mu-scale,x<=mu+scale)),(0,True))
        m._cdf = lambda x:H+(x-mu)/(2*scale)+S.sin(S.pi*(x-mu)/scale)/(2*S.pi)
        var = scale*scale*(S.Rational(1,3)-2/S.pi**2)
        m._moment = lambda n:mu if n == 1 else mu*mu+var if n == 2 else None
    elif family == 'UniformSum':
        n, = values;_integer(n,1);m.lower,m.upper = S.Integer(0),n
        def irwin(x,power):
            k = S.Dummy('k',integer=True)
            return sum_expr((-1)**k*S.binomial(n,k)*(x-k)**power,k,S.Integer(0),S.floor(x))/S.factorial(power)
        m._pdf = lambda x:irwin(x,n-1)
        m._cdf = lambda x:irwin(x,n)
        def irwin_cdf(x):
            k = S.Dummy('k')
            return S.Piecewise((0,x < 0),(1/S.factorial(n)*S.Sum((-1)**k*S.binomial(n,k)*(x-k)**n,(k,0,S.floor(x))),x <= n),(1,True))
        m._cdf_form = irwin_cdf
        m._moment = lambda k:S.stirling(n+k,n,kind=2)/S.binomial(n+k,n) if n.is_Integer and k.is_Integer else None
        m._mgf = lambda t:((S.exp(t)-1)/t)**n
        m._sampler = lambda r,z:r.uniform(size=tuple(z)+(int(n),)).sum(axis=-1)
    elif family == 'VonMises':
        mu,k = values;positive(k);m.lower,m.upper = S.Integer(0),2*S.pi
        m._pdf = lambda x:S.exp(k*S.cos(x-mu))/(2*S.pi*S.besseli(0,k))
        m._sampler = lambda r,z:np.mod(r.vonmises(float(mu),float(k),z),2*np.pi)
    elif family == 'WignerSemicircle':
        radius, = values;positive(radius);m.lower,m.upper = -radius,radius
        m._pdf = lambda x:2*S.sqrt(radius*radius-x*x)/(S.pi*radius*radius)
        m._cdf = lambda x:H+(x*S.sqrt(radius*radius-x*x)+radius*radius*S.asin(x/radius))/(S.pi*radius*radius)
        m._moment = lambda n:S.Integer(0) if n.is_odd else S.catalan(n/2)*(radius/2)**n
        m._sampler = lambda r,z:float(radius)*(2*r.beta(1.5,1.5,z)-1)
    return m

FINITE = {
 'DiscreteUniform':('items',()), 'Die':('sides',(6,)),
 'Bernoulli':('p succ fail',(1,0)), 'Coin':('p',(H,)),
 'Binomial':('n p succ fail',(1,0)), 'BetaBinomial':('n alpha beta',()),
 'Hypergeometric':('N m n',()), 'Rademacher':('',()),
 'IdealSoliton':('k',()), 'RobustSoliton':('k delta c',()),
}
DISCRETE = {'FlorySchulz':('a',()),'Geometric':('p',()),'Hermite':('a1 a2',()),
            'Logarithmic':('p',()),'NegativeBinomial':('r p',()),'Poisson':('lamda',()),
            'Skellam':('mu1 mu2',()),'YuleSimon':('rho',()),'Zeta':('s',())}


def finite_model(family,params):
    if family == 'DiscreteUniform':
        items = tuple(map(scalar,params[0]))
        if not items or len(set(items)) != len(items):raise ValueError('uniform items must be a nonempty set of distinct outcomes')
        m = Distribution();m.args = items;m.family = family;m.kind = 'finite'
        m._table = {value:S.Rational(1,len(items)) for value in items}
        return m
    if family == 'Coin':return finite_model('Bernoulli',(params[0],S.Symbol('H'),S.Symbol('T')))
    values = tuple(map(scalar,params))
    m = Distribution(*values);m.family = family;m.kind = 'finite';m._argnames = tuple(FINITE[family][0].split())
    if family == 'Die':
        sides, = values;_integer(sides,1);m.lower,m.upper = S.Integer(1),sides
        m._pdf = lambda k:S.Piecewise((1/sides,S.And(k>=1,k<=sides,S.Contains(k,S.Integers))),(0,True))
        if sides.is_Integer:m._table = {S.Integer(k):1/sides for k in range(1,int(sides)+1)}
        else:m._held_density = True;m._symbolic_mass = lambda k:1/sides  # SymPy: density(Dn).dict is Density(DieDistribution(n))
        m._moment = lambda n:(S.bernoulli(n+1,sides+1)-S.bernoulli(n+1))/(sides*(n+1)) if n.is_Integer else None
        m._mgf = lambda t:S.exp(t)*(1-S.exp(sides*t))/(sides*(1-S.exp(t)))
        m._sampler = lambda r,z:r.integers(1,int(sides)+1,z)
    elif family == 'Bernoulli':
        p,success,failure = values;_probability(p)
        table = {failure:1-p};table[success] = table.get(success,S.Integer(0))+p
        m._table = table
        m._moment = lambda n:p*success**n+(1-p)*failure**n
        m._mgf = lambda t:p*S.exp(t*success)+(1-p)*S.exp(t*failure)
    elif family == 'Binomial':
        n,p,success,failure = values;_integer(n);_probability(p)
        m.lower,m.upper = S.Integer(0),n
        mass = lambda k:S.binomial(n,k)*p**k*(1-p)**(n-k)
        if n.is_Integer:
            table = {}
            for j in range(int(n)+1):
                value = j*success+(n-j)*failure
                table[value] = table.get(value,S.Integer(0))+mass(S.Integer(j))
            m._table = table
        else:
            m._held_density = True  # SymPy: a held Density(BinomialDistribution(...))
            if success == 1 and failure == 0:m._symbolic_mass = mass
        if success == failure:
            m._pdf = lambda k:S.Piecewise((1,S.Eq(k,n*failure)),(0,True))
        else:
            def probability(value):
                k = (value-n*failure)/(success-failure)
                return S.Piecewise((mass(k),S.And(k>=0,k<=n,S.Contains(k,S.Integers))),(0,True))
            m._pdf = probability
        m.lower,m.upper = S.Min(n*success,n*failure),S.Max(n*success,n*failure)
        def raw(k):
            if not k.is_Integer:return None
            base = lambda q:S.Add(*(S.stirling(q,j,kind=2)*S.ff(n,j)*p**j for j in range(int(q)+1)))
            return _shift(k,n*failure,lambda q:(success-failure)**q*base(q))
        m._moment = raw
        m._mgf = lambda t:((1-p)*S.exp(t*failure)+p*S.exp(t*success))**n
        m._sampler = lambda r,z:float(failure)*int(n)+(float(success)-float(failure))*r.binomial(int(n),float(p),z)
    elif family == 'BetaBinomial':
        n,a,b = values;_integer(n);positive(a,b);m.lower,m.upper = S.Integer(0),n
        m._pdf = lambda k:S.binomial(n,k)*S.beta(k+a,n-k+b)/S.beta(a,b)
        if n.is_Integer:m._table = {S.Integer(k):m._pdf(S.Integer(k)) for k in range(int(n)+1)}
        m._moment = lambda k:S.Add(*(S.stirling(k,j,kind=2)*S.ff(n,j)*S.rf(a,j)/S.rf(a+b,j) for j in range(int(k)+1))) if k.is_Integer else None
        m._sampler = lambda r,z:r.binomial(int(n),r.beta(float(a),float(b),z),z)
    elif family == 'Hypergeometric':
        N,M,n = values;_integer(N);_integer(M);_integer(n);check(S.And(M<=N,n<=N),'sample and marked population must not exceed population')
        m.lower,m.upper = S.Max(0,n+M-N),S.Min(n,M)
        m._pdf = lambda k:S.binomial(M,k)*S.binomial(N-M,n-k)/S.binomial(N,n)
        if m.lower.is_Integer and m.upper.is_Integer:m._table = {S.Integer(k):m._pdf(S.Integer(k)) for k in range(int(m.lower),int(m.upper)+1)}
        m._moment = lambda k:S.Add(*(S.stirling(k,j,kind=2)*S.ff(n,j)*S.ff(M,j)/S.ff(N,j) for j in range(min(int(k),int(n) if n.is_Integer else int(k))+1))) if k.is_Integer else None
        m._sampler = lambda r,z:r.hypergeometric(int(M),int(N-M),int(n),z)
    elif family == 'Rademacher':
        m._table = {S.Integer(-1):H,S.Integer(1):H}
        m._moment = lambda n:S.Integer(0) if n.is_odd else S.Integer(1)
        m._mgf = lambda t:S.cosh(t)
        m._sampler = lambda r,z:2*r.integers(0,2,z)-1
    elif family in ('IdealSoliton','RobustSoliton'):
        k = values[0];_integer(k,1);m.lower,m.upper = S.Integer(1),k
        if family == 'RobustSoliton':
            delta,c = values[1:];check(S.And(delta>0,delta<=1),'delta must lie in (0,1]');positive(c)
        if not k.is_Integer or (family == 'RobustSoliton' and (c*S.log(k/delta)).free_symbols):
            # SymPy keeps the density of a symbolic soliton as a held
            # Density object that evaluates once the parameters are numbers.
            m._held_density = True
            if family == 'IdealSoliton':
                def held_pmf(x):
                    if x.is_integer is False:return S.Integer(0)
                    return S.Piecewise((1/k,S.Eq(x,1)),(1/(x*(x-1)),S.And(x >= 1,x <= k)),(0,True))
                m._pdf = held_pmf
            else:
                def held_pmf(x):raise NotImplementedError('the robust soliton mass function requires numeric parameters')
                m._pdf = held_pmf
            return m
        table = {S.Integer(1):1/k}
        table.update({S.Integer(j):S.Rational(1,j*(j-1)) for j in range(2,int(k)+1)})
        if family == 'RobustSoliton':
            R = c*S.log(k/delta)*S.sqrt(k)
            if R.free_symbols:raise NotImplementedError('robust soliton requires numeric cutoff parameters')
            cutoff = int(round(float((k/R).evalf())))
            for j in range(1,int(k)+1):
                tau = R/(k*j) if j < cutoff else R*S.log(R/delta)/k if j == cutoff else S.Integer(0)
                table[S.Integer(j)] += tau
            normalizer = S.Add(*table.values())
            table = {v:p/normalizer for v,p in table.items()}
            m.R,m.Z = R,normalizer
        m._table = table
    return m


def discrete_model(family,params):
    values = tuple(map(scalar,params));m = Distribution(*values);m.family = family;m.kind = 'discrete'
    m._argnames = tuple(DISCRETE[family][0].split());m.lower = S.Integer(0)
    if family == 'Poisson':
        lamda, = values;check(lamda>=0,'Poisson rate must be nonnegative')
        m._pdf = lambda k:S.exp(-lamda)*lamda**k/S.factorial(k)
        m._cdf = lambda x:S.uppergamma(S.floor(x)+1,lamda)/S.gamma(S.floor(x)+1)
        m._moment = lambda n:S.Add(*(S.stirling(n,j,kind=2)*lamda**j for j in range(int(n)+1))) if n.is_Integer else None
        m._mgf = lambda t:S.exp(lamda*(S.exp(t)-1))
        m._sampler = lambda r,z:r.poisson(float(lamda),z)
    elif family in ('Geometric','NegativeBinomial','FlorySchulz'):
        if family == 'Geometric':p, = values;shape,offset = S.Integer(1),S.Integer(1);check(p>0,'geometric probability must be positive');_probability(p)
        elif family == 'FlorySchulz':p, = values;_probability(p,False);shape,offset = S.Integer(2),S.Integer(1)
        else:shape,p = values;positive(shape);_probability(p,False);offset = S.Integer(0)
        m.lower = offset
        m._pdf = lambda k:S.binomial(k-offset+shape-1,k-offset)*(1-p)**(k-offset)*p**shape
        if family == 'Geometric':m._pdf = lambda k:(1-p)**(k-1)*p
        elif family == 'FlorySchulz':m._pdf = lambda k:p**2*k*(1-p)**(k-1)
        m._cdf = lambda x:S.betainc_regularized(shape,S.floor(x)-offset+1,0,p)
        m._moment = lambda n:_shift(n,offset,lambda k:S.Add(*(S.stirling(k,j,kind=2)*S.rf(shape,j)*((1-p)/p)**j for j in range(int(k)+1))))
        m._mgf = lambda t:S.exp(offset*t)*(p/(1-(1-p)*S.exp(t)))**shape
        m._sampler = lambda r,z:r.negative_binomial(float(shape),float(p),z)+int(offset)
    elif family == 'Hermite':
        a,b = values;check(S.And(a>=0,b>=0),'Hermite rates must be nonnegative')
        def mass(k):
            j = S.Dummy('j',integer=True)
            return S.exp(-a-b)*sum_expr(a**(k-2*j)*b**j/(S.factorial(k-2*j)*S.factorial(j)),j,S.Integer(0),S.floor(k/2))
        m._pdf = mass
        m._mgf = lambda t:S.exp(a*(S.exp(t)-1)+b*(S.exp(2*t)-1))
        m._moment = lambda n:(a+2*b) if n == 1 else (a+2*b)**2+a+4*b if n == 2 else None
        m._sampler = lambda r,z:r.poisson(float(a),z)+2*r.poisson(float(b),z)
    elif family == 'Logarithmic':
        p, = values;_probability(p,False);m.lower = S.Integer(1)
        m._pdf = lambda k:-p**k/(k*S.log(1-p))
        m._moment = lambda n:S.polylog(1-n,p)/(-S.log(1-p))
        m._mgf = lambda t:S.log(1-p*S.exp(t))/S.log(1-p)
        m._sampler = lambda r,z:r.logseries(float(p),z)
    elif family == 'Skellam':
        a,b = values;check(S.And(a>=0,b>=0),'Skellam rates must be nonnegative');m.lower = -S.oo;m._support = S.Integers
        if a == 0 and b == 0:
            m.lower = m.upper = S.Integer(0);m._support = S.FiniteSet(0)
            m._pdf = lambda k:S.Piecewise((1,S.Eq(k,0)),(0,True))
        elif b == 0:
            m.lower = S.Integer(0);m._support = S.Range(0,S.oo)
            m._pdf = lambda k:S.exp(-a)*a**k/S.factorial(k)
        elif a == 0:
            m.upper = S.Integer(0);m._support = S.Range(-S.oo,1)
            m._pdf = lambda k:S.exp(-b)*b**(-k)/S.factorial(-k)
        else:m._pdf = lambda k:S.exp(-a-b)*(a/b)**(k/2)*S.besseli(k,2*S.sqrt(a*b))
        m._moment = lambda n:a-b if n == 1 else (a-b)**2+a+b if n == 2 else None
        m._mgf = lambda t:S.exp(a*(S.exp(t)-1)+b*(S.exp(-t)-1))
        m._sampler = lambda r,z:r.poisson(float(a),z)-r.poisson(float(b),z)
    elif family == 'YuleSimon':
        rho, = values;positive(rho);m.lower = S.Integer(1)
        m._pdf = lambda k:rho*S.beta(k,rho+1)
        m._cdf = lambda x:1-S.floor(x)*S.beta(S.floor(x),rho+1)
        m._moment = lambda n:_exists(rho/(rho-1),rho>1,n) if n == 1 else _exists(rho*rho/((rho-1)*(rho-2)),rho>2,n) if n == 2 else None
        m._mgf = lambda t:rho*S.exp(t)*S.hyper((1,1),(rho+2,),S.exp(t))/(rho+1)
        m._sampler = lambda r,z:r.geometric(r.beta(float(rho),1,z),z)
    elif family == 'Zeta':
        exponent, = values;check(exponent>1,'zeta exponent must exceed one');m.lower = S.Integer(1)
        m._pdf = lambda k:1/(k**exponent*S.zeta(exponent))
        m._moment = lambda n:_exists(S.zeta(exponent-n)/S.zeta(exponent),exponent>n+1,n)
        m._mgf = lambda t:S.polylog(exponent,S.exp(t))/S.zeta(exponent)
        m._sampler = lambda r,z:r.zipf(float(exponent),z)
    return m


def _constructor(family,spec,kind):
    names,defaults = spec;names = names.split();mandatory = len(names)-len(defaults)
    parameters = [inspect.Parameter('name',inspect.Parameter.POSITIONAL_OR_KEYWORD)]
    parameters += [inspect.Parameter(name,inspect.Parameter.POSITIONAL_OR_KEYWORD,
                   default=inspect.Parameter.empty if i < mandatory else defaults[i-mandatory]) for i,name in enumerate(names)]
    signature = inspect.Signature(parameters)
    def create(*args,**kwargs):
        bound = signature.bind(*args,**kwargs);bound.apply_defaults()
        name = bound.arguments['name'];params = [bound.arguments[p] for p in names]
        if family == 'Normal' and isinstance(params[0],(list,tuple,S.Matrix)):
            from ._joint import MultivariateNormal
            return MultivariateNormal(name,*params)
        if kind == 'continuous':model = continuous_model(family,params)
        elif kind == 'finite':model = finite_model(family,params)
        else:model = discrete_model(family,params)
        model._params = tuple(params)
        model._remake = lambda values:create(name,*values).distribution
        return variable(name,model)
    create.__name__ = family;create.__qualname__ = family;create.__signature__ = signature
    create.__doc__ = 'Construct a %s random variable with native symbolic parameters.' % family
    return create


for _kind,_specs in [('continuous',CONTINUOUS),('finite',FINITE),('discrete',DISCRETE)]:
    for _family,_spec in _specs.items():globals()[_family] = _constructor(_family,_spec,_kind)


def _laplace(scalar_laplace):
    from functools import wraps
    @wraps(scalar_laplace)
    def Laplace(name,mu,b):
        # SymPy: list or matrix location and scale give the multivariate law.
        if isinstance(mu,(list,S.MatrixBase)) and isinstance(b,(list,S.MatrixBase)):
            from ._joint import MultivariateLaplace
            return MultivariateLaplace(name,mu,b)
        return scalar_laplace(name,mu,b)
    return Laplace
Laplace = _laplace(Laplace)


class FiniteDistributionHandmade(Distribution):
    kind = 'finite';family = 'FiniteHandmade'
    def __init__(self,density):
        super().__init__()
        self._table = {scalar(value):scalar(probability) for value,probability in dict(density).items()}
        if not self._table:raise ValueError('a finite distribution needs at least one outcome')
        for p in self._table.values():check(S.And(p>=0,p<=1),'probabilities must lie in [0,1]')
        check(S.simplify(S.Add(*self._table.values())-1).is_zero is not False,'probabilities must sum to one')


def FiniteRV(name,density,**kwargs):
    if kwargs.pop('check',True) is False:
        model = Distribution();model.kind = 'finite';model.family = 'FiniteHandmade';model._table = {scalar(k):scalar(v) for k,v in dict(density).items()}
    else:model = FiniteDistributionHandmade(density)
    if kwargs:raise TypeError('unknown finite random variable options')
    return variable(name,model)


class ContinuousDistributionHandmade(Distribution):
    def __init__(self,pdf,set=None):
        if set is None:set = S.Interval(-S.oo,S.oo)
        if not isinstance(set,S.Interval):raise NotImplementedError('handmade continuous support must be an interval')
        super().__init__();self.family = 'ContinuousHandmade';self._support = set;self.lower,self.upper = set.start,set.end
        if getattr(pdf,'is_Lambda',False):
            x = pdf.variables[0];expression = pdf.expr
        else:
            expression = scalar(pdf);free = expression.free_symbols
            if len(free) != 1:raise ValueError('a handmade density needs a univariate Lambda')
            x = next(iter(free))
        self.args = expression,x,self.lower,self.upper
        self._bound_symbols = {x}
        self._pdf = lambda z:expression.subs(x,z)


def _integer_run(support):
    """A finite support of consecutive integers (given as a Python set or a
    ``FiniteSet``) as the equal ``Range``; any other support is returned
    unchanged."""
    import builtins
    if isinstance(support,(builtins.set,frozenset,list,tuple)):support = S.FiniteSet(*support)
    if getattr(support,'is_FiniteSet',False) and len(support) > 0:
        values = list(support)
        if all(getattr(v,'is_Integer',False) for v in values):
            low,high = min(values),max(values)
            if high-low+1 == len(values):return S.Range(low,high+1)
    return support


class DiscreteDistributionHandmade(Distribution):
    kind = 'discrete'
    def __init__(self,pdf,set=None):
        if set is None:set = S.Integers
        set = _integer_run(set)
        super().__init__();self.family = 'DiscreteHandmade';self._support = set
        if isinstance(set,S.Range):self.lower,self.upper = set.start,set.stop-1
        elif set == S.Naturals:self.lower,self.upper = S.Integer(1),S.oo
        elif set == S.Naturals0:self.lower,self.upper = S.Integer(0),S.oo
        elif set == S.Integers:self.lower,self.upper = -S.oo,S.oo
        else:raise NotImplementedError('handmade discrete support must be an integer range')
        if getattr(pdf,'is_Lambda',False):x,expression = pdf.variables[0],pdf.expr
        else:
            expression = scalar(pdf);free = expression.free_symbols
            if len(free) != 1:raise ValueError('a handmade mass function needs a univariate Lambda')
            x = next(iter(free))
        self.args = expression,x,self.lower,self.upper
        self._bound_symbols = {x}
        self._pdf = lambda z:expression.subs(x,z)


def ContinuousRV(symbol,density,set=None,**kwargs):
    x = S.Symbol(symbol) if isinstance(symbol,str) else scalar(symbol)
    if kwargs:raise TypeError('unknown continuous random variable options')
    return variable(str(x),ContinuousDistributionHandmade(S.Lambda(x,density),set))


def DiscreteRV(symbol,density,set=None,**kwargs):
    x = S.Symbol(symbol,integer=True) if isinstance(symbol,str) else scalar(symbol)
    if kwargs:raise TypeError('unknown discrete random variable options')
    return variable(str(x),DiscreteDistributionHandmade(S.Lambda(x,density),set))


__all__ = list(CONTINUOUS)+list(FINITE)+list(DISCRETE)+['FiniteRV','ContinuousRV','DiscreteRV',
    'FiniteDistributionHandmade','ContinuousDistributionHandmade','DiscreteDistributionHandmade']
