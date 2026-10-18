"""Typed probability spaces connected to the native expectation engine."""
from collections import defaultdict
import symox as S
from ._base import Distribution,RandomSymbol,PSpace,variable,scalar
from ._domains import *


class _DistributionType(type):
    def __instancecheck__(cls,obj):return isinstance(obj,Distribution) and obj.kind==cls._distribution_kind
class ContinuousDistribution(Distribution,metaclass=_DistributionType):_distribution_kind=kind='continuous'
class SingleContinuousDistribution(ContinuousDistribution):
    def compute_cdf(self,**kwargs):
        x=S.Dummy('x',real=True);return S.Lambda(x,self.cdf(x))
    def compute_characteristic_function(self,**kwargs):
        t=S.Dummy('t',real=True);return S.Lambda(t,self.mgf(S.I*t))
    def compute_moment_generating_function(self,**kwargs):
        t=S.Dummy('t',real=True);return S.Lambda(t,self.mgf(t))
    def compute_quantile(self,**kwargs):
        p=S.Dummy('p',real=True);return S.Lambda(p,self.quantile(p))
class SingleFiniteDistribution(Distribution,metaclass=_DistributionType):_distribution_kind=kind='finite'


class FiniteDensity(dict):
    def __call__(self,item):return self.get(scalar(item),S.Integer(0))
    dict=property(lambda self:dict(self))


class _SpaceQueries:
    def _ordinary(self,expr):
        return ordinary(expr)
    def _random(self,expr):
        return scalar(expr).xreplace({v.symbol:v._atom for v in self.variables})
    def where(self,condition):
        domain_cls=ConditionalFiniteDomain if self.domain.is_Finite else ConditionalContinuousDomain
        return domain_cls(self.domain,self._ordinary(condition))
    def compute_density(self,expr,**kwargs):
        from ._queries import density
        return density(self._random(expr),**kwargs)
    def compute_cdf(self,expr,**kwargs):
        from ._queries import cdf
        return cdf(self._random(expr),**kwargs)
    def compute_quantile(self,expr,**kwargs):
        from ._queries import quantile
        return quantile(self._random(expr),**kwargs)
    def compute_characteristic_function(self,expr,**kwargs):
        from ._queries import characteristic_function
        return characteristic_function(self._random(expr),**kwargs)
    def compute_moment_generating_function(self,expr,**kwargs):
        from ._queries import moment_generating_function
        return moment_generating_function(self._random(expr),**kwargs)


class SinglePSpace(_SpaceQueries,PSpace):
    def __init__(self,symbol,distribution):self._initialize(variable(str(symbol),distribution))
    @classmethod
    def _from_variable(cls,rv):
        instance=object.__new__(cls);instance._initialize(rv);return instance
    def _initialize(self,rv):
        PSpace.__init__(self,(rv,))
        if rv.distribution.kind=='continuous':self.domain=SingleContinuousDomain(rv.symbol,rv.distribution.set)
        elif rv.distribution.kind=='finite':self.domain=SingleFiniteDomain(rv.symbol,rv.distribution.set)
    value=property(lambda self:self.variables[0])
    symbols=property(lambda self:self.domain.symbols)
    pdf=property(lambda self:self.distribution.pdf(self.symbol))
    density=property(lambda self:self.distribution)
    def compute_expectation(self,expr,rvs=None,evaluate=False,**kwargs):
        from ._queries import E
        return E(self._random(expr),**kwargs)
    def probability(self,condition,**kwargs):
        from ._queries import P
        return P(self._random(condition),**kwargs)


class ProductPSpace(_SpaceQueries,PSpace):
    def __init__(self,*spaces):
        self.spaces=tuple(spaces)
        PSpace.__init__(self,[rv for space in spaces for rv in space.variables])
        domains=[space.domain for space in spaces]
        cls=ProductFiniteDomain if all(d.is_Finite for d in domains) else ProductContinuousDomain if all(d.is_Continuous for d in domains) else ProductDomain
        self.domain=cls(*domains)
    symbols=property(lambda self:self.domain.symbols)
    pdf=property(lambda self:S.Mul(*(space.pdf for space in self.spaces)))
    def compute_expectation(self,expr,rvs=None,evaluate=False,**kwargs):
        from ._queries import E
        return E(self._random(expr),**kwargs)
    def probability(self,condition,**kwargs):
        from ._queries import P
        return P(self._random(condition),**kwargs)
class IndependentProductPSpace(ProductPSpace):
    def __init__(self,*spaces):
        if len(spaces)==1 and isinstance(spaces[0],(list,tuple,set,frozenset)):spaces=tuple(spaces[0])
        groups=[{v.group for v in space.variables} for space in spaces]
        if any(a&b for i,a in enumerate(groups) for b in groups[i+1:]):raise ValueError('product spaces must be independent')
        super().__init__(*spaces)


class ContinuousPSpace(_SpaceQueries,PSpace):
    is_Continuous=True
    def __init__(self,domain,density):
        self.domain=domain;symbols=sorted(domain.symbols,key=str)
        expression=density(*symbols) if callable(density) else scalar(density)
        self.density=S.Lambda(tuple(symbols),expression)
        variables=[]
        for symbol in symbols:
            marginal=domain.compute_expectation(expression,set(symbols)-{symbol}).doit() if len(symbols)>1 else expression
            support=next(d.set for d in getattr(domain,'domains',(domain,)) if symbol in d.symbols)
            model=Distribution();model._support=support
            if isinstance(support,S.Interval):model.lower,model.upper=support.start,support.end
            model._pdf=lambda z,p=marginal,x=symbol:p.subs(x,z)
            variables.append(RandomSymbol(str(symbol),model,group=self))
        PSpace.__init__(self,variables);self.domain=domain
    symbols=property(lambda self:self.domain.symbols)
    pdf=property(lambda self:self.density.expr)
    def compute_expectation(self,expr,rvs=None,evaluate=False,**kwargs):
        expression=self.domain.compute_expectation(self.pdf*self._ordinary(expr))
        return expression.doit() if evaluate else expression
    def expectation(self,expr,variables):return self.compute_expectation(expr,evaluate=True)
    def joint_moment(self,factors):return self.compute_expectation(S.Mul(*(rv.symbol**power for rv,power in factors.items())),evaluate=True)
    def probability(self,condition,**kwargs):
        domain=self.where(condition)
        return domain.compute_expectation(self.pdf).doit()
    def sample_mapping(self,rng,variables):
        values={};mapping={};domains=getattr(self.domain,'domains',(self.domain,))
        ordered=sorted(self.variables,key=lambda rv:str(rv.symbol))
        for position,rv in enumerate(ordered):
            remaining={v.symbol for v in ordered[position+1:]}
            marginal=self.domain.compute_expectation(self.pdf,remaining).doit().subs(values,simultaneous=True)
            support=next(d.set for d in domains if rv.symbol in d.symbols)
            if not isinstance(support,S.Interval):raise NotImplementedError('conditional inverse sampling requires interval support')
            normalization=S.integrate(marginal,(rv.symbol,support.start,support.end))
            model=Distribution();model.lower,model.upper=support.start,support.end
            model._pdf=lambda z,e=marginal/normalization,x=rv.symbol:e.subs(x,z)
            draw=model.sample(rng,());values[rv.symbol]=S.Float(draw);mapping[rv._atom]=draw
        return {rv._atom:mapping[rv._atom] for rv in variables}


class SingleContinuousPSpace(SinglePSpace,ContinuousPSpace):is_Continuous=True


class FinitePSpace(_SpaceQueries,PSpace):
    is_Finite=True
    def __init__(self,domain,density):
        symbols=sorted(domain.symbols,key=str);density=dict(density);worlds=[]
        for state in domain:
            values=dict(state)
            key=values[symbols[0]] if len(symbols)==1 else frozenset(state)
            weight=scalar(density.get(frozenset(state),density.get(key,0)))
            worlds.append((values,weight))
        self._worlds=worlds;variables=[]
        for symbol in symbols:
            table=defaultdict(lambda:S.Integer(0))
            for state,weight in worlds:table[state[symbol]]+=weight
            model=Distribution();model.kind='finite';model._table=dict(table)
            variables.append(RandomSymbol(str(symbol),model,group=self))
        PSpace.__init__(self,variables);self.domain=domain
    symbols=property(lambda self:self.domain.symbols)
    def finite_worlds(self,variables):
        return [({rv._atom:state[rv.symbol] for rv in variables},weight) for state,weight in self._worlds]
    def sample_mapping(self,rng,variables):
        import numpy as np
        probabilities=np.asarray([float(p) for state,p in self._worlds])
        state=self._worlds[int(rng.choice(len(self._worlds),p=probabilities))][0]
        return {rv._atom:state[rv.symbol] for rv in variables}
    def prob_of(self,state):
        return sum((p for values,p in self._worlds if values==dict(state)),S.Integer(0))
    def compute_expectation(self,expr,rvs=None,evaluate=False,**kwargs):
        expr=self._ordinary(expr)
        return S.simplify(sum((p*expr.subs(state,simultaneous=True) for state,p in self._worlds),S.Integer(0)))
    def probability(self,condition,**kwargs):
        condition=self._ordinary(condition)
        return S.simplify(sum((p for state,p in self._worlds if bool(condition.subs(state,simultaneous=True))),S.Integer(0)))


class SingleFinitePSpace(SinglePSpace,FinitePSpace):is_Finite=True
class ProductFinitePSpace(IndependentProductPSpace,FinitePSpace):is_Finite=True


def make_space(variables):
    variables=tuple(variables)
    if len(variables)==1:
        rv=variables[0]
        cls=SingleContinuousPSpace if rv.distribution.kind=='continuous' else SingleFinitePSpace if rv.distribution.kind=='finite' else SinglePSpace
        return cls._from_variable(rv)
    if len({rv.group for rv in variables})==len(variables):
        cls=ProductFinitePSpace if all(v.distribution.kind=='finite' for v in variables) else IndependentProductPSpace
        return cls(*(make_space((rv,)) for rv in variables))
    return PSpace(variables)
