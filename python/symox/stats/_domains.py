"""Symbolic probability domains, separate from their probability measures."""
from itertools import product
import symox as S


def ordinary(expression):
    from ._base import REGISTRY,scalar
    expression=scalar(expression)
    return expression.xreplace({atom:rv.symbol for atom,rv in REGISTRY.items() if atom in expression.free_symbols})


class RandomDomain:
    is_ProductDomain=is_Finite=is_Continuous=is_Discrete=False
    def __init__(self,symbols,set):self.symbols=frozenset(symbols);self.set=set
    def as_boolean(self):
        symbols=tuple(self.symbols)
        return S.Contains(symbols[0],self.set) if len(symbols)==1 else S.Contains(S.Tuple(*symbols),self.set)
    def __contains__(self,state):
        state=dict(state)
        if set(state)!=set(self.symbols):return False
        values=[state[s] for s in self.symbols]
        return bool(self.set.contains(values[0] if len(values)==1 else S.Tuple(*values)))


class SingleDomain(RandomDomain):
    def __init__(self,symbol,set):
        self.symbol=S.Symbol(symbol) if isinstance(symbol,str) else symbol
        self.symbols=frozenset([self.symbol]);self.set=set;self.args=(self.symbol,set)
    def as_boolean(self):return self.set.as_relational(self.symbol) if hasattr(self.set,'as_relational') else S.Contains(self.symbol,self.set)

class MatrixDomain(SingleDomain):pass

class ProductDomain(RandomDomain):
    is_ProductDomain=True
    def __init__(self,*domains):
        flattened=[]
        for domain in domains:flattened.extend(domain.domains if isinstance(domain,ProductDomain) else [domain])
        if sum(len(d.symbols) for d in flattened)!=len(set().union(*(set(d.symbols) for d in flattened))):raise ValueError('product domain symbols must be disjoint')
        self.domains=tuple(flattened);self.args=self.domains
        self.symbols=frozenset().union(*(d.symbols for d in self.domains))
        self._support=S.ProductSet(*(d.set for d in self.domains))
    set=property(lambda self:self._support)
    def as_boolean(self):return S.And(*(d.as_boolean() for d in self.domains))
    def __contains__(self,state):
        state=dict(state)
        return set(state)==set(self.symbols) and all({s:state[s] for s in d.symbols} in d for d in self.domains)


class ConditionalDomain(RandomDomain):
    def __init__(self,fulldomain,condition):
        self.fulldomain,self.condition=fulldomain,ordinary(condition)
        self.symbols=fulldomain.symbols;self.args=(fulldomain,self.condition)
    def as_boolean(self):return S.And(self.fulldomain.as_boolean(),self.condition)
    def __contains__(self,state):return state in self.fulldomain and bool(self.condition.subs(dict(state),simultaneous=True))
    @property
    def set(self):
        if len(self.symbols)!=1:raise NotImplementedError('a multivariate conditional support is represented by its domain predicate')
        symbol=next(iter(self.symbols))
        return S.ConditionSet(symbol,self.condition,self.fulldomain.set)


class ContinuousDomain(RandomDomain):
    is_Continuous=True
    def compute_expectation(self,expr,variables=None,**kwargs):
        raise NotImplementedError('continuous integration requires concrete single or product domains')

class SingleContinuousDomain(SingleDomain,ContinuousDomain):
    def compute_expectation(self,expr,variables=None,**kwargs):
        selected=self.symbols if variables is None else frozenset(variables)
        if not selected:return expr
        if selected!=self.symbols:raise ValueError('integration variables must belong to the domain')
        if isinstance(self.set,S.Interval):return S.Integral(expr,(self.symbol,self.set.start,self.set.end))
        if isinstance(self.set,S.Union):return S.Add(*(SingleContinuousDomain(self.symbol,piece).compute_expectation(expr) for piece in self.set.args))
        raise NotImplementedError('continuous integration requires interval support')

class ProductContinuousDomain(ProductDomain,ContinuousDomain):
    def compute_expectation(self,expr,variables=None,**kwargs):
        selected=self.symbols if variables is None else frozenset(variables)
        for domain in self.domains:
            subset=domain.symbols&selected
            if subset:expr=domain.compute_expectation(expr,subset,**kwargs)
        return expr

class ConditionalContinuousDomain(ConditionalDomain,ContinuousDomain):
    @property
    def set(self):
        if len(self.symbols)!=1:return super().set
        symbol=next(iter(self.symbols))
        return S.Intersection(self.fulldomain.set,S.solve_univariate_inequality(self.condition,symbol,relational=False))
    def compute_expectation(self,expr,variables=None,**kwargs):
        if len(self.symbols)==1:
            symbol=next(iter(self.symbols))
            return SingleContinuousDomain(symbol,self.set).compute_expectation(expr,variables)
        return self.fulldomain.compute_expectation(S.Piecewise((expr,self.condition),(0,True)),variables)


class FiniteDomain(RandomDomain):
    is_Finite=is_Discrete=True
    def __init__(self,elements):
        self._elements=tuple(frozenset(state) for state in elements)
        self.symbols=frozenset(s for state in self._elements for s,v in state);self.args=(self._elements,)
    def __iter__(self):return iter(self._elements)
    @property
    def elements(self):return frozenset(self)
    @property
    def dict(self):return tuple(dict(state) for state in self)
    def __contains__(self,state):return frozenset(dict(state).items()) in self.elements
    def as_boolean(self):return S.Or(*(S.And(*(S.Eq(s,v) for s,v in state)) for state in self))
    @property
    def set(self):
        symbols=sorted(self.symbols,key=str)
        if len(symbols)==1:return S.FiniteSet(*(dict(state)[symbols[0]] for state in self))
        return S.FiniteSet(*(S.Tuple(*(dict(state)[s] for s in symbols)) for state in self))

class SingleFiniteDomain(FiniteDomain,SingleDomain):
    def __init__(self,symbol,set):
        self.symbol=S.Symbol(symbol) if isinstance(symbol,str) else symbol
        support=set if isinstance(set,S.Set) else S.FiniteSet(*set)
        if support.is_finite_set is not True:raise ValueError('finite domains require finite support')
        self._support=support;self.symbols=frozenset([self.symbol]);self.args=(self.symbol,support)
    set=property(lambda self:self._support)
    def __iter__(self):return (frozenset([(self.symbol,v)]) for v in self.set)

class ProductFiniteDomain(ProductDomain,FiniteDomain):
    def __init__(self,*domains):ProductDomain.__init__(self,*domains)
    def __iter__(self):
        count=1
        for domain in self.domains:count*=len(domain.elements)
        if count>1000000:raise ValueError('finite product domain exceeds one million states')
        return (frozenset().union(*states) for states in product(*(iter(domain) for domain in self.domains)))

class ConditionalFiniteDomain(ConditionalDomain,FiniteDomain):
    def __iter__(self):return (state for state in self.fulldomain if bool(self.condition.subs(dict(state),simultaneous=True)))
    def __contains__(self,state):return ConditionalDomain.__contains__(self,state)
    @property
    def set(self):return FiniteDomain.set.fget(self)
    def as_boolean(self):return FiniteDomain.as_boolean(self)
