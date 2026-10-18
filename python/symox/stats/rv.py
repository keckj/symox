"""Probability spaces, random scalars and exact scalar queries."""
import symox as S
from ._base import Distribution,RandomSymbol,PSpace,Domain,scalar
from ._queries import (P,E,cdf,characteristic_function,density,dependent,given,
    independent,moment_generating_function,pspace,quantile,random_symbols,
    sample,sample_iter,sampling_density,where)
from ._processes import sample_stochastic_process

probability = P
expectation = E


def is_random(expr):return bool(random_symbols(expr))


def pspace_independent(a,b):
    return not {v.group for v in random_symbols(a)} & {v.group for v in random_symbols(b)}


def rv_subs(expr,symbols=None):
    variables = random_symbols(expr) if symbols is None else symbols
    return scalar(expr).subs({v._atom:v.symbol for v in variables},simultaneous=True)


def rs_swap(a,b):
    replacements = {}
    for source in a:
        matches = [target for target in b if source.symbol == target.symbol]
        if not matches:raise IndexError('no random variable with the same underlying symbol')
        replacements[source] = matches[0]
    return replacements


def sumsets(sets):return frozenset().union(*sets)


def sampling_E(expr,condition=None,numsamples=1,**kwargs):return E(expr,condition,numsamples=numsamples,**kwargs)
def sampling_P(expr,condition=None,numsamples=1,**kwargs):return P(expr,condition,numsamples=numsamples,**kwargs)


__all__ = ['Distribution','RandomSymbol','PSpace','cdf','characteristic_function','density',
 'dependent','expectation','given','independent','is_random','moment_generating_function',
 'probability','pspace','pspace_independent','quantile','random_symbols','rs_swap','rv_subs',
 'sample','sample_iter','sample_stochastic_process','sampling_E','sampling_P','sampling_density','sumsets','where']

from ._domains import RandomDomain,SingleDomain,MatrixDomain,ProductDomain,ConditionalDomain
from ._spaces import SinglePSpace,ProductPSpace,IndependentProductPSpace

class NamedArgsMixin:
    _argnames=()
    def __getattr__(self,name):
        if name in self._argnames:return self.args[self._argnames.index(name)]
        raise AttributeError(name)

class Density:
    def __init__(self,expr,condition=None):self.expr,self.condition=expr,condition;self.args=(expr,condition)
    def doit(self,**kwargs):return density(self.expr,self.condition,**kwargs)


def sample_iter_lambdify(expr,condition=None,size=(),numsamples=S.oo,seed=None,**kwargs):
    return sample_iter(expr,condition,size=size,numsamples=numsamples,seed=seed,**kwargs)
def sample_iter_subs(expr,condition=None,size=(),numsamples=S.oo,seed=None,**kwargs):
    return sample_iter(expr,condition,size=size,numsamples=numsamples,seed=seed,**kwargs)

__all__ += ['RandomDomain','SingleDomain','MatrixDomain','ProductDomain','ConditionalDomain','SinglePSpace',
 'ProductPSpace','IndependentProductPSpace','NamedArgsMixin','Density','sample_iter_lambdify','sample_iter_subs']
