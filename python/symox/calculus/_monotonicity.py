"""Monotonicity proofs on connected real intervals using native sets."""
from .. import S
from ._analysis import scalar, variable, real_domain, inequality, zeros, real_constants


def _subset(domain, allowed):
    missing=domain-allowed
    if missing.is_empty is not None:return missing.is_empty
    return domain.is_subset(allowed)


def _context(expression, interval, symbol):
    from .. import Dummy, Interval
    from .util import continuous_domain
    expression=scalar(expression)
    if symbol is None:
        free=expression.free_symbols
        if len(free)>1:raise NotImplementedError('Supply the variable for a multivariate monotonicity check')
        symbol=next(iter(free)) if free else Dummy('constant_x',real=True)
    symbol=variable(symbol)
    interval=real_domain(interval)
    if symbol in interval.free_symbols:
        raise ValueError('Monotonicity interval bounds must be independent of the variable')
    if interval==S.EmptySet:return True,None,None,None
    singleton=interval.is_FiniteSet and len(interval.args)==1
    if interval!=S.Reals and not isinstance(interval,Interval) and not singleton:
        raise NotImplementedError('Monotonicity currently requires a connected real interval')
    parameters=expression.free_symbols-{symbol}
    if any(getattr(parameter,'is_real',None) is not True for parameter in parameters):
        raise NotImplementedError('Monotonicity parameters must be known to be real')
    real=Dummy('monotonic_x',real=True)
    expression=expression.xreplace({symbol:real})
    real_constants(expression,real)
    continuous=continuous_domain(expression,real,interval)
    defined=_subset(interval,continuous)
    if defined is not True:return defined,None,None,None
    if singleton:return True,None,None,None
    interior=Interval.open(interval.start,interval.end) if isinstance(interval,Interval) else interval
    derivative=expression.diff(real)
    differentiable=continuous_domain(derivative,real,interior)
    if _subset(interior,differentiable) is not True:
        # A missing derivative alone cannot disprove monotonicity.
        return None,None,None,None
    return None,derivative,real,interior


def _direction(derivative, symbol, interior, increasing, strict):
    from .. import Interval
    if derivative==0:return not strict
    condition=derivative>=0 if increasing else derivative<=0
    if condition==S.true:
        nonnegative=True
    elif condition==S.false:
        nonnegative=False
    elif not condition.has(symbol):
        return None
    else:
        allowed=inequality(condition,symbol,interior)
        nonnegative=_subset(interior,allowed)
    if nonnegative is not True or not strict:return nonnegative
    stationary=zeros(derivative,symbol,interior)
    if stationary==S.EmptySet or stationary.is_FiniteSet:return True
    # A continuous nonnegative derivative with finitely many zeros gives a
    # strict change on every nontrivial subinterval. A zero interval gives
    # a plateau; other zero-set representations remain undecided.
    components=stationary.args if stationary.is_Union else (stationary,)
    if any(isinstance(component,Interval) and component.is_empty is False for component in components):
        return False
    return None


def _check(expression, interval, symbol, increasing, strict):
    result,derivative,real,interior=_context(expression,interval,symbol)
    if derivative is None:return result
    return _direction(derivative,real,interior,increasing,strict)


def is_increasing(expression, interval=S.Reals, symbol=None):
    return _check(expression,interval,symbol,True,False)


def is_strictly_increasing(expression, interval=S.Reals, symbol=None):
    return _check(expression,interval,symbol,True,True)


def is_decreasing(expression, interval=S.Reals, symbol=None):
    return _check(expression,interval,symbol,False,False)


def is_strictly_decreasing(expression, interval=S.Reals, symbol=None):
    return _check(expression,interval,symbol,False,True)


def is_monotonic(expression, interval=S.Reals, symbol=None):
    result,derivative,real,interior=_context(expression,interval,symbol)
    if derivative is None:return result
    increasing=_direction(derivative,real,interior,True,False)
    if increasing is True:return True
    decreasing=_direction(derivative,real,interior,False,False)
    if decreasing is True:return True
    return False if increasing is False and decreasing is False else None


__all__=['is_increasing','is_strictly_increasing','is_decreasing',
    'is_strictly_decreasing','is_monotonic']
