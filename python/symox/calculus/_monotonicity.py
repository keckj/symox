"""Monotonicity proofs on connected real intervals using native sets."""
from .. import S
from ._analysis import scalar, variable, real_domain, inequality, zeros, real_constants


def _subset(domain, allowed):
    missing=domain-allowed
    if missing.is_empty is not None:return missing.is_empty
    return domain.is_subset(allowed)


class _Reference(Exception):
    """The native proof declines this input; use the reference algorithm."""


def _reference_variable(expression, symbol):
    from .. import Symbol, sympify
    expression=sympify(expression)
    free=expression.free_symbols
    if symbol is None and len(free)>1:
        raise NotImplementedError('The function has not yet been implemented'
            ' for all multivariate expressions.')
    return expression,symbol or (free.pop() if free else Symbol('x'))


def _reference_helper(expression, predicate, interval, symbol):
    # Reference monotonicity_helper: compare the interval with the solution
    # set of the derivative predicate over the reals.
    from .. import solveset
    expression,var=_reference_variable(expression,symbol)
    derivative=expression.diff(var)
    predicate_interval=solveset(predicate(derivative),var,S.Reals)
    return interval.is_subset(predicate_interval)


def _reference_monotonic(expression, interval, symbol):
    from .. import solveset
    expression,var=_reference_variable(expression,symbol)
    turning_points=solveset(expression.diff(var),var,interval)
    return interval.intersection(turning_points) is S.EmptySet


def _context(expression, interval, symbol):
    from .. import Dummy, Interval
    from .util import continuous_domain
    expression=scalar(expression)
    if symbol is None:
        free=expression.free_symbols
        if len(free)>1:raise NotImplementedError('Supply the variable for a multivariate monotonicity check')
        symbol=next(iter(free)) if free else Dummy('constant_x',real=True)
    symbol=variable(symbol)
    try:
        interval=real_domain(interval)
    except NotImplementedError:
        raise _Reference from None
    if symbol in interval.free_symbols:
        # The reference accepts bounds depending on the variable.
        raise _Reference
    if interval==S.EmptySet:return True,None,None,None
    singleton=interval.is_FiniteSet and len(interval.args)==1
    if interval!=S.Reals and not isinstance(interval,Interval) and not singleton:
        # Disconnected or non-interval domains: the reference decides them
        # through the solution set of the derivative predicate.
        raise _Reference
    parameters=expression.free_symbols-{symbol}
    if any(getattr(parameter,'is_real',None) is not True for parameter in parameters):
        raise _Reference
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


_PREDICATES={
    (True,False):lambda d:d>=0,
    (True,True):lambda d:d>0,
    (False,False):lambda d:d<=0,
    (False,True):lambda d:d<0,
}


def _check(expression, interval, symbol, increasing, strict):
    try:
        result,derivative,real,interior=_context(expression,interval,symbol)
        if derivative is not None:
            result=_direction(derivative,real,interior,increasing,strict)
    except _Reference:
        result=None
    if result is None:
        # Undecided natively: give the reference answer.
        result=_reference_helper(expression,_PREDICATES[increasing,strict],interval,symbol)
    return result


def is_increasing(expression, interval=S.Reals, symbol=None):
    return _check(expression,interval,symbol,True,False)


def is_strictly_increasing(expression, interval=S.Reals, symbol=None):
    return _check(expression,interval,symbol,True,True)


def is_decreasing(expression, interval=S.Reals, symbol=None):
    return _check(expression,interval,symbol,False,False)


def is_strictly_decreasing(expression, interval=S.Reals, symbol=None):
    return _check(expression,interval,symbol,False,True)


def _native_monotonic(expression, interval, symbol):
    result,derivative,real,interior=_context(expression,interval,symbol)
    if derivative is None:return result
    increasing=_direction(derivative,real,interior,True,False)
    if increasing is True:return True
    decreasing=_direction(derivative,real,interior,False,False)
    if decreasing is True:return True
    return False if increasing is False and decreasing is False else None


def is_monotonic(expression, interval=S.Reals, symbol=None):
    try:
        result=_native_monotonic(expression,interval,symbol)
    except _Reference:
        result=None
    if result is None:
        result=_reference_monotonic(expression,interval,symbol)
    return result


__all__=['is_increasing','is_strictly_increasing','is_decreasing',
    'is_strictly_decreasing','is_monotonic']
