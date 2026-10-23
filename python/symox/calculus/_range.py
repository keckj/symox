"""Images of elementary continuous real interval components."""
from .. import S
from ._analysis import scalar, variable, real_domain, zeros, real_constants
from ._extrema import _components, minimum, maximum, MAX_CANDIDATES


def _complete_period(domain, period):
    from .. import Interval
    for part in _components(domain):
        if part==S.Reals:return True
        if isinstance(part,Interval):
            if part.start==-S.Infinity or part.end==S.Infinity:return True
            width=part.end-part.start
            if (width>period)==S.true:return True
            if width==period and not (part.left_open and part.right_open):return True
    return False


def _attained(expression, symbol, value, domain):
    if value in (S.Infinity,-S.Infinity):return False
    solutions=zeros(expression-value,symbol,domain)
    if solutions.is_empty is not None:return not solutions.is_empty
    if solutions.is_Union:
        states=[part.is_empty for part in solutions.args]
        if False in states:return True
        if all(state is True for state in states):return False
    raise NotImplementedError('Range endpoint attainment is undecided')


def _real_line_as_interval(value):
    # The reference builds the whole line as Interval(-oo, oo), not Reals.
    from .. import Interval
    if value==S.Reals:return Interval(-S.Infinity,S.Infinity)
    return value


def function_range(f, symbol, domain):
    from .. import Dummy, Interval, FiniteSet, Union
    from .util import continuous_domain, periodicity
    f=scalar(f)
    symbol=variable(symbol)
    domain=real_domain(domain)
    if symbol in domain.free_symbols:raise ValueError('Range bounds must be independent of the variable')
    if domain==S.EmptySet:return domain
    if any(getattr(parameter,'is_real',None) is not True for parameter in f.free_symbols-{symbol}):
        raise NotImplementedError('Range parameters must be known to be real')
    real=Dummy('range_x',real=True)
    f=f.xreplace({symbol:real})
    real_constants(f,real)
    period=periodicity(f,real)
    if period is not None and period.is_positive is True and period.is_finite is True and _complete_period(domain,period):
        domain=Interval(0,period)
    continuous=continuous_domain(f,real,domain)
    if continuous==S.EmptySet:return continuous
    if f==real:return _real_line_as_interval(continuous)
    if not f.has(real) and continuous.is_empty is False:return FiniteSet(f)
    images=[]
    points=set()
    for part in _components(continuous):
        if part.is_FiniteSet:
            points.update(f.subs(real,point) for point in part.args)
            if len(points)>MAX_CANDIDATES:raise NotImplementedError('Range exceeds 2048 finite values')
            continue
        low=minimum(f,real,part)
        high=maximum(f,real,part)
        if low==high:
            if low in (S.Infinity,-S.Infinity):
                raise NotImplementedError('An infinite constant range is unresolved')
            images.append(FiniteSet(low))
            continue
        compact=(isinstance(part,Interval) and not part.left_open and not part.right_open
                 and part.start.is_finite is True and part.end.is_finite is True)
        lower_attained=compact or _attained(f,real,low,part)
        upper_attained=compact or _attained(f,real,high,part)
        images.append(Interval(low,high,not lower_attained,not upper_attained))
    result=Union(*images,FiniteSet(*points))
    return _real_line_as_interval(result)


__all__=['function_range']
