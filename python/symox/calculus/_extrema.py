"""Infimum and supremum candidates on elementary real domains."""
from .. import S
from ._analysis import scalar, variable, real_domain, zeros, walk, name, real_constants

MAX_CANDIDATES=2048


def _components(domain):
    from .. import Interval
    pending=[domain]
    components=[]
    while pending:
        part=pending.pop()
        if part==S.EmptySet:continue
        if part.is_Union:pending.extend(part.args)
        elif part==S.Reals or isinstance(part,Interval) or part.is_FiniteSet:
            components.append(part)
        else:raise NotImplementedError('Extrema require explicit real intervals or finite sets')
        if len(components)+len(pending)>MAX_CANDIDATES:
            raise NotImplementedError('Extrema exceed 2048 domain components')
    return components


def _full_period(domain, period):
    from .. import Interval
    for part in _components(domain):
        if part==S.Reals:return True
        if isinstance(part,Interval):
            if part.start==-S.Infinity or part.end==S.Infinity:return True
            if (part.end-part.start>=period)==S.true:return True
    return False


def _finite_points(result):
    if result==S.EmptySet:return ()
    if result.is_FiniteSet:
        if len(result.args)>MAX_CANDIDATES:raise NotImplementedError('Extrema exceed 2048 critical points')
        return result.args
    if result.is_Union:
        return tuple(point for part in result.args for point in _finite_points(part))
    raise NotImplementedError('Extrema require a finite, complete critical-point set')


def _stationary(derivative, symbol, domain, depth=0):
    from .. import Union
    if depth>256:raise NotImplementedError('Critical-point factorization exceeds depth 256')
    if derivative.is_zero is False or derivative.is_nonzero is True:
        return S.EmptySet
    if derivative.is_Mul:
        return Union(*(_stationary(factor,symbol,domain,depth+1) for factor in derivative.args))
    if derivative.is_Pow:
        if derivative.exp.is_negative is True:return S.EmptySet
        if derivative.exp.is_positive is True:
            return _stationary(derivative.base,symbol,domain,depth+1)
    head=name(derivative)
    if head=='exp':return S.EmptySet
    if head=='sign':return zeros(derivative.args[0],symbol,domain)
    return zeros(derivative,symbol,domain)


def _values(f, symbol, domain, upper):
    from .. import Dummy, Interval, trigsimp
    from .util import continuous_domain, periodicity
    from .singularities import singularities
    f=scalar(f)
    symbol=variable(symbol)
    domain=real_domain(domain)
    if domain==S.EmptySet:raise ValueError('Extrema are not defined for an empty domain')
    if symbol in domain.free_symbols:raise ValueError('Extrema bounds must be independent of the variable')
    if any(getattr(parameter,'is_real',None) is not True for parameter in f.free_symbols-{symbol}):
        raise NotImplementedError('Extrema parameters must be known to be real')
    real=Dummy('extrema_x',real=True)
    f=f.xreplace({symbol:real})
    real_constants(f,real)
    period=periodicity(f,real)
    if period is not None and period.is_positive is True and period.is_finite is True and _full_period(domain,period):
        domain=Interval(0,period)
    continuous=continuous_domain(f,real,domain)
    if continuous==S.EmptySet:raise ValueError('The requested domain has no real function values')
    if not f.has(real) and continuous.is_empty is False:return {f}
    derivative=trigsimp(f.diff(real))
    values=set()
    def add(value):
        value=scalar(value)
        if value.has(real) or value.is_Limit or value==S.NaN or value==S.ComplexInfinity:
            raise NotImplementedError('An extremum candidate or endpoint limit is unresolved')
        if value.has(S.Infinity,-S.Infinity) and value not in (S.Infinity,-S.Infinity):
            raise NotImplementedError('Extrema require a resolved sign for infinite endpoint limits')
        values.add(value.xreplace({real:symbol}))
        if len(values)>MAX_CANDIDATES:raise NotImplementedError('Extrema exceed 2048 candidate values')
    intervals=[]
    absorbing=S.Infinity if upper else -S.Infinity
    for part in _components(continuous):
        if part.is_FiniteSet:
            for point in part.args:
                add(f.subs(real,point))
                if absorbing in values:return values
            continue
        if part==S.Reals:
            left,right=-S.Infinity,S.Infinity
            left_open=right_open=True
        else:
            left,right=part.start,part.end
            left_open,right_open=part.left_open,part.right_open
            if (left<right)!=S.true:
                raise NotImplementedError('Extrema require an ordered nonempty interval')
        if derivative==0:
            sample=(left+right)/2 if left.is_finite is True and right.is_finite is True else (
                left+1 if left.is_finite is True else right-1 if right.is_finite is True else S.Zero)
            add(f.subs(real,sample))
            continue
        add(f.limit(real,left,dir='+') if left_open else f.subs(real,left))
        add(f.limit(real,right,dir='-') if right_open else f.subs(real,right))
        if absorbing in values:return values
        intervals.append(part)
    for part in intervals:
        critical=_stationary(derivative,real,part)
        # Formal derivative poles can be removable at an extremum of f.
        # Include such points whenever f itself is continuous there.
        extras=[singularities(derivative,real,part)]
        extras.extend(zeros(node.args[0],real,part) for node in walk(f)
                      if name(node)=='Abs' and node.has(real))
        from .. import Union
        for point in _finite_points(Union(critical,*extras)):
            add(f.subs(real,point))
    return values


def minimum(f, symbol, domain=S.Reals):
    from .. import Min
    return Min(*_values(f,symbol,domain,False))


def maximum(f, symbol, domain=S.Reals):
    from .. import Max
    return Max(*_values(f,symbol,domain,True))


__all__=['minimum','maximum']
