"""Real continuity, stationary points and structural periodicity."""
from ._analysis import scalar, variable, walk, name, zeros, real_domain, inequality
from .. import S as _S


def lcim(numbers):
    """Least common integral multiple of commensurable real constants."""
    import math
    from .. import sympify,cancel,simplify,Abs,Rational,S
    values=[sympify(value) for value in numbers]
    if not values:return S.One
    if any(value.free_symbols or value.is_finite is not True or value.is_real is not True
           for value in values):return None
    if any(value==0 for value in values):
        return S.Zero if all(value.is_Rational for value in values) else None
    unit=values[0]
    ratios=[simplify(cancel(value/unit)) for value in values]
    if any(not value.is_Rational for value in ratios):return None
    numerator=math.lcm(*(abs(int(value.p)) for value in ratios))
    denominator=math.gcd(*(int(value.q) for value in ratios))
    return simplify(Abs(unit)*Rational(numerator,denominator))


def is_convex(f, *syms, domain=_S.Reals):
    """Prove univariate convexity through a nondecreasing first derivative."""
    if not syms:raise ValueError('Supply the variable for a convexity check')
    if len(syms)!=1:raise NotImplementedError('Multivariate convexity is not implemented')
    from ._monotonicity import _context, is_increasing
    result,derivative,real,interior=_context(f,domain,syms[0])
    if derivative is None:return result
    return is_increasing(derivative,interior,real)


def _continuous_domain(f, symbol, domain, *, continuous=True):
    from .. import S, Complement, cos, sin, sinh
    f=scalar(f)
    symbol=variable(symbol)
    domain=real_domain(domain)
    result=domain
    entire={'sin','cos','sinh','cosh','tanh','sech','exp','Abs','atan','erf','erfc','asinh'}
    for node in walk(f):
        if not node.has(symbol):
            if node.is_real is False:
                raise NotImplementedError('Expression is not known to be real on this domain')
            continue
        if node.is_Symbol or node.is_Add or node.is_Mul:
            continue
        if node.is_Pow:
            exponent=node.exp
            if exponent.has(symbol):
                if node.base.is_positive is not True:
                    raise NotImplementedError('Variable powers require a positive real base')
            elif exponent.is_integer is True:
                if exponent.is_negative is True:
                    result=Complement(result,zeros(node.base,symbol,domain))
                elif exponent.is_nonnegative is not True:
                    raise NotImplementedError('Integer power sign is undecidable')
            elif exponent.is_Rational:
                # Native powers follow the principal complex branch, so
                # a negative real base is outside the real domain even
                # when the rational denominator is odd.
                condition=node.base>0 if exponent.is_negative is True else node.base>=0
                result=result.intersect(inequality(condition,symbol,domain))
            elif exponent.is_real is True:
                if exponent.is_positive is True:condition=node.base>=0
                elif exponent.is_negative is True:condition=node.base>0
                else:raise NotImplementedError('Real power domain requires a known exponent sign')
                result=result.intersect(inequality(condition,symbol,domain))
            else:
                raise NotImplementedError('Power domain requires a real, known exponent')
            continue
        head=name(node)
        if not continuous and head in ('sign','floor','ceiling','Heaviside','acot'):
            # Their real arguments can have poles or restricted real domains;
            # the walk still checks those subexpressions below.
            continue
        if head=='Abs' and node.args[0].is_real is not True:
            raise NotImplementedError('Absolute-value domains with complex intermediate values are unsupported')
        if head in entire:
            continue
        if head=='log':condition=node.args[0]>0
        elif head in ('asin','acos'):
            result=result.intersect(inequality(node.args[0]>=-1,symbol,domain))
            condition=node.args[0]<=1
        elif head=='acosh':condition=node.args[0]>=1
        elif head=='atanh':
            result=result.intersect(inequality(node.args[0]>-1,symbol,domain))
            condition=node.args[0]<1
        elif head=='acoth':
            from .. import Union
            good=Union(inequality(node.args[0]<-1,symbol,domain),
                       inequality(node.args[0]>1,symbol,domain))
            result=result.intersect(good)
            continue
        elif head in ('tan','sec','cot','csc','coth','csch','acot'):
            kernel=cos if head in ('tan','sec') else sin if head in ('cot','csc') else sinh
            base=node.args[0] if head=='acot' else kernel(node.args[0])
            result=Complement(result,zeros(base,symbol,domain))
            continue
        else:
            raise NotImplementedError('Continuous domain is unsupported for %s' % (head or node.func))
        result=result.intersect(inequality(condition,symbol,domain))
    return result


def continuous_domain(f, symbol, domain):
    """Real domains of elementary scalar expressions on native sets."""
    from .. import Dummy, S
    f=scalar(f)
    symbol=variable(symbol)
    domain=real_domain(domain)
    if domain==S.EmptySet:return domain
    if symbol in domain.free_symbols:
        raise ValueError('Analysis domain bounds must be independent of the variable')
    real_symbol=Dummy('domain_x',real=True)
    expression=f.xreplace({symbol:real_symbol})
    return _continuous_domain(expression,real_symbol,domain).xreplace({real_symbol:symbol})


def stationary_points(f, symbol, domain=None):
    from .. import S, Dummy
    f=scalar(f)
    symbol=variable(symbol)
    domain=S.Reals if domain is None else real_domain(domain)
    if domain==S.EmptySet:return domain
    if symbol in domain.free_symbols:
        raise ValueError('Analysis domain bounds must be independent of the variable')
    real_symbol=Dummy('stationary_x',real=True)
    expression=f.xreplace({symbol:real_symbol})
    continuous=_continuous_domain(expression,real_symbol,domain)
    derivative=expression.diff(real_symbol)
    if any(name(node) in ('Abs','sign') and node.has(real_symbol) for node in walk(derivative)):
        raise NotImplementedError('Stationary points across nondifferentiable branches are unsupported')
    return zeros(derivative,real_symbol,continuous).xreplace({real_symbol:symbol})


def _common_period(periods):
    from .. import S, cancel
    nonzero=[period for period in periods if period!=0]
    if not nonzero:return S.Zero
    if any(period is None for period in nonzero):return None
    period=nonzero[0]
    for candidate in nonzero[1:]:
        ratio=cancel(candidate/period)
        if ratio.is_Rational is not True:
            return None
        numerator,_=ratio.as_numer_denom()
        period *= abs(int(numerator))
    return period


def periodicity(f, symbol, check=False):
    """Find a real period; return None when this structural search finds none."""
    from .. import S, Dummy, pi, Abs, im, simplify, trigsimp
    from .._core import Expr
    f=scalar(f)
    symbol=variable(symbol)
    # Periodicity is with respect to a real independent variable. This
    # also exposes simplifications such as Abs(exp(I*x)) == 1.
    real_symbol=Dummy('period_x',real=True)
    expression=f.xreplace({symbol:real_symbol})
    expression=trigsimp(expression)
    list(walk(expression))
    cache={}
    def affine(argument):
        slope=argument.diff(real_symbol)
        return slope if not slope.has(real_symbol) else None
    def candidate(node):
        if not isinstance(node,Expr):
            return S.Zero if hasattr(node,'has') and not node.has(real_symbol) else None
        if node in cache:return cache[node]
        if not node.has(real_symbol):return S.Zero
        head=name(node)
        period=None
        if node.is_Add or node.is_Mul:
            period=_common_period([candidate(arg) for arg in node.args])
        elif node.is_Pow:
            if not node.exp.has(real_symbol):
                period=candidate(node.base)
                if period and node.exp.is_even is True:
                    shifted=node.base.subs(real_symbol,real_symbol+period/2)
                    if trigsimp(shifted+node.base)==0:period/=2
            elif not node.base.has(real_symbol):
                period=candidate(node.exp)
        elif head in ('sin','cos','tan','cot','sec','csc','sinh','cosh','tanh','coth','sech','csch'):
            slope=affine(node.args[0])
            if slope is not None and slope.is_zero is False:
                multiplier=1 if head in ('tan','cot','tanh','coth') else 2
                if head in ('sin','cos','tan','cot','sec','csc') and slope.is_real is True:
                    period=multiplier*pi/Abs(slope)
                elif head in ('sinh','cosh','tanh','coth','sech','csch') and slope.is_imaginary is True:
                    period=multiplier*pi/Abs(im(slope))
            if period is None:period=candidate(node.args[0])
        elif head=='exp':
            slope=affine(node.args[0])
            if slope is not None and slope.is_imaginary is True:
                period=2*pi/Abs(im(slope))
            else:period=candidate(node.args[0])
        elif head=='Mod':
            numerator,modulus=node.args
            slope=affine(numerator)
            if (slope is not None and slope.is_real is True and slope.is_zero is False
                    and not modulus.has(real_symbol) and modulus.is_real is True
                    and modulus.is_zero is False):
                period=Abs(modulus/slope)
        elif node.is_Derivative:
            period=candidate(node.expr)
        elif node.is_Function or node.is_UndefinedFunction:
            period=_common_period([candidate(arg) for arg in node.args])
            if head=='Abs' and period:
                original=node.args[0]
                if trigsimp(original.subs(real_symbol,real_symbol+period/2)+original)==0:period/=2
        cache[node]=period
        return period
    period=candidate(expression)
    if check and period is not None and period!=0:
        difference=expression.subs(real_symbol,real_symbol+period)-expression
        if simplify(trigsimp(difference))!=0:
            raise NotImplementedError('The candidate period could not be verified symbolically')
    return period


from ._extrema import minimum, maximum
from ._range import function_range
from ._notempty import not_empty_in

__all__=['continuous_domain','stationary_points','periodicity','is_convex','minimum','maximum','function_range','not_empty_in','lcim']
