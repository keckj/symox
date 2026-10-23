"""Real continuity, stationary points and structural periodicity."""
from ._analysis import scalar, variable, walk, name, zeros, real_domain, inequality
from .. import S as _S
S = _S


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


def _trig_period(f, symbol):
    """SymPy's ``TrigonometricFunction._period``; None when it raises."""
    from .. import S, pi, expand_mul
    general = pi if name(f) in ('tan', 'cot') else 2*pi
    arg = expand_mul(f.args[0])
    if not arg.has(symbol):
        return S.Zero
    if arg == symbol:
        return general
    if symbol in arg.free_symbols:
        if arg.is_Mul:
            g, h = arg.as_independent(symbol)
            if h == symbol:
                return general/abs(g)
        if arg.is_Add:
            a, h = arg.as_independent(symbol)
            g, h = h.as_independent(symbol, as_Add=False)
            if h == symbol:
                return general/abs(g)
    return None


def _lcim(numbers):
    """SymPy's ``lcim``: least common integral multiple, or None."""
    from .. import lcm_list
    if all(num.is_irrational for num in numbers):
        factors_num = [num.factor().as_coeff_Mul() for num in numbers]
        term = factors_num[0][1]
        if all(factor == term for coeff, factor in factors_num):
            return lcm_list([coeff for coeff, factor in factors_num])*term
        return None
    if all(num.is_rational for num in numbers):
        return lcm_list(numbers)
    return None


_TRIG = ('sin', 'cos', 'tan', 'cot', 'sec', 'csc')


def periodicity(f, symbol, check=False):
    """SymPy's ``periodicity``: the period of ``f`` in ``symbol``, ``S.Zero``
    for a constant and None when no period is found."""
    from .. import (S, Dummy, Abs, Pow, Mod, Piecewise, Relational, sin, re,
                    im, expand_mul, degree, sympify)
    from ..core.numbers import equal_valued
    from .. import cos, tan, cot, sec, csc
    trig_heads = (sin, cos, tan, cot, sec, csc)
    f = sympify(f)
    symbol = variable(symbol)
    if getattr(symbol, 'kind', None) is not None and str(symbol.kind) != 'NumberKind':
        raise NotImplementedError("Cannot use symbol of kind %s" % symbol.kind)
    temp = Dummy('x', real=True)
    f = f.subs(symbol, temp)
    symbol = temp

    def _check(orig_f, period):
        new_f = orig_f.subs(symbol, symbol + period)
        if new_f.equals(orig_f):
            return period
        raise NotImplementedError(
            "The period of the given function cannot be verified. When `%s` "
            "was replaced with `%s + %s` in `%s`, the result was `%s` which "
            "was not recognized as being the same as the original function. "
            "So either the period was wrong or the two forms were not "
            "recognized as being equal. Set check=False to obtain the value."
            % (symbol, symbol, period, orig_f, new_f))

    orig_f = f
    period = None
    if isinstance(f, Relational):
        f = f.lhs - f.rhs
    simplified = f.simplify()
    # SymPy's simplify replaces a trigonometric expression by its deep
    # trigsimp form (kept unless 1.7 times longer than the input); the
    # native simplify does not always apply that step (sin(x)*cos(x)).
    if simplified.has(*trig_heads):
        from .. import trigsimp, count_ops
        candidate = trigsimp(simplified, deep=True)
        if count_ops(candidate) <= 1.7*count_ops(f):
            simplified = candidate
    f = simplified
    if symbol not in f.free_symbols:
        return S.Zero
    if name(f) in _TRIG:
        period = _trig_period(f, symbol)
    if isinstance(f, Abs):
        arg = f.args[0]
        if name(arg) in ('sec', 'csc', 'cos'):
            arg = sin(arg.args[0])
        period = periodicity(arg, symbol)
        if period is not None and name(arg) == 'sin':
            orig_f = Abs(arg)
            try:
                return _check(orig_f, period/2)
            except NotImplementedError as err:
                if check:
                    raise NotImplementedError(err)
    if name(f) == 'exp' or (f.is_Pow and f.base == S.Exp1):
        f = Pow(S.Exp1, expand_mul(f.exp))
        if im(f) != 0:
            period_real = periodicity(re(f), symbol)
            period_imag = periodicity(im(f), symbol)
            if period_real is not None and period_imag is not None:
                period = _lcim([period_real, period_imag])
    if f.is_Pow and f.base != S.Exp1:
        base, expo = f.args
        base_has_sym = base.has(symbol)
        expo_has_sym = expo.has(symbol)
        if base_has_sym and not expo_has_sym:
            period = periodicity(base, symbol)
        elif expo_has_sym and not base_has_sym:
            period = periodicity(expo, symbol)
        else:
            period = _periodicity(f.args, symbol)
    elif f.is_Mul:
        coeff, g = f.as_independent(symbol, as_Add=False)
        if name(g) in _TRIG or not equal_valued(coeff, 1):
            period = periodicity(g, symbol)
        else:
            period = _periodicity(g.args, symbol)
    elif f.is_Add:
        k, g = f.as_independent(symbol)
        if k is not S.Zero:
            return periodicity(g, symbol)
        period = _periodicity(g.args, symbol)
    elif isinstance(f, Mod):
        a, n = f.args
        if a == symbol:
            period = n
        elif name(a) in _TRIG:
            period = periodicity(a, symbol)
        elif (a.is_polynomial(symbol) and degree(a, symbol) == 1 and
                symbol not in n.free_symbols):
            period = Abs(n/a.diff(symbol))
    elif isinstance(f, Piecewise):
        pass
    elif period is None:
        from ..solvers.decompogen import compogen, decompogen
        g_s = decompogen(f, symbol)
        num_of_gs = len(g_s)
        if num_of_gs > 1:
            for index, g in enumerate(reversed(g_s)):
                start_index = num_of_gs - 1 - index
                g = compogen(g_s[start_index:], symbol)
                if g not in (orig_f, f):
                    period = periodicity(g, symbol)
                    if period is not None:
                        break
    if period is not None:
        if check:
            return _check(orig_f, period)
        return period
    return None


def _periodicity(args, symbol):
    periods = []
    for f in args:
        period = periodicity(f, symbol)
        if period is None:
            return None
        if period is not S.Zero and period != 0:
            periods.append(period)
    if len(periods) > 1:
        return _lcim(periods)
    if periods:
        return periods[0]
    return None


from ._extrema import minimum, maximum
from ._range import function_range
from ._notempty import not_empty_in

__all__=['continuous_domain','stationary_points','periodicity','is_convex','minimum','maximum','function_range','not_empty_in','lcim']
