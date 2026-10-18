"""Sequence limits with exact continuous limits and discrete differences."""
import operator

from .. import difference_delta


def _continuous(expression,n):
    from .. import limit,oo,S,simplify,combsimp
    from ..solvers._differential import nodes
    try:
        expression=simplify(combsimp(expression))
        unsupported={'factorial','factorial2','subfactorial','gamma','binomial','RisingFactorial','FallingFactorial'}
        if any(node.has(n) and (node.is_UndefinedFunction or
               getattr(node.func,'__name__','') in unsupported) for node in nodes(expression)):
            return None
        result=limit(expression,n,oo)
    except (NotImplementedError,ValueError):
        return None
    return None if result.is_Limit or result in (S.NaN,S.ComplexInfinity) else result


def _limit_seq(expression,n,trials):
    from .. import S,oo,cancel,simplify,combsimp,Dummy,Abs
    seen=set()
    for _ in range(trials):
        expression=simplify(combsimp(expression))
        if expression in seen:return None
        seen.add(expression)
        result=_continuous(expression,n)
        if result is not None:return result
        step_ratio=combsimp(expression.subs(n,n+1)/expression)
        ratio=_continuous(Abs(step_ratio),n)
        if ratio is not None and ratio.is_nonnegative is True and (1-ratio).is_positive is True:
            # The ratio test proves decay even when the continuous engine
            # cannot evaluate factorials or other hypergeometric sequences.
            return S.Zero
        if ratio is not None and (ratio==oo or (ratio-1).is_positive is True) and step_ratio.is_positive is True:
            initial=expression.subs(n,1).doit()
            if initial.is_finite is True:
                if initial.is_positive is True:return oo
                if initial.is_negative is True:return -oo
        numerator,denominator=expression.as_numer_denom()
        if not numerator.has(n) or not denominator.has(n):
            evaluated=expression.doit()
            return _continuous(evaluated,n) if evaluated!=expression else None
        # Stolz--Cesaro is used only with a certified divergent denominator.
        # Positive increasing finite sums provide the common discrete case.
        delta_den=combsimp(difference_delta(denominator,n))
        den_limit=_continuous(denominator,n)
        diverges=den_limit in (oo,-oo)
        positive_sum=False
        if not diverges and denominator.is_Sum:
            limits=denominator.limits
            if len(limits)==1 and len(limits[0])==3 and limits[0][2]==n:
                # A positive summand bounded away from zero forces divergence.
                summand=denominator.function
                k=limits[0][0]
                term_limit=_continuous(summand,k)
                diverges=(term_limit==oo or term_limit is not None and term_limit.is_positive is True)
                lower=limits[0][1]
                if lower.is_Integer:
                    offset=Dummy('sequence_sum_offset',integer=True,nonnegative=True)
                    ratio=combsimp(summand.subs(k,k+1)/summand)
                    positive_sum=(summand.subs(k,lower).is_positive is True
                                  and ratio.subs(k,lower+offset).is_positive is True)
                    growth=_continuous(ratio,k)
                    if positive_sum and growth is not None and (growth-1).is_positive is True:
                        diverges=True
        increment_limit=_continuous(delta_den,n)
        monotone=(positive_sum or delta_den.is_positive is True or delta_den.is_negative is True
                  or increment_limit in (oo,-oo)
                  or increment_limit is not None and (increment_limit.is_positive is True or increment_limit.is_negative is True))
        if not diverges or not monotone or delta_den==0:
            return None
        expression=cancel(combsimp(difference_delta(numerator,n)/delta_den))
    return _continuous(expression,n)


def limit_seq(expr,n=None,trials=5):
    from .. import sympify,Dummy,S,simplify,Abs,Min,Max,sin,cos,Sum
    from ..calculus.accumulationbounds import AccumBounds
    from ..solvers._differential import nodes
    expression=sympify(expr)
    if n is None:
        free=expression.free_symbols
        if not free:return expression
        if len(free)!=1:
            raise ValueError('Specify the index for an expression with several free symbols')
        n=next(iter(free))
    n=sympify(n)
    if not n.is_Symbol:raise ValueError('A sequence index must be a symbol')
    if n not in expression.free_symbols:return expression
    trials=operator.index(trials)
    if not 0<trials<=64:raise ValueError('Sequence limit trials must lie between 1 and 64')
    positive=Dummy('sequence_limit_index',integer=True,positive=True)
    expression=expression.xreplace({n:positive})
    # Exact Binet formulas preserve integral indices and their alternating
    # component, which is then handled by the even/odd subsequence analysis.
    from .. import sqrt
    phi=(1+sqrt(5))/2
    binet={}
    for node in nodes(expression):
        head=getattr(node.func,'__name__','')
        if head in ('fibonacci','lucas') and len(node.args)==1 and node.args[0].is_integer is True:
            index=node.args[0]
            second=(-1)**index*phi**(-index)
            binet[node]=(phi**index-second)/sqrt(5) if head=='fibonacci' else phi**index+second
    expression=expression.xreplace(binet)
    alternating=any(node.is_Pow and node.base.is_negative is True and node.exp.has(positive)
                    for node in nodes(expression)) or expression.has(sin,cos)
    if alternating:
        even=_limit_seq(simplify(expression.xreplace({positive:2*positive})),positive,trials)
        odd=_limit_seq(simplify(expression.xreplace({positive:2*positive+1})),positive,trials)
        if even is not None and odd is not None:
            from .. import cancel
            if even==odd or simplify(cancel(even-odd))==0:return simplify(cancel(even))
            if even.is_extended_real is True and odd.is_extended_real is True:
                return AccumBounds(Min(even,odd),Max(even,odd))
            return None
    result=_limit_seq(expression,positive,trials)
    if result is not None:return result
    if expression.is_Add:
        terms=[limit_seq(term,positive,trials) for term in expression.args]
        if all(term is not None for term in terms):
            value=sum(terms,S.Zero)
            return None if value in (S.NaN,S.ComplexInfinity) else value
    if not expression.has(Sum):
        absolute=_limit_seq(Abs(expression),positive,trials)
        if absolute==0:return S.Zero
    return None


__all__=['difference_delta','limit_seq']
