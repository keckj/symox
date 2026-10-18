"""Real parameter domains where a finite symbolic set has a member."""
from ._analysis import scalar, variable, real_constants, inequality, zeros


def _constant_witness(specification, symbol, depth=0):
    from .. import S, Intersection
    if depth>256:raise NotImplementedError('Nonempty-domain nesting exceeds 256 levels')
    if specification.is_Union:
        return any(_constant_witness(branch,symbol,depth+1) for branch in specification.args)
    factors=specification.args if specification.is_Intersection else (specification,)
    for finite in factors:
        if not finite.is_FiniteSet:continue
        remaining=list(factors)
        remaining.remove(finite)
        target=Intersection(S.Reals,*remaining)
        for candidate in finite.args:
            if not candidate.has(symbol) and candidate.is_real is True:
                if target.contains(candidate)==S.true:return True
    return False


def _definition(expression, symbol):
    from .. import S
    from .util import _continuous_domain
    expression=scalar(expression)
    real_constants(expression,symbol)
    return _continuous_domain(expression,symbol,S.Reals,continuous=False)


def _condition(condition, symbol):
    from .. import S, true, false
    if condition==true:return S.Reals
    if condition==false:return S.EmptySet
    return inequality(condition,symbol,S.Reals)


def _preimage(expression, target, symbol, depth=0):
    from .. import S, Union, Intersection, Complement, Interval
    if depth>256:raise NotImplementedError('Nonempty-domain nesting exceeds 256 levels')
    if target==S.EmptySet:return S.EmptySet
    if target==S.Reals:return _definition(expression,symbol)
    if target.is_Union:
        return Union(*[_preimage(expression,part,symbol,depth+1) for part in target.args])
    if target.is_Intersection:
        return Intersection(*[_preimage(expression,part,symbol,depth+1) for part in target.args])
    if target.is_Complement:
        return Complement(_preimage(expression,target.args[0],symbol,depth+1),
                          _preimage(expression,target.args[1],symbol,depth+1))
    defined=_definition(expression,symbol)
    if target.is_FiniteSet:
        return Union(*[zeros(expression-value,symbol,
                            defined.intersect(_definition(value,symbol))) for value in target.args])
    if isinstance(target,Interval):
        result=defined
        if target.start!=S.NegativeInfinity:
            result=result.intersect(_definition(target.start,symbol))
            condition=expression>target.start if target.left_open else expression>=target.start
            result=result.intersect(_condition(condition,symbol))
        if target.end!=S.Infinity:
            result=result.intersect(_definition(target.end,symbol))
            condition=expression<target.end if target.right_open else expression<=target.end
            result=result.intersect(_condition(condition,symbol))
        return result
    raise NotImplementedError('Unsupported target set for finite-set nonemptiness')


def _domain(specification, symbol, depth=0):
    from .. import S, Set, Union, Intersection
    if depth>256:raise NotImplementedError('Nonempty-domain nesting exceeds 256 levels')
    if not isinstance(specification,Set):raise ValueError('Supply a finite set or its intersection with real sets')
    if specification==S.EmptySet:return S.EmptySet
    if specification.is_Union:
        if _constant_witness(specification,symbol):return S.Reals
        result=S.EmptySet
        for branch in specification.args:
            result=Union(result,_domain(branch,symbol,depth+1))
            if result==S.Reals:return result
        return result
    factors=specification.args if specification.is_Intersection else (specification,)
    finite=[factor for factor in factors if factor.is_FiniteSet]
    if not finite:raise ValueError('A FiniteSet must be supplied')
    candidates=min(finite,key=lambda factor:len(factor.args))
    remaining=list(factors)
    remaining.remove(candidates)
    target=Intersection(S.Reals,*remaining)
    # A constant witness proves nonemptiness even if another candidate cannot
    # be analysed. Inspect these before asking the solver about varying terms.
    if _constant_witness(specification,symbol):return S.Reals
    result=S.EmptySet
    for candidate in candidates.args:
        result=Union(result,_preimage(candidate,target,symbol))
        if result==S.Reals:return result
    return result


def not_empty_in(finset_intersection, *syms):
    """Return real values of one variable for which a finite set is nonempty.

    Members must be real where defined. Native equation and inequality
    solvers handle interval, finite, union and complement membership targets.
    Unsupported functions and undecidable real coefficients raise explicitly.
    """
    from .. import Dummy, Set, sympify, _wrap_set
    if not syms:raise ValueError('Supply one variable for a nonempty-domain calculation')
    if len(syms)!=1:raise NotImplementedError('Multivariate nonempty domains are unsupported')
    symbol=variable(syms[0])
    specification=finset_intersection
    if not isinstance(specification,Set):specification=_wrap_set(sympify(specification))
    if not isinstance(specification,Set):raise ValueError('Supply a finite set or its intersection with real sets')
    real=Dummy('nonempty_x',real=True)
    return _domain(specification.xreplace({symbol:real}),real).xreplace({real:symbol})


__all__=['not_empty_in']
