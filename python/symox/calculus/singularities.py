"""Isolated singularity candidates, using the native equation solver."""
from ._analysis import scalar, variable, walk, name, zeros, inequality


def singularities(expression, symbol, domain=None):
    from .. import S, Set, cos, sin, sinh, cosh, Union
    expression=scalar(expression)
    symbol=variable(symbol)
    domain=(S.Reals if symbol.is_real is True else S.Complexes) if domain is None else domain
    if not isinstance(domain,Set):
        raise TypeError('Singularity domain must be a native symbolic set')
    if domain==S.EmptySet:return domain
    candidates=[]
    reciprocal={'tan':cos,'sec':cos,'cot':sin,'csc':sin,
                'tanh':cosh,'sech':cosh,'coth':sinh,'csch':sinh}
    regular={'exp','sin','cos','sinh','cosh','Abs','sign','floor','ceiling',
             're','im','conjugate','arg','asin','acos','atan','acot',
             'asinh','acosh','atan2','Max','Min'}
    for node in walk(expression):
        if not node.has(symbol):
            continue
        if node.is_UndefinedFunction or node.is_Derivative:
            raise NotImplementedError('Singularities of undefined functions and held derivatives are unknown')
        head=name(node)
        if node.is_Pow:
            if node.exp.is_infinite is True:
                raise NotImplementedError('Singularities of infinite powers are unsupported')
            if node.exp.is_negative is True:
                candidates.append(zeros(node.base,symbol,domain))
            elif node.exp.is_nonnegative is not True and node.base.is_nonzero is not True:
                raise NotImplementedError('Singularities require a known exponent sign')
        elif head in reciprocal:
            candidates.append(zeros(reciprocal[head](node.args[0]),symbol,domain))
        elif head in ('log','asech','acsch'):
            candidates.append(zeros(node.args[0],symbol,domain))
        elif head in ('atanh','acoth'):
            candidates.extend(zeros(node.args[0]+offset,symbol,domain) for offset in (-1,1))
        elif node.is_Atom or node.is_Add or node.is_Mul or head in regular:
            continue
        else:
            raise NotImplementedError('Singularities are unsupported for %s' % (head or node.func))
    return Union(*candidates) if candidates else S.EmptySet


from ._monotonicity import (is_increasing, is_strictly_increasing,
    is_decreasing, is_strictly_decreasing, is_monotonic)


def monotonicity_helper(expression,predicate,interval=None,symbol=None):
    """Prove that a derivative satisfies a predicate throughout an interval."""
    from .. import S,Dummy,sympify
    from ._monotonicity import _context,_subset
    interval=S.Reals if interval is None else interval
    result,derivative,real,interior=_context(expression,interval,symbol)
    if derivative is None:
        if result is not True or interval==S.EmptySet:return result
        expression=scalar(expression)
        if symbol is None:
            symbol=next(iter(expression.free_symbols),Dummy('constant_x',real=True))
        derivative=expression.diff(symbol).subs(symbol,interval.args[0])
        condition=sympify(predicate(derivative))
        return True if condition==S.true else False if condition==S.false else None
    condition=sympify(predicate(derivative))
    if condition==S.true:return True
    if condition==S.false:return False
    if not condition.has(real):return None
    def satisfying(value):
        from .. import Union
        if value==S.true:return interval
        if value==S.false:return S.EmptySet
        head=getattr(value.func,'__name__','')
        if head=='And':
            result=interval
            for part in value.args:result=result.intersect(satisfying(part))
            return result
        if head=='Or':return Union(*(satisfying(part) for part in value.args))
        if head=='Not':return interval-satisfying(value.args[0])
        return inequality(value,real,interval)
    allowed=satisfying(condition)
    return _subset(interval,allowed)

__all__=['singularities','is_increasing','is_strictly_increasing',
    'is_decreasing','is_strictly_decreasing','is_monotonic','monotonicity_helper']
