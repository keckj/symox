"""Absorb redundant arbitrary integration constants using native algebra."""


def _children(value):
    from ..._core import Expr
    from ..._matrixexpr_methods import _children as children
    return [item for item in children(value) if isinstance(item,Expr)]


def _constant_parts(expression, constants):
    from ... import Add, Mul
    from ..._core import Expr
    pending, seen, result = [expression], set(), []
    while pending:
        value = pending.pop()
        if value in seen:
            continue
        seen.add(value)
        if len(seen)>10000:
            raise NotImplementedError('Constant simplification exceeds 10000 graph nodes')
        symbols = value.free_symbols
        if symbols and symbols <= constants:
            result.append(value)
            continue
        if getattr(value.func,'__name__','') == 'exp':
            value = value.expand(mul=True)
        if value.is_Add or value.is_Mul:
            independent = [arg for arg in value.args if arg.free_symbols <= constants]
            if len(independent)>1:
                part = (Add if value.is_Add else Mul)(*independent)
                if part.free_symbols:
                    result.append(part)
        elif value.is_Integral:
            if symbols <= constants and all(len(limit)==3 for limit in value.limits):
                result.append(value)
        pending.extend(reversed(_children(value)))
    return result


def _linear_redundancies(expression, constants):
    from ... import Add, Mul, Equality, Eq, S
    from .._differential import nodes
    counts = {constant:expression.count(constant) for constant in constants}
    constants = sorted((constant for constant in constants if counts[constant]),key=str)
    replacements = {}
    for original in nodes(expression):
        value = original.xreplace(replacements)
        if value.is_Add:
            groups = {}
            for constant in constants:
                if value.count(constant)==counts[constant] and value.diff(constant,2)==0:
                    groups.setdefault(value.diff(constant),[]).append(constant)
            discard = {constant:S.Zero for group in groups.values() for constant in group[1:]}
            if discard:
                value = value.xreplace(discard)
        if value != original:
            replacements[original] = value
    expression = expression.xreplace(replacements)
    if not isinstance(expression,Equality):
        return expression
    lhs,rhs = expression.lhs,expression.rhs
    if lhs.is_Symbol and lhs in constants:
        lhs,rhs = rhs,lhs
    additive = lambda value:value.is_Add or value.is_Symbol
    multiplicative = lambda value:value.is_Mul or value.is_Symbol
    constant_term = lambda value:bool(value.is_Number) or value in constants
    if additive(lhs) and additive(rhs):
        left,right = Add.make_args(lhs),Add.make_args(rhs)
        lhs = Add(*(value for value in left if not constant_term(value)))-Add(*(value for value in right if not constant_term(value)))
        rhs = Add(*(value for value in right if constant_term(value)))-Add(*(value for value in left if constant_term(value)))
    elif multiplicative(lhs) and multiplicative(rhs):
        factors = Mul.make_args(lhs)
        constant = [value for value in factors if constant_term(value)]
        if constant:
            lhs = Mul(*(value for value in factors if not constant_term(value)))
            rhs = rhs/Mul(*constant)
    return Eq(lhs,rhs)


def _factor(expression):
    from ... import Equality, Eq, terms_gcd, Mul
    if isinstance(expression,Equality):
        return Eq(_factor(expression.lhs),_factor(expression.rhs))
    factored = terms_gcd(expression,clear=False,deep=True,expand=False)
    if factored.is_Mul:
        exponential = lambda value:getattr(value.func,'__name__','')=='exp'
        outside = any(exponential(value) for value in factored.args)
        inside = any(exponential(factor) for value in factored.args if value.is_Add
                     for term in value.args for factor in Mul.make_args(term))
        if outside and inside:
            return expression
    return factored


def constantsimp(expr, constants):
    """Remove redundant arbitrary constants without changing their family.

    Correlated appearances of a constant are checked against the entire
    expression before a subexpression is absorbed. Parameter symbols outside
    ``constants`` are retained. Lists and tuples follow the reference API's
    elementwise convention.
    """
    from ... import sympify, Symbol, cse
    from ..._core import Expr
    if isinstance(expr,(list,tuple,set,frozenset)):
        constants = tuple(constants)
        # SymPy's vectorized interface materializes nested sequences as lists.
        return [constantsimp(item,constants) for item in expr]
    expression = sympify(expr)
    if not isinstance(expression,Expr):
        raise TypeError('Constant simplification requires a scalar expression or equation')
    constants = {sympify(value) for value in constants}
    if any(not value.is_Symbol for value in constants):
        raise ValueError('Arbitrary constants must be symbols')
    if not constants or not expression.free_symbols.intersection(constants):
        return expression
    seen = set()
    for _ in range(64):
        if expression in seen:
            return expression
        seen.add(expression)
        original = expression
        for part in _constant_parts(expression,constants):
            selected = sorted(part.free_symbols,key=str)
            if selected and all(expression.count(value)==part.count(value) for value in selected):
                expression = expression.subs(part,selected[0])
        # A repeated nonlinear function of one constant can be renamed once,
        # while a constant also used elsewhere must retain its correlation.
        common, reduced = cse(expression)
        if reduced:
            reduced = reduced[0]
            for temporary,value in reversed(common):
                symbols = value.atoms(Symbol)
                if (len(symbols)==1 and next(iter(symbols)) in constants
                        and not reduced.has(next(iter(symbols)))
                        and not any(other.has(next(iter(symbols))) for name,other in common if name!=temporary)):
                    reduced = reduced.xreplace({temporary:next(iter(symbols))})
                else:
                    reduced = reduced.xreplace({temporary:value})
            expression = reduced
        expression = _factor(_linear_redundancies(expression,constants))
        if expression == original:
            return expression
    raise NotImplementedError('Constant simplification did not converge within 64 transformations')


__all__ = ['constantsimp']
