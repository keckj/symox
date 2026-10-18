"""Shared bounded differential-equation inspection, without foreign CAS calls."""


def nodes(expression):
    from .._structural_queries import _graph
    _, records, order = _graph(expression)
    if len(records) > 10000:
        raise NotImplementedError('Differential analysis exceeds 10000 graph nodes')
    return [records[key][0] for key in order]


def equation(eq):
    from .. import sympify, Equality
    from ..calculus._analysis import scalar
    eq = sympify(eq)
    return scalar(eq.lhs-eq.rhs if isinstance(eq, Equality) else eq)


def dependent(expression, func=None, dimension=1):
    from .. import sympify
    if func is None:
        candidates = {node for node in nodes(expression)
                      if getattr(node, 'is_UndefinedFunction', False)
                      and len(node.args) == dimension and all(arg.is_Symbol for arg in node.args)}
        # Prefer functions occurring under a derivative; a forcing term may
        # itself be an undefined function of the same independent variables.
        differentiated = {node.expr for node in nodes(expression)
                          if node.is_Derivative and node.expr in candidates}
        candidates = differentiated or candidates
        if len(candidates) != 1:
            raise ValueError('Specify exactly one dependent function')
        func = candidates.pop()
    else:
        func = sympify(func)
    if (not getattr(func, 'is_UndefinedFunction', False) or len(func.args) != dimension
            or not all(arg.is_Symbol for arg in func.args) or len(set(func.args)) != dimension):
        raise ValueError('Expected an applied undefined function of distinct independent symbols')
    return func


def replace_function(expression, func, value):
    from .. import diff
    from .deutils import ode_order
    if ode_order(expression, func) > 256:
        raise NotImplementedError('Differential order exceeds 256')
    replacements = {func: value}
    for node in nodes(expression):
        if node.is_Derivative and node.expr == func:
            result = value
            for variable, count in node.variable_count:
                result = diff(result, variable, int(count))
            replacements[node] = result
    return expression.xreplace(replacements)


def linear_data(expression, coordinates, clear_denominators=True):
    from .. import Dummy, expand, cancel, S
    placeholders = [Dummy('differential_coordinate') for _ in coordinates]
    expression = expression.xreplace(dict(zip(coordinates, placeholders)))
    numerator = cancel(expression)
    if clear_denominators:
        numerator = numerator.as_numer_denom()[0]
    numerator = expand(numerator)
    coefficients = [cancel(numerator.coeff(symbol)) for symbol in placeholders]
    if any(coefficient.has(*placeholders) for coefficient in coefficients):
        return None
    rest = cancel(expand(numerator-sum((c*s for c, s in zip(coefficients, placeholders)), S.Zero)))
    if rest.has(*placeholders):
        return None
    return coefficients, rest


def fresh_constant(expression):
    from .. import Symbol
    names = {symbol.name for symbol in expression.free_symbols}
    index = 1
    while 'C'+str(index) in names:
        index += 1
    return Symbol('C'+str(index))
