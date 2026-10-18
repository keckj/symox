"""Native scalar analysis primitives shared by real calculus adapters."""
from .. import _core

MAX_NODES=10_000


def scalar(value):
    from .. import sympify, Set, _wrap_set
    from .._structural_queries import _is_expr
    value=sympify(value)
    if not isinstance(value,_core.Expr):
        raise TypeError('Scalar analysis requires a native scalar expression')
    if getattr(value.func,'__name__','')=='Tuple':
        raise TypeError('Scalar analysis does not accept tuple expressions')
    if not _is_expr(value) or isinstance(_wrap_set(value),Set):
        raise TypeError('Scalar analysis requires an arithmetic expression')
    return value


def real_constants(expression, symbol):
    """Require independent coefficient chunks to have a known real value."""
    pending=[expression]
    seen=set()
    while pending:
        node=pending.pop()
        if not isinstance(node,_core.Expr) or node in seen:continue
        seen.add(node)
        if len(seen)>MAX_NODES:raise NotImplementedError('Real coefficient checks exceed 10000 nodes')
        if not node.has(symbol):
            if node.is_real is not True:
                raise NotImplementedError('Independent coefficients must be known to be real')
        else:pending.extend(node.args)


def variable(value):
    value=scalar(value)
    if not value.is_Symbol:
        raise ValueError('Analysis variable must be a symbol')
    return value


def walk(expression):
    pending=[expression]
    seen=set()
    while pending:
        node=pending.pop()
        if not isinstance(node,_core.Expr) or node in seen:
            continue
        seen.add(node)
        if len(seen)>MAX_NODES:
            raise NotImplementedError('Scalar analysis exceeds 10000 distinct nodes')
        yield node
        pending.extend(node.args)


def name(expression):
    return '' if expression.is_UndefinedFunction else getattr(expression.func,'__name__','')


def zeros(expression, symbol, domain):
    from .. import solveset
    return solveset(expression,symbol,domain=domain)


def real_domain(domain):
    from .. import S, Set
    if not isinstance(domain,Set):
        raise TypeError('Domain must be a native symbolic set')
    if domain.is_subset(S.Reals) is not True:
        raise NotImplementedError('Real analysis requires a domain known to be real')
    return domain


def inequality(condition, symbol, domain):
    from .. import solve_univariate_inequality
    return solve_univariate_inequality(condition,symbol,relational=False,domain=domain)
