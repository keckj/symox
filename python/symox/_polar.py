"""Polar surface conversions using native expressions and explicit winding."""
from . import _core


def _s():
    import symox
    return symox


def polar_lift(arg):
    return _core._apply('polar_lift', arg)


def periodic_argument(ar, period):
    return _core._apply('periodic_argument', ar, period)


def principal_branch(x, period):
    return _core._apply('principal_branch', x, period)


def unbranched_argument(arg):
    return periodic_argument(arg, _s().oo)


for _head in (polar_lift, periodic_argument, principal_branch):
    _head._symox_head_name = _head.__name__


def polarify(eq, subs=True, lift=False):
    """Lift numbers to the logarithm surface; optionally replace symbols.

    With subs=True return (expression, reverse_substitution). With lift=True
    lift symbols and sums explicitly and return just the expression.
    """
    s = _s()
    expr = s.sympify(eq)
    def visit(value, pause=False):
        if value.is_polar is True:
            return value
        if not pause and (value.is_number or (lift and value.is_Symbol)):
            return polar_lift(value)
        if not value.args:
            return value
        if value.is_Add:
            result = value.func(*(visit(a, True) for a in value.args))
            return polar_lift(result) if lift else result
        if value.class_name in ('Integral', 'Sum', 'Product'):
            limits = [tuple([bound.args[0]]+[visit(a) for a in bound.args[1:]]) for bound in value.args[1:]]
            return value.func(visit(value.args[0]), *limits)
        if value.is_Function or (value.is_Pow and value.args[0] == s.E):
            pause = False
        return value.func(*(visit(a, pause) for a in value.args))
    result = visit(expr)
    if not subs or lift:
        return result
    mapping = {symbol: s.Dummy(symbol.name, polar=True) for symbol in result.free_symbols}
    return result.xreplace(mapping), {v:k for k,v in mapping.items()}


# Heads SymPy marks ``unbranched``: re, im, Abs, ExpBase, the trigonometric,
# hyperbolic, airy, error, Fresnel, gamma, beta and Mathieu families.
_UNBRANCHED = {'exp','exp_polar','sin','cos','tan','cot','sec','csc','sinc',
               'sinh','cosh','tanh','coth','sech','csch','gamma','multigamma',
               'erf','erfc','erfi','fresnels','fresnelc','Abs','re','im',
               'airyai','airybi','airyaiprime','airybiprime','beta','betainc',
               'betainc_regularized','mathieus','mathieuc','mathieusprime',
               'mathieucprime'}


def _projects_children(value):
    """Add, Mul, Boolean and the relations SymPy's _unpolarify recurses into."""
    if value.is_Add or value.is_Mul or getattr(value, 'is_Boolean', False):
        return True
    if getattr(value, 'is_Relational', False):
        if value.rel_op in ('==', '!='):
            return any(a == 0 for a in value.args)
        return True
    return False


def unpolarify(eq, subs=None, exponents_only=False):
    """Project polar factors only where the surrounding operation is unbranched.

    A polar base of a noninteger power or argument of log stays lifted.
    exponents_only restricts projection to exponents and unbranched arguments.
    """
    s = _s()
    if isinstance(eq, bool):
        return eq
    expr = s.sympify(eq)
    if subs is not None:
        # The reference recurses as unpolarify(eq.subs(subs)): the
        # substitution drops exponents_only.
        expr = expr.subs(subs)
        exponents_only = False
    def visit(value, pause):
        name = value.class_name
        if not pause:
            if name == 'exp_polar':
                return s.exp(visit(value.args[0], False))
            if name == 'polar_lift':
                return visit(value.args[0], False)
            if name == 'principal_branch' and value.args[1] == 2*s.pi:
                return visit(value.args[0], False)
        if not value.args:
            return value
        if value.is_Pow:
            base, exponent = value.args
            exponent = visit(exponent, False)
            return value.func(visit(base, pause or exponent.is_integer is not True), exponent)
        if name in _UNBRANCHED:
            child_pause = bool(exponents_only)
        elif not pause and _projects_children(value):
            child_pause = False
        else:
            child_pause = True
        return value.func(*(visit(a, child_pause) for a in value.args))
    while True:
        updated = visit(expr, bool(exponents_only))
        if updated == expr:
            break
        expr = updated
    return expr.xreplace({s.exp_polar(0):s.S.One, polar_lift(0):s.S.Zero})


__all__ = ['polar_lift','periodic_argument','principal_branch',
           'unbranched_argument','polarify','unpolarify']
