"""The ``eval`` branches of the reference functions that return
accumulation bounds (sympy/functions/elementary/trigonometric.py,
exponential.py, integers.py and hyperbolic.py, SymPy 1.14).

``sin(oo)`` and ``cos(oo)`` are ``AccumBounds(-1, 1)``, ``tan(oo)``,
``cot(oo)``, ``sec(oo)`` and ``csc(oo)`` are ``AccumBounds(-oo, oo)``,
``frac(oo)`` is ``AccumBounds(0, 1)``, ``atanh(zoo)`` and ``asech(zoo)`` are
``I*AccumBounds(-pi/2, pi/2)``, and these functions, ``exp`` and ``log`` map
accumulation bounds to accumulation bounds.
"""
from .. import _core

HEADS = frozenset(('sin', 'cos', 'tan', 'cot', 'sec', 'csc', 'exp', 'log', 'frac',
                   'atanh', 'asech', 'atan', 'sinh', 'cosh', 'tanh', 'coth', 'sech', 'csch'))


_MUL_HEADS = frozenset(('exp', 'sinh', 'cosh', 'tanh', 'coth', 'sech', 'csch'))


def _bounds():
    from .accumulationbounds import AccumBounds
    return AccumBounds


def _sin(arg, AB):
    from .. import S, pi, floor, Rational, Min, Max, FiniteSet, sin, oo
    if arg in (S.Infinity, S.NegativeInfinity):
        return AB(-1, 1)
    if isinstance(arg, AB):
        min, max = arg.min, arg.max
        d = floor(min/(2*pi))
        if min is not S.NegativeInfinity:
            min = min - d*2*pi
        if max is not S.Infinity:
            max = max - d*2*pi
        if AB(min, max).intersection(FiniteSet(pi/2, pi*Rational(5, 2))) \
                is not S.EmptySet and \
                AB(min, max).intersection(FiniteSet(pi*Rational(3, 2),
                    pi*Rational(7, 2))) is not S.EmptySet:
            return AB(-1, 1)
        elif AB(min, max).intersection(FiniteSet(pi/2, pi*Rational(5, 2))) \
                is not S.EmptySet:
            return AB(Min(sin(min), sin(max)), 1)
        elif AB(min, max).intersection(FiniteSet(pi*Rational(3, 2), pi*Rational(8, 2))) \
                is not S.EmptySet:
            return AB(-1, Max(sin(min), sin(max)))
        else:
            return AB(Min(sin(min), sin(max)),
                      Max(sin(min), sin(max)))
    return None


def _cos(arg, AB):
    from .. import S, pi, sin
    if arg in (S.Infinity, S.NegativeInfinity):
        return AB(-1, 1)
    if isinstance(arg, AB):
        return sin(arg + pi/2)
    if _real_infinity(arg):
        return AB(-1, 1)
    return None


def _real_infinity(arg):
    # ``arg.is_extended_real and arg.is_finite is False`` of the reference;
    # the kernel leaves both unknown for ``oo*sign(s)`` with a nonzero s.
    from .. import S
    if arg.is_extended_real and arg.is_finite is False:
        return True
    if arg.is_Mul:
        coeff, rest = arg.as_coeff_Mul()
        return (coeff in (S.Infinity, S.NegativeInfinity) and rest.is_extended_real is True
                and rest.is_zero is False)
    return False


def _tan(arg, AB):
    from .. import S, pi, floor, Rational, FiniteSet, tan, oo
    if arg in (S.Infinity, S.NegativeInfinity):
        return AB(-oo, oo)
    if isinstance(arg, AB):
        min, max = arg.min, arg.max
        d = floor(min/pi)
        if min is not S.NegativeInfinity:
            min = min - d*pi
        if max is not S.Infinity:
            max = max - d*pi
        if AB(min, max).intersection(FiniteSet(pi/2, pi*Rational(3, 2))):
            return AB(-oo, oo)
        else:
            return AB(tan(min), tan(max))
    return None


def _cot(arg, AB):
    from .. import S, pi, tan, oo
    if arg in (S.Infinity, S.NegativeInfinity):
        return AB(-oo, oo)
    if isinstance(arg, AB):
        return -tan(arg + pi/2)
    return None


def _reciprocal(of):
    # ReciprocalTrigonometricFunction.eval: 1/t for an evaluated t.
    def evaluate(arg, AB):
        t = of(arg, AB)
        return None if t is None else 1/t
    return evaluate


def _exp(arg, AB):
    from .. import exp, S, re, im
    if isinstance(arg, AB):
        return AB(exp(arg.min), exp(arg.max))
    if arg.is_Mul:
        # exp.eval: an infinite coefficient of a number, exp(I*oo) is nan.
        coeff, terms = arg.as_coeff_Mul()
        if coeff in (S.Infinity, S.NegativeInfinity) and terms.is_number:
            if coeff is S.NegativeInfinity or coeff == S.NegativeInfinity:
                terms = -terms
            if re(terms).is_zero and terms != S.Zero:
                return S.NaN
            if re(terms).is_positive and im(terms) != S.Zero:
                return S.ComplexInfinity
            if re(terms).is_negative:
                return S.Zero
    return None


def _log(arg, AB):
    from .. import S, log, oo
    if isinstance(arg, AB):
        if arg.min.is_positive:
            return AB(log(arg.min), log(arg.max))
        elif arg.min.is_zero:
            return AB(S.NegativeInfinity, log(arg.max))
        return None
    return None


def _frac(arg, AB):
    from .. import S
    if arg in (S.Infinity, S.NegativeInfinity):
        return AB(0, 1)
    return None


def _zoo_half_pi(arg, AB):
    from .. import S, I, pi
    if arg is S.ComplexInfinity or arg == S.ComplexInfinity:
        return I*AB(-pi/2, pi/2)
    return None


def _atan(arg, AB):
    from .. import S, pi
    if arg is S.ComplexInfinity or arg == S.ComplexInfinity:
        return AB(-pi/2, pi/2)
    return None


def _hyperbolic(name):
    # HyperbolicFunction.eval of an imaginary infinity: ``sinh(I*oo)`` is
    # ``I*sin(oo)``, ``cosh(I*oo)`` is ``cos(oo)``, and so on.
    def evaluate(arg, AB):
        from .. import S, I, sin, cos, tan, cot, Mul, oo
        if not arg.is_Mul:
            return None
        coeff = arg/I
        if coeff not in (S.Infinity, S.NegativeInfinity):
            return None
        return {
            'sinh': lambda c: I*sin(c),
            'cosh': lambda c: cos(c),
            'tanh': lambda c: I*tan(c),
            'coth': lambda c: -I*cot(c),
            'sech': lambda c: 1/cos(c),
            # an odd reciprocal: csch(-I*oo) is -csch(I*oo), and 1/(I*sin(oo))
            # keeps I as a factor of the bounds.
            'csch': lambda c: Mul(I, AB(-oo, oo), evaluate=False) if c is S.Infinity
            or c == S.Infinity else Mul(-1, I, AB(-oo, oo), evaluate=False),
        }[name](coeff)
    return evaluate


_EVAL = {
    'sin': _sin, 'cos': _cos, 'tan': _tan, 'cot': _cot,
    'sec': _reciprocal(_cos), 'csc': _reciprocal(_sin),
    'exp': _exp, 'log': _log, 'frac': _frac,
    'atanh': _zoo_half_pi, 'asech': _zoo_half_pi, 'atan': _atan,
    'sinh': _hyperbolic('sinh'), 'cosh': _hyperbolic('cosh'), 'tanh': _hyperbolic('tanh'),
    'coth': _hyperbolic('coth'), 'sech': _hyperbolic('sech'), 'csch': _hyperbolic('csch'),
}


def evaluate(name, arg):
    """The reference value of ``name(arg)`` when it is (or involves)
    accumulation bounds, else None."""
    if not isinstance(arg, _core.Expr):
        return None
    AB = _bounds()
    # Only infinite arguments and the bounds themselves reach the branches.
    if not isinstance(arg, AB) and not arg.is_infinite and arg.is_finite is not False:
        if not arg.is_Mul or not arg.args or not arg.args[0].is_infinite:
            return None
        if name not in _MUL_HEADS and not (name in ('cos', 'sec') and _real_infinity(arg)):
            return None
    return _EVAL[name](arg, AB)


def _illegal_values():
    from .. import S
    return (S.Infinity, S.NegativeInfinity, S.ComplexInfinity)


def _needs(e):
    """Whether ``e`` holds an application the reference evaluates to
    accumulation bounds (or one of their functions)."""
    stack, seen = [e], set()
    while stack:
        node = stack.pop()
        args = getattr(node, 'args', ())
        if not args or node in seen:
            continue
        seen.add(node)
        if len(args) == 1 and getattr(node, 'is_Function', False):
            name = getattr(node.func, '__name__', None)
            if name in HEADS and evaluate(name, args[0]) is not None:
                return True
        if node.is_Pow or node.is_Add or node.is_Mul:
            from .accumulationbounds import AccumBounds
            if any(isinstance(a, AccumBounds) for a in args):
                return True
        stack.extend(a for a in args if isinstance(a, _core.Expr))
    return False


def reevaluate(e, cache=None):
    """Rebuild ``e`` bottom-up through the public constructors and the
    bounds arithmetic, as the reference rebuilds after a substitution."""
    if cache is None:
        cache = {}
    try:
        return cache[e]
    except (KeyError, TypeError):
        pass
    value = _reevaluate(e, cache)
    try:
        cache[e] = value
    except TypeError:
        pass
    return value


def _reevaluate(e, cache):
    import operator
    from functools import reduce
    import symox
    from .accumulationbounds import AccumBounds
    args = getattr(e, 'args', ())
    if not args or not isinstance(e, _core.Expr) or isinstance(e, AccumBounds):
        return e
    new = [reevaluate(a, cache) if isinstance(a, _core.Expr) else a for a in args]
    if all(a is b for a, b in zip(new, args)):
        if not (len(args) == 1 and getattr(e, 'is_Function', False)
                and getattr(e.func, '__name__', None) in HEADS):
            if not ((e.is_Add or e.is_Mul or e.is_Pow) and any(isinstance(a, AccumBounds) for a in args)):
                return e
    if e.is_Add:
        return reduce(operator.add, new)
    if e.is_Mul:
        return reduce(operator.mul, new)
    if e.is_Pow:
        return new[0]**new[1]
    if getattr(e, 'is_Function', False):
        name = getattr(e.func, '__name__', None)
        head = getattr(symox, name, None) if name in HEADS else None
        if head is not None:
            return head(*new)
    return e.func(*new)


def _substitutes_illegal(args):
    from .accumulationbounds import AccumBounds
    illegal = _illegal_values()
    if len(args) == 2:
        values = [args[1]]
    elif len(args) == 1:
        sequence = args[0]
        if isinstance(sequence, dict):
            values = list(sequence.values())
        else:
            try:
                values = [pair[1] for pair in sequence]
            except (TypeError, IndexError, KeyError):
                return False
    else:
        return False
    return any(isinstance(v, AccumBounds) or (isinstance(v, _core.Expr) and v in illegal)
               for v in values)


def install_subs():
    """``subs`` of an infinity or of bounds evaluates the functions that the
    reference ``eval`` maps to accumulation bounds: ``sin(x).subs(x, oo)``
    is ``AccumBounds(-1, 1)``."""
    from functools import wraps
    from .accumulationbounds import _CREATED
    original = _core.Expr.subs

    @wraps(original)
    def subs(self, *args, **kwargs):
        result = original(self, *args, **kwargs)
        if (args and isinstance(result, _core.Expr) and (_CREATED[0] or _substitutes_illegal(args))
                and _needs(result)):
            return reevaluate(result)
        return result
    _core.Expr.subs = subs


def install_doit():
    """``Limit(...).doit()`` and ``Integral(...).doit()`` take the bounds
    route of ``limit`` and ``integrate`` when the native evaluation keeps
    them held."""
    from functools import wraps
    original = _core.Expr.doit

    @wraps(original)
    def doit(self, *args, **kwargs):
        result = original(self, *args, **kwargs)
        if getattr(result, 'is_Limit', False):
            return _limit_doit(result)
        if getattr(result, 'is_Integral', False):
            return _integral_doit(result)
        return result
    _core.Expr.doit = doit


def _limit_doit(held):
    from ._accum_limit import oscillating_limit
    try:
        e, z, z0, direction = held.args
        value = oscillating_limit(e, z, z0, str(direction))
    except (TypeError, ValueError, NotImplementedError, AttributeError):
        return held
    return held if value is None else value


def _integral_doit(held):
    from ._accum_limit import oscillating_integral
    from .. import _native_integrate
    try:
        limits = held.limits
        if len(limits) != 1 or len(limits[0]) != 3:
            return held
        value = oscillating_integral(_native_integrate, None, held.function, *limits[0])
    except (TypeError, ValueError, NotImplementedError, AttributeError):
        return held
    return held if value is None else value
