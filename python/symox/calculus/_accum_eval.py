"""The ``eval`` branches of the reference functions that return
accumulation bounds (sympy/functions/elementary/trigonometric.py,
exponential.py, integers.py and hyperbolic.py, SymPy 1.14).

``sin(oo)`` and ``cos(oo)`` are ``AccumBounds(-1, 1)``, ``tan(oo)``,
``cot(oo)``, ``sec(oo)`` and ``csc(oo)`` are ``AccumBounds(-oo, oo)``,
``frac(oo)`` is ``AccumBounds(0, 1)``, ``atanh(zoo)`` and ``asech(zoo)`` are
``I*AccumBounds(-pi/2, pi/2)``, and these functions, ``exp`` and ``log`` map
accumulation bounds to accumulation bounds. The branches are native
(``symox::limits::accum``); this module keeps the admission and the
re-evaluation of substituted expressions.
"""
from .. import _core

HEADS = frozenset(('sin', 'cos', 'tan', 'cot', 'sec', 'csc', 'exp', 'log', 'frac',
                   'atanh', 'asech', 'atan', 'sinh', 'cosh', 'tanh', 'coth', 'sech', 'csch'))



def evaluate(name, arg):
    """The reference value of ``name(arg)`` when it is (or involves)
    accumulation bounds, else None. The branches are native
    (``symox::limits::accum::evaluate``)."""
    if not isinstance(arg, _core.Expr) or name not in HEADS:
        return None
    return _core._accum_evaluate(name, arg)


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
    """The oscillating route for a limit the native ``doit`` keeps held; if
    it does not answer either, the reference ``NotImplementedError`` of the
    native decision is raised (as from ``Limit.doit``)."""
    from ._accum_limit import oscillating_limit
    try:
        e, z, z0, direction = held.args
    except (TypeError, ValueError):
        return held
    try:
        value = oscillating_limit(e, z, z0, str(direction))
    except (TypeError, ValueError, NotImplementedError, AttributeError):
        value = None
    if value is not None:
        return value
    try:
        _core._limit(e, z, z0, str(direction), True)
    except NotImplementedError:
        raise
    except (TypeError, ValueError, AttributeError):
        pass
    return held


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
