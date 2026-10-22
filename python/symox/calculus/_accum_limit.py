"""Limits and definite integrals whose value is an accumulation bounds.

The reference reaches ``limit(sin(x), x, oo) = AccumBounds(-1, 1)`` through
the leading term of ``sin(1/x)`` at 0 (``sin._eval_as_leading_term`` returns
``AccumBounds(-1, 1)`` for an infinite argument) and combines the leading
terms of sums, products and powers with the bounds arithmetic. A definite
integral then subtracts the limits of the antiderivative at the endpoints
(``ExprWithLimits._eval_interval``), so ``integrate(sin(x), (x, 0, oo))`` is
``AccumBounds(0, 2)``.

The native limit engine has no bounds. When it leaves a limit unevaluated,
the expression is evaluated here term by term: every subexpression without
an oscillating function goes to the native engine, an oscillating function
of an argument tending to an infinity becomes the bounds of its reference
``eval``, and sums, products, powers and functions combine their parts with
the bounds arithmetic. The result is an enclosure of the accumulation
points, as the reference bounds are. Any part left unevaluated, any nan and
any ``oo - oo`` keeps the limit unevaluated.
"""
from .. import _core
from ._accum_eval import _real_infinity

_OSCILLATING = ('sin', 'cos', 'tan', 'cot', 'sec', 'csc', 'frac')


class _Fail(Exception):
    pass


def _oscillates(e, z):
    from .. import sin, cos, tan, cot, sec, csc, frac
    return any(f.has(z) for f in e.atoms(sin, cos, tan, cot, sec, csc, frac))


def _bounded_value(value):
    from .. import S
    from .accumulationbounds import AccumBounds
    if value is None or value == S.NaN or value.has(S.NaN):
        raise _Fail
    if any(getattr(node, 'is_Limit', False) for node in _nodes(value)):
        raise _Fail
    return value


def _native(e, z, z0, dir):
    r = _core._limit(e, z, z0, dir)
    if r is None:
        raise _Fail
    return _bounded_value(r)


def _number_like(v):
    from .accumulationbounds import AccumBounds
    return isinstance(v, AccumBounds) or v.is_number and v.is_extended_real is True


def _combine(values, op, rebuild):
    """Combine numbers and bounds with the bounds arithmetic, as the
    reference flatten does for numeric terms; symbolic parts stay in the
    native constructor ``rebuild``."""
    from functools import reduce
    numeric = [v for v in values if _number_like(v)]
    symbolic = [v for v in values if not _number_like(v)]
    if not numeric:
        return rebuild(*symbolic)
    combined = reduce(op, numeric)
    if not symbolic:
        return combined
    return op(combined, rebuild(*symbolic)) if _number_like(combined) and not \
        _is_bounds(combined) else rebuild(combined, *symbolic)


def _is_bounds(v):
    from .accumulationbounds import AccumBounds
    return isinstance(v, AccumBounds)


def _walk(e, z, z0, dir, used):
    from .. import S, Add, Mul
    import operator
    if not e.has(z):
        return e
    if e == z:
        return z0
    if not _oscillates(e, z):
        return _native(e, z, z0, dir)
    if e.is_Add:
        values = [_walk(a, z, z0, dir, used) for a in e.args]
        if S.Infinity in values and S.NegativeInfinity in values:
            raise _Fail
        return _bounded_value(_combine(values, operator.add, Add))
    if e.is_Mul:
        values = [_walk(a, z, z0, dir, used) for a in e.args]
        return _bounded_value(_combine(values, operator.mul, Mul))
    if e.is_Pow:
        base = _walk(e.base, z, z0, dir, used)
        exp = _walk(e.exp, z, z0, dir, used)
        return _bounded_value(base**exp)
    if e.is_Function and len(e.args) == 1:
        name = e.func.__name__
        arg = _walk(e.args[0], z, z0, dir, used)
        head = _function(name)
        if head is None:
            raise _Fail
        if name in _OSCILLATING and not _is_bounds(arg) and arg.is_infinite is not True \
                and arg.is_finite is not False and not _real_infinity(arg):
            # A finite limit of the argument: the native engine knows
            # better (continuity at the point), unless it was bounds.
            return _native(e, z, z0, dir)
        value = head(arg)
        if _contains_bounds(value):
            used.append(True)
        if name in _OSCILLATING and value.func.__name__ == name:
            # left unevaluated, e.g. sin(oo*sign(s))
            raise _Fail
        return _bounded_value(value)
    raise _Fail


def _contains_bounds(value):
    from .accumulationbounds import AccumBounds
    return any(isinstance(node, AccumBounds) for node in _nodes(value))


def _function(name):
    import symox
    head = getattr(symox, name, None)
    return head if callable(head) else None


def oscillating_limit(e, z, z0, dir):
    """The limit of ``e`` with accumulation bounds, or None."""
    from .. import S, sympify
    e, z, z0 = sympify(e), sympify(z), sympify(z0)
    if dir not in ('+', '-', '+-') or not z.is_Symbol or z0 is S.ComplexInfinity \
            or z0 == S.ComplexInfinity or z0 == S.NaN:
        return None
    if not isinstance(e, _core.Expr) or not _oscillates(e, z):
        return None
    if dir == '+-':
        if z0.is_infinite:
            return None
        right, left = oscillating_limit(e, z, z0, '+'), oscillating_limit(e, z, z0, '-')
        return right if right is not None and right == left else None
    used = []
    try:
        value = _walk(e, z, z0, dir, used)
    except (_Fail, TypeError, ValueError, NotImplementedError):
        return None
    if not used or not isinstance(value, _core.Expr):
        return None
    if value.has(S.NaN) or value == S.ComplexInfinity:
        return None
    return value


def oscillating_integral(native_integrate, limit, f, x, a, b):
    """``integrate(f, (x, a, b))`` from the limits of the antiderivative at
    infinite endpoints, as ``ExprWithLimits._eval_interval``, when one of
    them is an accumulation bounds; else None."""
    from .. import S, sympify
    f, x, a, b = (sympify(v) for v in (f, x, a, b))
    infinite = (S.Infinity, S.NegativeInfinity)
    if not x.is_Symbol or (a not in infinite and b not in infinite):
        return None
    if not _oscillates(f, x) or a.has(x) or b.has(x):
        return None
    F = native_integrate(f, x)
    if not isinstance(F, _core.Expr) or _unsafe(F):
        return None
    values = []
    used = False
    for point, side in ((b, '-'), (a, '+')):
        if point in infinite:
            value = _core._limit(F, x, point, side)
            if value is None:
                value = oscillating_limit(F, x, point, side)
                if value is None:
                    return None
                used = True
        else:
            value = F.subs(x, point)
        if value.has(S.NaN) or value.has(S.ComplexInfinity) or not isinstance(value, _core.Expr):
            return None
        values.append(value)
    if not used:
        return None
    result = values[0] - values[1]
    if result.has(S.NaN):
        return None
    return result


def _unsafe(F):
    """Antiderivatives whose endpoint limits the bounds route cannot take:
    held integrals, Piecewise and parameter conditions."""
    return any(getattr(node, 'is_Integral', False) or getattr(node, 'is_Piecewise', False)
               or getattr(node, 'is_Limit', False) for node in _nodes(F))


def _nodes(e):
    stack, seen = [e], set()
    while stack:
        node = stack.pop()
        try:
            if node in seen:
                continue
            seen.add(node)
        except TypeError:
            pass
        yield node
        stack.extend(getattr(node, 'args', ()))
