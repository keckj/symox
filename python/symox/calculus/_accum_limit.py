"""Limits and definite integrals whose value is an accumulation bounds.

The reference reaches ``limit(sin(x), x, oo) = AccumBounds(-1, 1)`` through
the leading term of ``sin(1/x)`` at 0 (``sin._eval_as_leading_term`` returns
``AccumBounds(-1, 1)`` for an infinite argument) and combines the leading
terms of sums, products and powers with the bounds arithmetic. A definite
integral then subtracts the limits of the antiderivative at the endpoints
(``ExprWithLimits._eval_interval``), so ``integrate(sin(x), (x, 0, oo))`` is
``AccumBounds(0, 2)``.

The native limit engine has no bounds. When it leaves a limit unevaluated,
the expression is evaluated natively term by term
(``symox::limits::accum``): every subexpression without
an oscillating function goes to the native engine, an oscillating function
of an argument tending to an infinity becomes the bounds of its reference
``eval``, and sums, products, powers and functions combine their parts with
the bounds arithmetic. The result is an enclosure of the accumulation
points, as the reference bounds are. Any part left unevaluated, any nan and
any ``oo - oo`` keeps the limit unevaluated.
"""
from .. import _core


def _oscillates(e, z):
    return _core._accum_oscillates(e, z)


def oscillating_limit(e, z, z0, dir):
    """The limit of ``e`` with accumulation bounds, or None (native:
    ``symox::limits::accum::oscillating_limit``)."""
    from .. import sympify
    e, z, z0 = sympify(e), sympify(z), sympify(z0)
    if not isinstance(e, _core.Expr) or not isinstance(z, _core.Expr) \
            or not isinstance(z0, _core.Expr):
        return None
    return _core._oscillating_limit(e, z, z0, str(dir))


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
