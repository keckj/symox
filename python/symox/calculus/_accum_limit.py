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
        if _is_bounds(base) and exp == S.Zero and e.exp.has(z):
            # The reference leading term keeps the vanishing exponent:
            # AccumBounds(a, b)**t with t -> 0+ has the bounds of a**t and
            # b**t, where a zero endpoint stays 0 (oo for t -> 0-), an infinite
            # endpoint stays oo (0 for t -> 0-), so
            # frac(x)**(1/x) at oo is AccumBounds(0, 1), not 1.
            side = _native(1/e.exp, z, z0, dir)
            if side not in (S.Infinity, S.NegativeInfinity) or base.min.is_nonnegative is not True:
                raise _Fail
            pos = side is S.Infinity
            ends = [(S.Zero if pos else S.Infinity) if v == S.Zero
                    else (S.Infinity if pos else S.Zero) if v is S.Infinity
                    else S.One
                    for v in (base.min, base.max)]
            from .accumulationbounds import AccumBounds
            return _bounded_value(AccumBounds(min(ends), max(ends)))
        if _is_bounds(exp) or (_is_bounds(base) and exp.is_rational is not True
                               and base.min.is_nonnegative is not True):
            # Admission only: a base that may be negative under an
            # irrational, infinite or oscillating exponent is not real, and
            # the reference keeps such limits unevaluated (issues 14276,
            # 18442: sin(x)**log(x), sin(x)**cos(x), tan(x)**(2**sqrt(pi))).
            # A rational exponent stays a bounds power (sqrt(AccumBounds(-1,
            # 1)), whose product with a vanishing factor is 0).
            raise _Fail
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
        if _held_oscillation(arg):
            # exp(sin(oo*sign(s))): the reference leaves it unevaluated
            raise _Fail
        value = head(arg)
        if _contains_bounds(value):
            used.append(True)
        if name in _OSCILLATING and _held_oscillation(value):
            # Left unevaluated, e.g. sin(oo*sign(s)). The reference
            # heuristics keep it when the argument limit is known to be
            # infinite (a nonzero s): sin(oo*sign(s)) stands for
            # AccumBounds(-1, 1) once the sign is known. An argument that
            # may be finite (oo*sign(s) for a possibly zero s) is undecided.
            if arg.is_finite is False and _real_infinity_times_sign(arg):
                used.append(True)
                return value
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
    if not isinstance(e, _core.Expr):
        return None
    if not _oscillates(e, z):
        return _bounded_affine(e, z, z0, dir)
    # The native decomposition (affine in bounded sin/cos values, their
    # integer powers and rounding remainders) is tried first; the walker
    # below only takes the forms it declines.
    native = _bounded_affine(e, z, z0, dir)
    if native is not None:
        return native
    if dir == '+-':
        if z0.is_infinite:
            return None
        right, left = oscillating_limit(e, z, z0, '+'), oscillating_limit(e, z, z0, '-')
        return right if right is not None and right == left else None
    used = []
    try:
        value = _walk(e, z, z0, dir, used)
    except (_Fail, TypeError, ValueError, NotImplementedError):
        return _bounded_affine(e, z, z0, dir)
    if not used or not isinstance(value, _core.Expr):
        return None
    if value.has(S.NaN) or value == S.ComplexInfinity:
        return None
    return value


def _held_oscillation(value):
    """Whether ``value`` holds an oscillating function left unevaluated at an
    infinite argument, such as ``sin(oo*sign(s))``."""
    from .. import S
    for node in _nodes(value):
        func = getattr(node, 'func', None)
        if getattr(func, '__name__', None) in _OSCILLATING and getattr(node, 'args', None) \
                and node.args[0].has(S.Infinity, S.NegativeInfinity):
            return True
    return False


def _real_infinity_times_sign(arg):
    """``oo*sign(s)`` (up to a sign): the reference heuristics value of the
    argument ``s*x`` at an infinity."""
    from .. import S, sign
    for factor in arg.as_ordered_factors():
        if factor in (S.Infinity, S.NegativeInfinity, S.NegativeOne):
            continue
        if not isinstance(factor, sign):
            return False
    return True


def _bounded_affine(e, z, z0, dir):
    """The limit of an expression affine in oscillating sin/cos parts once
    the tractable rewrites of Si, Ci and the Fresnel integrals are made:
    ``c0 + sum(c_i*AccumBounds(-1, 1))`` as the reference leading terms
    give (``x*Ci(x)`` at ``oo`` is ``AccumBounds(-1, 1)``); a floor, ceiling
    or Mod of an argument tending to an infinity contributes the bounds
    ``AccumBounds(-1, 0)`` or ``AccumBounds(0, 1)`` of its remainder. None when the
    native decomposition does not apply."""
    from .accumulationbounds import AccumBounds
    if dir not in ('+', '-'):
        return None
    try:
        parts = _core._limit_bounded_affine(e, z, z0, dir)
    except (TypeError, ValueError, NotImplementedError):
        return None
    if parts is None:
        return None
    from .. import Add, Mul
    c0, coeffs = parts
    # The reference leading terms combine through Mul and Add flatten:
    # x*(fresnels(x) - 1/2) is AccumBounds(-1, 1)/pi, and the floor remainder
    # AccumBounds(-1, 0) of floor(x + 1/2) - floor(x) gives 1/2 + AccumBounds(-1, 1).
    value = Add(c0, *[Mul(c, AccumBounds(lo, hi)) for c, lo, hi in coeffs])
    return _checked(value)


def _checked(value):
    from .. import S
    if not isinstance(value, _core.Expr) or value.has(S.NaN) or value == S.ComplexInfinity:
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
