"""The ``risch``, ``manual`` and ``heurisch`` methods of ``integrate``.

The native integrator runs ``Integral.doit`` and calls :func:`_native_host`
for the stages of the reference ``Integral._eval_integral`` that are ported
to Python: ``risch_integrate`` on the whole integrand, the stages of the term loop
(``ratint`` to ``heurisch``, then ``manualintegrate``) in the default chain,
and
``_eval_integral`` itself when ``manual=True`` or ``heurisch=True`` selects
one method.
"""

from symox import S, Add, DiracDelta, Integral, Piecewise, Poly, Wild, Ne, log


def _native_host(method, f, x, conds):
    """One stage of ``_eval_integral`` for the native integrator.

    ``None`` when the method does not apply, an antiderivative, or for
    ``risch`` a pair ``(elementary, integrand)`` when the Risch algorithm
    proved ``Integral(integrand, x)`` non-elementary.
    """
    if method == 'risch':
        from .risch import risch_integrate
        try:
            result, i = risch_integrate(f, x, separate_integral=True, conds=conds)
        except NotImplementedError:
            return None
        if i:
            return (result, i.function)
        return result
    if method == 'term':
        return _term(f, x, conds)
    if method == 'manual':
        return _manual_stage(f, x, conds)
    if method == 'manual_only':
        return _eval_integral(f, x, 'manual', conds)
    if method == 'heurisch_only':
        return _eval_integral(f, x, 'heurisch', conds)
    raise ValueError(method)


def _linear_power(g, x, conds):
    """The shortcut of the term loop for ``(a*x + b)**n``."""
    if not (g.is_Pow and not g.exp.has(x)):
        return None
    a = Wild('a', exclude=[x])
    b = Wild('b', exclude=[x])
    M = g.base.match(a*x + b)
    if M is None:
        return None
    if g.exp == -1:
        h = log(g.base)
    elif conds != 'piecewise':
        h = g.base**(g.exp + 1) / (g.exp + 1)
    else:
        h1 = log(g.base)
        h2 = g.base**(g.exp + 1) / (g.exp + 1)
        h = Piecewise((h2, Ne(g.exp, -1)), (h1, True))
    return h / M[a]


def _term(g, x, conds):
    """The default term loop on one term ``g`` (its coefficient taken out)
    up to the Meijer G method: the shortcuts, ``ratint``,
    ``trigintegrate``, ``deltaintegrate``, ``singularityintegrate``, then
    ``risch_integrate`` and ``heurisch`` unless they are turned off."""
    from symox import _core
    from .rationaltools import ratint
    from .trigonometry import trigintegrate
    from .deltafunctions import deltaintegrate
    from .singularityfunctions import singularityintegrate
    risch, _, heurisch = _core._integration_methods()
    if g == S.One:
        return x
    h = _linear_power(g, x, conds)
    if h is not None:
        return h
    if g.is_rational_function(x):
        return ratint(g, x)
    # The native delta rules stand in for deltaintegrate: they ran on the
    # whole integrand and decline the singular products (DiracDelta(x)/x,
    # DiracDelta(x)**2, a derivative at an endpoint, ...) that
    # deltaintegrate assigns a value.
    stages = [lambda: trigintegrate(g, x, conds=conds)]
    if not g.has(DiracDelta):
        stages.append(lambda: deltaintegrate(g, x))
    stages.append(lambda: singularityintegrate(g, x))
    for stage in stages:
        h = stage()
        if h is not None:
            return h
    if risch is not False:
        from .risch import risch_integrate, NonElementaryIntegral
        try:
            h, i = risch_integrate(g, x, separate_integral=True, conds=conds)
        except NotImplementedError:
            pass
        else:
            if i:
                rest = i.doit(risch=False)
                if rest.is_Integral and type(rest).__name__ == 'Integral':
                    rest = NonElementaryIntegral(rest.function, *rest.limits)
                h = h + rest
            return h
    if heurisch is not False:
        return _heurisch(g, x, conds)
    return None


def _heurisch(g, x, conds):
    """The ``heurisch`` stage of the term loop."""
    from symox.polys.polyerrors import PolynomialError
    from .heurisch import heurisch, heurisch_wrapper
    try:
        if conds == 'piecewise':
            return heurisch_wrapper(g, x, hints=[])
        return heurisch(g, x, hints=[])
    except PolynomialError:
        return None


def _manualintegrate(g, x):
    """``manualintegrate`` as the term loop takes it: ``None`` for no
    result, an unevaluated integral, ValueError or PolynomialError."""
    from symox.polys.polyerrors import PolynomialError
    from .manualintegrate import manualintegrate
    try:
        result = manualintegrate(g, x)
    except (ValueError, PolynomialError):
        return None
    if result is None or result.is_Integral:
        return None
    return result


def _manual_stage(g, x, conds):
    """The last stage of the default term loop: the integrals left in the
    ``manualintegrate`` result are done with ``manual=False``, and the
    result counts only when none remains."""
    from symox import _core
    result = _manualintegrate(g, x)
    if result is None:
        return None
    if result.has(Integral):
        risch, _, heurisch = _core._integration_methods()
        hints = dict(manual=False, final=False, conds=conds)
        if risch is not None:
            hints['risch'] = risch
        if heurisch is not None:
            hints['heurisch'] = heurisch

        def done(arg):
            if not arg.has(Integral):
                return arg
            return arg.xreplace({i: i.doit(**hints) for i in arg.atoms(Integral)})

        result = result.func(*[done(arg) for arg in result.args]).expand(
            multinomial=False, log=False, power_exp=False, power_base=False)
    if result.has(Integral):
        return None
    return result


def _eval_integral(f, x, method, conds):
    """The reference ``Integral._eval_integral`` with ``manual=True`` or
    ``heurisch=True`` and the other methods off. ``None`` when no
    antiderivative is found."""
    from symox.polys.polyerrors import PolynomialError
    from symox.simplify.fu import sincos_to_sum
    manual = method == 'manual'
    if manual:
        from .manualintegrate import manualintegrate
        try:
            result = manualintegrate(f, x)
            if result is not None and not (result.is_Integral and type(result).__name__ == 'Integral'):
                return result
        except (ValueError, PolynomialError):
            pass
    hints = dict(manual=True) if manual else dict(heurisch=True)
    hints['conds'] = conds
    if isinstance(f, Poly) and not manual:
        return f.integrate(x)
    if isinstance(f, Piecewise):
        return f.piecewise_integrate(x, **hints)
    if not f.has(x):
        return f*x
    if not manual:
        poly = f.as_poly(x)
        if poly is not None:
            return poly.integrate().as_expr()
    parts = []
    args = Add.make_args(f)
    for g in args:
        coeff, g = g.as_independent(x)
        if g == S.One:
            parts.append(coeff*x)
            continue
        order_term = g.getO()
        if order_term is not None:
            h = _eval_integral(g.removeO(), x, method, conds)
            if h is not None:
                h_order = _eval_integral(order_term.expr, x, method, conds)
                if h_order is not None:
                    parts.append(coeff*(h + order_term.func(h_order, *order_term.variables)))
                    continue
            return None
        h = _linear_power(g, x, conds)
        if h is not None:
            parts.append(coeff * h)
            continue
        h = None
        if not manual:
            if g.is_rational_function(x):
                from .rationaltools import ratint
                parts.append(coeff * ratint(g, x))
                continue
            from .trigonometry import trigintegrate
            from .deltafunctions import deltaintegrate
            from .singularityfunctions import singularityintegrate
            h = trigintegrate(g, x, conds=conds)
            if h is not None:
                parts.append(coeff * h)
                continue
            h = deltaintegrate(g, x)
            if h is not None:
                parts.append(coeff * h)
                continue
            h = singularityintegrate(g, x)
            if h is not None:
                parts.append(coeff * h)
                continue
            h = _heurisch(g, x, conds)
        if h is None and manual:
            from .manualintegrate import manualintegrate
            try:
                result = manualintegrate(g, x)
                if result is not None and not result.is_Integral and not result.has(Integral):
                    parts.append(coeff * result)
                    continue
            except (ValueError, PolynomialError):
                pass
        if not h and len(args) == 1:
            f = sincos_to_sum(f).expand(mul=True, deep=False)
            if f.is_Add:
                return _eval_integral(f, x, method, conds)
        if h is not None:
            parts.append(coeff * h)
        else:
            return None
    return Add(*parts)
