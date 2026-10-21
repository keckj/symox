"""Periodic branch of ``solve_univariate_inequality`` (SymPy 1.14).

SymPy solves an inequality of a periodic expression over ``S.Reals`` on one
period ``[0, T)`` only (the convention of ``test_trig_inequalities``: "all
the inequalities are solved in a periodic interval").  For another domain it
intersects that one-period answer with the domain, which drops solutions
outside ``[0, T)`` (``sin(x) > 0`` on ``[2*pi, 3*pi]`` gives ``EmptySet``).
Here the one-period answer is used when it is exact on the domain; on any
other bounded domain the inequality is solved on the domain itself, and on
a half-unbounded domain it is refused (``ConditionSet``) instead of being
decided wrongly.

``window`` is the integration interval when the condition comes from
``Piecewise._intervals`` during a definite integral, ``False`` when the
caller needs a condition valid on the whole real line (an antiderivative,
``piecewise_simplify``), and ``None`` for a direct call.
"""


def _pt(start, end):
    from .. import S
    if not start.is_infinite and not end.is_infinite:
        return (start + end)/2
    if start.is_infinite and end.is_infinite:
        return S.Zero
    if (start.is_infinite and start.is_extended_positive is None or
            end.is_infinite and end.is_extended_positive is None):
        raise ValueError('cannot proceed with unsigned infinite values')
    if (end.is_infinite and end.is_extended_negative or
            start.is_infinite and start.is_extended_positive):
        start, end = end, start
    if end.is_infinite:
        if start.is_extended_positive:
            return start*2
        if start.is_extended_negative:
            return start*S.Half
        return start + 1
    if end.is_extended_positive:
        return end*S.Half
    if end.is_extended_negative:
        return end*2
    return end - 1


def _nsort_reals(points):
    """``_nsort(points, separated=True)[0]`` for exact real numbers."""
    from .. import S
    values = []
    for p in points:
        v = p.evalf(50)
        re, im = v.as_real_imag()
        if im != 0 and im.is_zero is not True:
            continue
        values.append((re, p))
    values.sort(key=lambda t: t[0])
    for (a, _), (b, _) in zip(values, values[1:]):
        if (b - a).is_zero is True or b == a:
            raise NotImplementedError('sorting of these roots is not supported')
    return [p for _, p in values]


def _cannot(expr, gen):
    from .. import Symbol
    from .._public_tools import filldedent
    return NotImplementedError(filldedent('''
        The inequality, %s, cannot be solved using
        solve_univariate_inequality.
        ''' % expr.subs(gen, Symbol('x'))))


def _solve_on(expr, e, gen, domain, continuous=False):
    """The general part of ``solve_univariate_inequality`` on ``domain``."""
    from .. import (S, FiniteSet, Interval, Union, Intersection, expand_mul,
                    ConditionSet)
    from ..calculus.util import continuous_domain
    from .solveset import solvify
    from ._algebraic import denoms
    _domain = domain
    n, d = e.as_numer_denom()
    try:
        if gen not in n.free_symbols and len(e.free_symbols) > 1:
            raise ValueError
        solns = solvify(e, gen, domain)
        if solns is None:
            raise ValueError
    except (ValueError, NotImplementedError):
        raise _cannot(expr, gen)
    expanded_e = expand_mul(e)

    def valid(x):
        v = expanded_e.subs(gen, expand_mul(x))
        try:
            r = expr.func(v, 0)
        except TypeError:
            r = S.false
        if r in (S.true, S.false):
            return r
        if v.is_extended_real is False:
            return S.false
        v = v.n(2)
        if v.is_comparable:
            return expr.func(v, 0)
        raise NotImplementedError('relationship did not evaluate: %s' % r)

    singularities = []
    for den in denoms(expr, gen):
        found = solvify(den, gen, domain)
        if found is None:
            raise _cannot(expr, gen)
        singularities.extend(found)
    if not continuous:
        domain = continuous_domain(expanded_e, gen, domain)
    include_x = '=' in expr.rel_op and expr.rel_op != '!='
    try:
        discontinuities = set(domain.boundary - FiniteSet(domain.inf, domain.sup))
        critical_points = FiniteSet(*(list(solns) + singularities + list(discontinuities))).intersection(
            Interval(domain.inf, domain.sup, domain.inf not in domain, domain.sup not in domain))
        if isinstance(critical_points, ConditionSet) or not critical_points.is_FiniteSet:
            raise NotImplementedError
        if not all(r.is_number for r in critical_points):
            raise NotImplementedError
        reals = _nsort_reals(list(critical_points))
    except (NotImplementedError, TypeError):
        raise NotImplementedError('sorting of these roots is not supported')
    if expanded_e.coeff(S.ImaginaryUnit) != S.Zero:
        # the imaginary-coefficient branch of the reference is not ported
        raise _cannot(expr, gen)
    sol_sets = [S.EmptySet]
    start = domain.inf
    if start in domain and valid(start) and start.is_finite:
        sol_sets.append(FiniteSet(start))
    for x in reals:
        end = x
        if valid(_pt(start, end)):
            sol_sets.append(Interval(start, end, True, True))
        if x in singularities:
            singularities.remove(x)
        else:
            if x in discontinuities:
                discontinuities.remove(x)
                _valid = valid(x)
            else:
                _valid = include_x
            if _valid:
                sol_sets.append(FiniteSet(x))
        start = end
    end = domain.sup
    if end in domain and valid(end) and end.is_finite:
        sol_sets.append(FiniteSet(end))
    if valid(_pt(start, end)):
        sol_sets.append(Interval.open(start, end))
    return Intersection(Union(*sol_sets), S.Reals, _domain)


def _range_decides(expr, e, gen, domain):
    from .. import S
    from ..calculus.util import function_range
    try:
        frange = function_range(e, gen, domain)
    except (NotImplementedError, ValueError):
        return None
    rel = expr.rel_op
    if rel in ('<', '<='):
        if expr.func(frange.sup, 0) == True:
            return domain
        if expr.func(frange.inf, 0) == False:
            return S.EmptySet
    elif rel in ('>', '>='):
        if expr.func(frange.inf, 0) == True:
            return domain
        if expr.func(frange.sup, 0) == False:
            return S.EmptySet
    return None


def _one_period(expr, e, gen, period, continuous):
    """The reference answer over ``S.Reals``: the solutions in ``[0, T)``."""
    from .. import S, Interval
    rv = _range_decides(expr, e, gen, S.Reals)
    if rv is not None:
        return rv
    return _solve_on(expr, e, gen, Interval(0, period, False, True), continuous)


def native_solve(expr, gen, domain, window=None):
    """Entry point of the native host: core expressions in and out."""
    from .. import _wrap_set, _unwrap
    if isinstance(window, tuple):
        window = tuple(_wrap_set(w) for w in window)
    rv = solve(_wrap_set(expr), _wrap_set(gen), _wrap_set(domain), window)
    return None if rv is None else _unwrap(rv)


def solve(expr, gen, domain, window=None, continuous=False):
    """Return the solution set, ``None`` when ``expr`` is not periodic in
    ``gen`` (the caller keeps its own path), or raise
    ``NotImplementedError``."""
    from .. import S, Dummy, Interval, sympify
    from ..calculus.util import periodicity
    expr = sympify(expr)
    if not getattr(expr, 'is_Relational', False) or expr.rel_op == '==':
        return None
    domain = S.Reals if domain is None else domain
    _gen = gen
    if gen.is_extended_real is False:
        return None
    if gen.is_extended_real is None:
        gen = Dummy('gen', extended_real=True)
        expr = expr.xreplace({_gen: gen})
    e = expr.lhs - expr.rhs
    try:
        period = periodicity(e, gen)
    except (NotImplementedError, ValueError, TypeError):
        return None
    if period is None or period == S.Zero or period.is_positive is not True:
        return None
    if window is False:
        rv = _range_decides(expr, e, gen, S.Reals)
        if rv is None:
            raise _cannot(expr, gen)
    elif domain == S.Reals:
        inside = True
        if window is not None:
            lo, hi = window
            inside = (lo >= 0) == True and (hi <= period) == True
        if inside:
            rv = _one_period(expr, e, gen, period, continuous)
        elif lo.is_number and hi.is_number and lo.is_finite and hi.is_finite and (lo < hi) == True:
            # only the values on the integration interval matter
            rv = _range_decides(expr, e, gen, S.Reals)
            if rv is None:
                rv = _solve_on(expr, e, gen, Interval(lo, hi), continuous)
        else:
            raise _cannot(expr, gen)
    else:
        if domain.is_subset(S.Reals) is not True:
            return None
        rv = _range_decides(expr, e, gen, domain)
        if rv is None:
            window_rv = _one_period(expr, e, gen, period, continuous).intersect(domain)
            if domain.is_subset(Interval(0, period, False, True)) is True:
                rv = window_rv
            elif domain.inf.is_finite and domain.sup.is_finite:
                rv = _solve_on(expr, e, gen, domain, continuous)
                if rv == window_rv:
                    rv = window_rv
            else:
                raise _cannot(expr, gen)
    return rv.xreplace({gen: _gen}) if gen != _gen else rv
