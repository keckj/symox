"""Hypergeometric summation (``eval_sum_hyper``), ported from SymPy 1.14
``sympy/concrete/summations.py``.

The native summation engine runs the polynomial, geometric and Gosper
stages of the reference ``eval_sum_symbolic``; a single-limit ``Sum`` it
leaves unevaluated is handed to the reference ``eval_sum_hyper`` stage here,
which writes the series as a ``hyper`` function and expands it with
``hyperexpand`` under the convergence condition of the series.
"""

from symox import S, oo


def radius_of_convergence(h):
    """``hyper.radius_of_convergence`` of the reference."""
    ap, bq = list(h.ap), list(h.bq)
    if any(a.is_integer and (a <= 0) == True for a in ap + bq):
        aints = [a for a in ap if a.is_Integer and (a <= 0) == True]
        bints = [a for a in bq if a.is_Integer and (a <= 0) == True]
        if len(aints) < len(bints):
            return S.Zero
        popped = False
        for b in bints:
            cancelled = False
            while aints:
                a = aints.pop()
                if a >= b:
                    cancelled = True
                    break
                popped = True
            if not cancelled:
                return S.Zero
        if aints or popped:
            return oo
    if len(ap) == len(bq) + 1:
        return S.One
    elif len(ap) <= len(bq):
        return oo
    else:
        return S.Zero


def convergence_statement(ap, bq, z):
    """``hyper.convergence_statement`` of the reference, from the parameters
    before the native ``hyper`` constructor cancels equal ones."""
    from symox import And, Or, Ne, re, Abs, Add
    R = _radius(ap, bq)
    if R == 0:
        return S.false
    if R == oo:
        return S.true
    e = Add(*[S(a) for a in ap]) - Add(*[S(b) for b in bq])
    c1 = And(re(e) < 0, Abs(z) <= 1)
    c2 = And(0 <= re(e), re(e) < 1, Abs(z) <= 1, Ne(z, 1))
    c3 = And(re(e) >= 1, Abs(z) < 1)
    return Or(c1, c2, c3)


class _Params:
    def __init__(self, ap, bq):
        self.ap, self.bq = ap, bq


def _radius(ap, bq):
    from symox import sympify
    return radius_of_convergence(_Params([sympify(a) for a in ap], [sympify(b) for b in bq]))


def _eval_sum_hyper(f, i, a):
    """Returns (res, cond). Sums from a to oo."""
    from symox import Dummy, Float, Poly, factor, hyper
    if a != 0:
        return _eval_sum_hyper(f.subs(i, i + a), i, 0)

    if f.subs(i, 0) == 0:
        from symox.simplify.simplify import simplify
        if simplify(f.subs(i, Dummy('i', integer=True, positive=True))) == 0:
            return S.Zero, True
        return _eval_sum_hyper(f.subs(i, i + 1), i, 0)

    from symox.simplify.simplify import hypersimp
    hs = hypersimp(f, i)
    if hs is None:
        return None

    if isinstance(hs, Float):
        from symox.simplify.simplify import nsimplify
        hs = nsimplify(hs)

    from symox.simplify.combsimp import combsimp
    from symox.simplify.hyperexpand import hyperexpand
    from symox.simplify.radsimp import fraction
    numer, denom = fraction(factor(hs))
    top, topl = numer.as_coeff_mul(i)
    bot, botl = denom.as_coeff_mul(i)
    ab = [top, bot]
    factors = [topl, botl]
    params = [[], []]
    for k in range(2):
        for fac in factors[k]:
            mul = 1
            if fac.is_Pow:
                mul = fac.exp
                fac = fac.base
                if not mul.is_Integer:
                    return None
            p = Poly(fac, i)
            if p.degree() != 1:
                return None
            m, n = p.all_coeffs()
            ab[k] *= m**mul
            params[k] += [n/m]*mul

    # Add "1" to numerator parameters, to account for implicit n! in
    # hypergeometric series.
    ap = params[0] + [1]
    bq = params[1]
    x = ab[0]/ab[1]
    h = hyper(ap, bq, x)
    f = combsimp(f)
    return f.subs(i, 0)*hyperexpand(h), convergence_statement(ap, bq, x)


def eval_sum_hyper(f, i_a_b):
    from symox import Sum, Piecewise, And, Dummy
    i, a, b = i_a_b

    if f.is_hypergeometric(i) is False:
        return

    if (b - a).is_Integer:
        return None

    old_sum = Sum(f, (i, a, b))

    if b != S.Infinity:
        if a is S.NegativeInfinity:
            res = _eval_sum_hyper(f.subs(i, -i), i, -b)
            if res is not None:
                return Piecewise(res, (old_sum, True))
        else:
            _illegal = (S.NaN, S.Infinity, S.NegativeInfinity, S.ComplexInfinity)
            n_illegal = lambda x: sum(x.count(_) for _ in _illegal)
            had = n_illegal(f)
            res1 = _eval_sum_hyper(f, i, a)
            if res1 is None or n_illegal(res1[0]) > had:
                return
            res2 = _eval_sum_hyper(f, i, b + 1)
            if res2 is None or n_illegal(res2[0]) > had:
                return
            (res1, cond1), (res2, cond2) = res1, res2
            cond = And(cond1, cond2)
            if cond == False:
                return None
            return Piecewise((res1 - res2, cond), (old_sum, True))

    if a is S.NegativeInfinity:
        res1 = _eval_sum_hyper(f.subs(i, -i), i, 1)
        res2 = _eval_sum_hyper(f, i, 0)
        if res1 is None or res2 is None:
            return None
        res1, cond1 = res1
        res2, cond2 = res2
        cond = And(cond1, cond2)
        if cond == False or cond.as_set() == S.EmptySet:
            return None
        return Piecewise((res1 + res2, cond), (old_sum, True))

    # Now b == oo, a != -oo
    res = _eval_sum_hyper(f, i, a)
    if res is not None:
        r, c = res
        if c == False:
            if r.is_number:
                f = f.subs(i, Dummy('i', integer=True, positive=True) + a)
                if f.is_positive or f.is_zero:
                    return S.Infinity
                elif f.is_negative:
                    return S.NegativeInfinity
            return None
        return Piecewise(res, (old_sum, True))


def _try_hyper(s):
    """The ``eval_sum_hyper`` value of a single-limit held ``Sum``, or
    None."""
    limits = s.limits
    if len(limits) != 1 or len(limits[0]) != 3:
        return None
    i, a, b = limits[0]
    f = s.function
    from symox import Sum, Integral, Product
    if f.has(Sum, Integral, Product) or not i.is_Symbol:
        return None
    if a.is_infinite and b.is_infinite and a == b:
        # Sum(f, (i, oo, oo)) has no agreed value (the reference gives nan
        # or oo depending on f): it stays held.
        return None
    if not (a.is_infinite or b.is_infinite):
        # Finite symbolic bounds: the native closed forms (Gosper) already
        # ran; the reference hyper stage answers only when both halves of
        # the split converge, which needs symbolic hyperexpand results the
        # native engine leaves alone.
        return None
    try:
        return eval_sum_hyper(f, (i, a, b))
    except (NotImplementedError, ValueError, TypeError, ZeroDivisionError,
            AttributeError, RecursionError):
        return None


def complete(result):
    """Apply the hypergeometric stage to the single-limit ``Sum`` nodes the
    native engine left in ``result``."""
    from symox import Sum
    if not getattr(result, 'has', None) or not result.has(Sum):
        return result
    cache = {}

    def visit(e):
        if e.is_Sum:
            if e in cache:
                return cache[e]
            value = _try_hyper(e)
            cache[e] = e if value is None else value
            return cache[e]
        if not e.args or e.is_Piecewise:
            # The fallback Sum of a native Piecewise closed form stays held.
            return e
        args = tuple(visit(a) for a in e.args)
        if args == e.args:
            return e
        return e.func(*args)
    try:
        return visit(result)
    except (TypeError, AttributeError):
        return result
