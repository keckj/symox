"""Limits of sequences (port of ``sympy.series.limitseq``)."""
import operator

from .. import difference_delta


def dominant(expr, n):
    """Finds the dominant term in a sum, that is a term that dominates
    every other term.

    If ``limit(a/b, n, oo)`` is ``oo`` then ``a`` dominates ``b``; if it is
    0 then ``b`` dominates ``a``; otherwise they are comparable. Returns
    ``None`` when there is no unique dominant term.

    >>> from symox import Sum
    >>> from symox.series.limitseq import dominant
    >>> from symox.abc import n, k
    >>> dominant(5*n**3 + 4*n**2 + n + 1, n)
    5*n**3
    >>> dominant(2**n + Sum(k, (k, 0, n)), n)
    2**n
    """
    from .. import Add, S
    global _dominant_depth
    if _dominant_depth >= _MAX_DOMINANT_DEPTH:
        # The reference recurses until gruntz decides each ratio; the native
        # engine leaves more ratios undecided, so the recursion is bounded.
        return None
    _dominant_depth += 1
    try:
        return _dominant(expr, n)
    finally:
        _dominant_depth -= 1


_dominant_depth = 0
_MAX_DOMINANT_DEPTH = 4


def _dominant(expr, n):
    from .. import Add, S
    terms = Add.make_args(expr.expand(func=True))
    term0 = terms[-1]
    comp = [term0]  # comparable terms
    for t in terms[:-1]:
        r = term0/t
        e = r.gammasimp()
        if e == r:
            e = r.factor()
        l = limit_seq(e, n)
        if l is None:
            return None
        elif l.is_zero:
            term0 = t
            comp = [term0]
        elif l not in [S.Infinity, S.NegativeInfinity]:
            comp.append(t)
    if len(comp) > 1:
        return None
    return term0


_HYPERGEOMETRIC = frozenset({'factorial', 'factorial2', 'subfactorial', 'gamma',
                             'binomial', 'RisingFactorial', 'FallingFactorial',
                             'loggamma', 'catalan'})


def _has_hypergeometric(expr, n):
    from ..solvers._differential import nodes
    return any(node.has(n) and (node.is_UndefinedFunction or
               getattr(node.func, '__name__', '') in _HYPERGEOMETRIC)
               for node in nodes(expr))


def _native_limit(expr, n):
    """The native limit at ``oo``, or ``None`` when it is not determined.

    The native gruntz engine is not trusted on Gamma-type terms (it misses
    their asymptotic expansion), so those are left to the ratio test."""
    from .. import Limit, S, PoleError
    if _has_hypergeometric(expr, n):
        return None
    try:
        result = Limit(expr, n, S.Infinity).doit(deep=False)
    except (NotImplementedError, PoleError, ValueError):
        return None
    if isinstance(result, Limit) or result in (S.NaN, S.ComplexInfinity):
        return None
    return result


def _ratio_limit(expr, n):
    """Exact ratio test: ``|a(n+1)/a(n)| -> r < 1`` proves ``a(n) -> 0``;
    ``r > 1`` with positive ratios and a positive (negative) first term
    proves ``a(n) -> oo`` (``-oo``)."""
    from .. import S, Abs, combsimp, simplify, expand_func
    try:
        step = combsimp(expand_func(expr.subs(n, n + 1)/expr))
        step = simplify(step)
        ratio = _native_limit(Abs(step), n)
        if ratio is None:
            ratio = _native_limit(simplify(Abs(step)), n)
    except (NotImplementedError, ValueError, ZeroDivisionError):
        return None
    if ratio is None:
        return None
    if ratio.is_nonnegative is True and (1 - ratio).is_positive is True:
        return S.Zero
    if (ratio is S.Infinity or (ratio - 1).is_positive is True) and step.is_positive is True:
        initial = expr.subs(n, 1).doit()
        if initial.is_finite is True:
            if initial.is_positive is True:
                return S.Infinity
            if initial.is_negative is True:
                return S.NegativeInfinity
    return None


def _limit_inf(expr, n):
    """``Limit(expr, n, oo).doit(deep=False)``, with ``None`` where the
    reference raises or where the native engine does not decide."""
    from .. import Sum, Product, Integral, combsimp, simplify
    # The reference gruntz raises NotImplementedError on the unknown mrv set
    # of an unevaluated Sum, Product or Integral.
    if expr.has(Product) and not expr.has(Sum, Integral):
        # A closed form of the product keeps the earlier symox capability
        # (the reference does not terminate on these).
        evaluated = expr.doit()
        if evaluated.has(Product):
            return None
        expr = evaluated
    if expr.has(Sum, Product, Integral):
        return None
    result = _native_limit(expr, n)
    if result is not None:
        return result
    try:
        simpler = simplify(combsimp(expr))
    except (NotImplementedError, ValueError):
        return None
    if simpler != expr:
        result = _native_limit(simpler, n)
        if result is not None:
            return result
    # The ratio test serves the terms the native engine does not handle:
    # Gamma-type factors and powers of non-real numbers.
    from .. import Pow
    if _has_hypergeometric(simpler, n) or any(
            p.base.is_number and p.base.is_extended_real is False and p.exp.has(n)
            for p in simpler.atoms(Pow)):
        return _ratio_limit(simpler, n)
    return None


def _limit_seq(expr, n, trials):
    from .. import Sum
    for i in range(trials):
        if not expr.has(Sum):
            result = _limit_inf(expr, n)
            if result is not None:
                return result

        num, den = expr.as_numer_denom()
        if not den.has(n) or not num.has(n):
            result = _limit_inf(expr.doit(), n)
            if result is not None:
                return result
            return None

        num, den = (difference_delta(t.expand(), n) for t in [num, den])
        expr = (num / den).gammasimp()

        if not expr.has(Sum):
            result = _limit_inf(expr, n)
            if result is not None:
                return result

        num, den = expr.as_numer_denom()

        num = dominant(num, n)
        if num is None:
            return None

        den = dominant(den, n)
        if den is None:
            return None

        expr = (num / den).gammasimp()


def _rewrite_golden_ratio(expr):
    """``expr.rewrite(fibonacci, GoldenRatio)``, as in the reference
    ``fibonacci._eval_rewrite_as_GoldenRatio``."""
    from .. import S, fibonacci, lucas, Function
    mapping = {}
    for f in expr.atoms(Function):
        if isinstance(f, fibonacci) and len(f.args) == 1:
            m = f.args[0]
            mapping[f] = (S.GoldenRatio**m - 1/(-S.GoldenRatio)**m)/(2*S.GoldenRatio - 1)
        elif isinstance(f, lucas) and len(f.args) == 1:
            # The reference does not rewrite lucas, and its gruntz then gets
            # lucas(n)/fibonacci(n) -> 0 wrong; the Binet form gives sqrt(5).
            m = f.args[0]
            mapping[f] = S.GoldenRatio**m + 1/(-S.GoldenRatio)**m
    if not mapping:
        return expr
    # Inner arguments first, as rewrite works bottom up.
    return expr.xreplace({f: v.xreplace(mapping) for f, v in mapping.items()})


def _parity_signs(expr):
    """``(-b)**e`` as ``b**e`` or ``-b**e`` for a number ``b > 0`` and an
    exponent of known parity, which the reference Pow constructor does for
    a symbol exponent of an odd or even index."""
    from .. import Pow, S
    mapping = {}
    for p in expr.atoms(Pow):
        b, e = p.args
        if b.is_number and b.is_negative is True and b is not S.NegativeOne:
            if e.is_even is True:
                mapping[p] = Pow(-b, e)
            elif e.is_odd is True:
                mapping[p] = -Pow(-b, e)
    return expr.xreplace(mapping) if mapping else expr



def limit_seq(expr, n=None, trials=5):
    """Finds the limit of a sequence as index ``n`` tends to infinity.

    Parameters
    ==========

    expr : Expr
        Expression for the ``n-th`` term of the sequence
    n : Symbol, optional
        The index of the sequence, an integer that tends to positive
        infinity. If None, inferred from the expression unless it has
        multiple symbols.
    trials: int, optional
        The algorithm is highly recursive. ``trials`` is a safeguard from
        infinite recursion in case the limit is not easily computed by the
        algorithm. Try increasing ``trials`` if the algorithm returns ``None``.

    Examples
    ========

    >>> from symox import limit_seq, Sum, binomial
    >>> from symox.abc import n, k, m
    >>> limit_seq((5*n**3 + 3*n**2 + 4) / (3*n**3 + 4*n - 5), n)
    5/3
    >>> limit_seq(binomial(2*n, n) / Sum(binomial(2*k, k), (k, 1, n)), n)
    3/4
    >>> limit_seq(Sum(k**2 * Sum(2**m/m, (m, 1, k)), (k, 1, n)) / (2**n*n), n)
    4

    References
    ==========

    .. [1] Computing Limits of Sequences - Manuel Kauers
    """
    from .. import (Sum, S, Dummy, Pow, Abs, Min, Max, cos, sin, Add,
                    factorial, subfactorial, gamma)
    from .. import sympify
    from ..calculus.accumulationbounds import AccumulationBounds

    expr = sympify(expr)
    if n is None:
        free = expr.free_symbols
        if len(free) == 1:
            n = free.pop()
        elif not free:
            return expr
        else:
            raise ValueError("Expression has more than one variable. "
                             "Please specify a variable.")
    elif n not in expr.free_symbols:
        return expr

    trials = operator.index(trials)
    expr = _rewrite_golden_ratio(expr)
    expr = expr.rewrite(factorial, subfactorial, gamma)

    n_ = Dummy("n", integer=True, positive=True)
    n1 = Dummy("n", odd=True, positive=True)
    n2 = Dummy("n", even=True, positive=True)

    # If there is a negative term raised to a power involving n, or a
    # trigonometric function, then consider even and odd n separately.
    powers = (p.as_base_exp() for p in expr.atoms(Pow))
    if (any(b.is_negative and e.has(n) for b, e in powers) or
            expr.has(cos, sin)):
        L1 = _limit_seq(_parity_signs(expr.xreplace({n: n1})), n1, trials)
        if L1 is not None:
            L2 = _limit_seq(_parity_signs(expr.xreplace({n: n2})), n2, trials)
            if L2 is None:
                return None
            if L1 != L2:
                if L1.is_comparable and L2.is_comparable:
                    return AccumulationBounds(Min(L1, L2), Max(L1, L2))
                else:
                    return None
    else:
        L1 = _limit_seq(expr.xreplace({n: n_}), n_, trials)
    if L1 is not None:
        return L1
    else:
        if expr.is_Add:
            limits = [limit_seq(term, n, trials) for term in expr.args]
            if any(result is None for result in limits):
                return None
            else:
                return Add(*limits)
        # Maybe the absolute value is easier to deal with (though not if
        # it has a Sum). If it tends to 0, the limit is 0.
        elif not expr.has(Sum):
            lim = _limit_seq(Abs(expr.xreplace({n: n_})), n_, trials)
            if lim is not None and lim.is_zero:
                return S.Zero


__all__ = ['difference_delta', 'dominant', 'limit_seq']
