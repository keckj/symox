"""Structured partial fraction decomposition (``sympy.polys.partfrac``).

``apart_list`` is Bronstein's full decomposition algorithm as implemented in
SymPy: each square-free factor ``d`` of the denominator yields terms
``RootSum(D(w), Lambda(a, b(a))/(x - a)**k)`` stored as
``(D(w), Lambda(a, b), Lambda(a, x - a), k)``. ``assemble_partfrac_list``
expands the roots of linear ``D`` (and explicit root lists); higher degree
``D`` would need ``RootSum``, which is not implemented.
"""
from .._core import Poly, _apart_full_decomposition
from .polytools import allowed_flags, parallel_poly_from_expr

__all__ = ['apart_list', 'assemble_partfrac_list']


def apart_list(f, x=None, dummies=None, **options):
    """Partial fraction decomposition as ``(common, polypart, terms)``."""
    from .. import Dummy, sympify
    allowed_flags(options, [])
    f = sympify(f)
    if f.is_Atom:
        return f
    P, Q = f.as_numer_denom()
    options = dict(options)
    options.setdefault('extension', True)
    gens = () if x is None else (x,)
    (P, Q), opt = parallel_poly_from_expr((P, Q), *gens, **options)
    if P.is_multivariate:
        raise NotImplementedError("multivariate partial fraction decomposition")
    common, P, Q = P.cancel(Q)
    poly, P = P.div(Q, auto=True)
    P, Q = P.rat_clear_denoms(Q)
    polypart = poly
    if dummies is None:
        def _dummies(name):
            d = Dummy(name)
            while True:
                yield d
        dummies = _dummies("w")
    rationalpart = apart_list_full_decomposition(P, Q, dummies)
    return (common, polypart, rationalpart)


def apart_list_full_decomposition(P, Q, dummygen):
    """Bronstein's full partial fraction decomposition of ``P/Q``.

    Over ``ZZ``/``QQ`` the coefficients are Taylor coefficients computed in
    ``QQ[a]/(d(a))`` natively (``symox::poly::partfrac``); other domains use
    SymPy's derivative-substitution formulation below.
    """
    from .. import Dummy, Function, Lambda, cancel
    factors = Q.sqf_list_include(all=True)
    native = _apart_full_decomposition(P, [(d, int(n)) for d, n in factors])
    if native is not None:
        x, a = P.gen, Dummy('a')
        return [(D.subs(x, next(dummygen)), Lambda(a, bb.as_expr().subs(x, a)),
                 Lambda(a, x - a), k) for D, bb, k in native]
    f, x, U = P.as_expr()/Q.as_expr(), P.gen, []
    u = Function('u')(x)
    a = Dummy('a')
    partial = []
    for d, n in factors:
        b = d.as_expr()
        U += [u.diff(x, n - 1)]
        h = cancel(f*b**n)/u**n
        H, subs = [h], []
        for j in range(1, n):
            H += [H[-1].diff(x)/j]
        for j in range(1, n + 1):
            subs += [(U[j - 1], b.diff(x, j)/j)]
        for j in range(0, n):
            Pj, Qj = cancel(H[j]).as_numer_denom()
            for i in range(0, j + 1):
                Pj = Pj.subs(*subs[j - i])
            Qj = Qj.subs(*subs[0])
            Pj = Poly(Pj, x)
            Qj = Poly(Qj, x)
            G = Pj.gcd(d)
            D = d.quo(G)
            B, g = Qj.half_gcdex(D)
            bb = (Pj*B.quo(g)).rem(D)
            Dw = D.subs(x, next(dummygen))
            numer = Lambda(a, bb.as_expr().subs(x, a))
            denom = Lambda(a, (x - a))
            exponent = n - j
            partial.append((Dw, numer, denom, exponent))
    return partial


def _root_sum(r, var, expr):
    """``RootSum(r, Lambda(var, expr), auto=False)`` when every irreducible
    factor of ``r`` is linear; otherwise ``NotImplementedError``."""
    from .. import Add, S
    deg = r.degree()
    if not expr.has(var):
        return deg*expr
    add_const = S.Zero
    if expr.is_Add:
        add_const, expr = expr.as_independent(var)
    mul_const = S.One
    if expr.is_Mul:
        mul_const, expr = expr.as_independent(var)
    terms = []
    for factor, k in r.factor_list()[1]:
        if factor.degree() != 1:
            raise NotImplementedError(
                "RootSum over the roots of %s is not implemented" % factor.as_expr())
        root = -factor.TC()/factor.LC()
        terms.append(k*expr.subs(var, root))
    return mul_const*Add(*terms) + deg*add_const


def assemble_partfrac_list(partial_list):
    """Rebuild an expression from the output of :func:`apart_list`."""
    common = partial_list[0]
    polypart = partial_list[1]
    pfd = polypart.as_expr()
    for r, nf, df, ex in partial_list[2]:
        if isinstance(r, Poly):
            an, nu = nf.variables, nf.expr
            ad, de = df.variables, df.expr
            de = de.subs(ad[0], an[0])
            pfd += _root_sum(r, an[0], nu/de**ex)
        else:
            for root in r:
                pfd += nf(root)/df(root)**ex
    return common*pfd
