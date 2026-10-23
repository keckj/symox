"""Structured partial fraction decomposition (``sympy.polys.partfrac``).

``apart_list`` is Bronstein's full decomposition algorithm as implemented in
SymPy: each square-free factor ``d`` of the denominator yields terms
``RootSum(D(w), Lambda(a, b(a))/(x - a)**k)`` stored as
``(D(w), Lambda(a, b), Lambda(a, x - a), k)``. ``assemble_partfrac_list``
sums over the roots of ``D`` with ``RootSum(D, ..., auto=False)``, which
expands linear factors and keeps the others unevaluated.
"""
from .._core import Poly, _apart_full_decomposition
from .polytools import allowed_flags, parallel_poly_from_expr

__all__ = ['apart', 'apart_list', 'assemble_partfrac_list']


def __getattr__(name):
    # ``apart`` is the package-level function, bound lazily because it is
    # rebound (noncommutative handling) after this module is imported.
    if name == 'apart':
        import symox
        return symox.apart
    raise AttributeError(name)


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


def assemble_partfrac_list(partial_list):
    """Rebuild an expression from the output of :func:`apart_list`."""
    from .. import Lambda
    from .rootoftools import RootSum
    common = partial_list[0]
    polypart = partial_list[1]
    pfd = polypart.as_expr()
    for r, nf, df, ex in partial_list[2]:
        if isinstance(r, Poly):
            an, nu = nf.variables, nf.expr
            ad, de = df.variables, df.expr
            de = de.subs(ad[0], an[0])
            func = Lambda(tuple(an), nu/de**ex)
            pfd += RootSum(r, func, auto=False, quadratic=False)
        else:
            for root in r:
                pfd += nf(root)/df(root)**ex
    return common*pfd


def _apart_full(f, x=None, **options):
    """The ``full=True`` branch of the reference ``apart`` for commutative
    univariate input (Bronstein's algorithm through :func:`apart_list` and
    :func:`assemble_partfrac_list`). ``None`` defers to the native ``apart``.
    """
    from .. import Add, sympify, factor
    from .polyerrors import PolynomialError
    from .rootoftools import RootSum
    f = sympify(f)
    if f.is_Atom or not f.is_commutative:
        return None
    P, Q = f.as_numer_denom()
    options = dict(options)
    options.setdefault('extension', True)
    gens = () if x is None else (x,)
    try:
        (P, Q), opt = parallel_poly_from_expr((P, Q), *gens, **options)
    except PolynomialError:
        return None
    if P.is_multivariate:
        return None
    common, P, Q = P.cancel(Q)
    poly, P = P.div(Q, auto=True)
    P, Q = P.rat_clear_denoms(Q)
    if Q.degree() <= 1:
        partial = P.as_expr()/Q.as_expr()
    else:
        partial = assemble_partfrac_list(apart_list(P.as_expr()/Q.as_expr(), P.gens[0]))
    terms = []
    for term in Add.make_args(partial):
        terms.append(term if term.has(RootSum) else factor(term))
    return common*(poly.as_expr() + Add(*terms))


def _apart_inexact(f, x=None, **options):
    """SymPy's ``apart`` for an expression with Float coefficients.

    Ported from ``sympy/polys/partfrac.py`` (BSD-3-Clause, see
    LICENSE-SYMPY): the decomposition stays over the inexact domain (RR or
    CC) by undetermined coefficients, so floats are not turned into exact
    binary rationals. Returns None when the native path should decide.
    """
    from .. import Add, Dummy, Float, S, factor, numbered_symbols, solve, sympify
    from .. import Expr
    if options:
        return None
    f = sympify(f)
    if not isinstance(f, Expr) or f.is_Atom or not f.is_commutative or not f.has(Float):
        return None
    P, Q = f.as_numer_denom()
    gens = () if x is None else (x,)
    try:
        (P, Q), opt = parallel_poly_from_expr((P, Q), *gens, extension=True)
    except Exception:
        return None
    if P.is_multivariate or P.get_domain().is_Exact:
        return None
    common, P, Q = P.cancel(Q)
    poly, P = P.div(Q, auto=True)
    P, Q = P.rat_clear_denoms(Q)
    if Q.degree() <= 1:
        partial = P.as_expr()/Q.as_expr()
    else:
        X = numbered_symbols(cls=Dummy)
        partial, symbols = [], []
        _, factors = Q.factor_list()
        for g, k in factors:
            n, q = g.degree(), Q
            for i in range(1, k + 1):
                coeffs = [next(X) for _ in range(n)]
                q = q.quo(g)
                partial.append((coeffs, q, g, i))
                symbols.extend(coeffs)
        dom = Q.get_domain().inject(*symbols)
        F = Poly(0, Q.gen, domain=dom)
        for i, (coeffs, q, g, k) in enumerate(partial):
            h = Poly(coeffs, Q.gen, domain=dom)
            partial[i] = (h, g, k)
            F += h*q.set_domain(dom)
        system = [coeff - P.nth(k) for (k,), coeff in F.terms()]
        solution = solve(system, symbols)
        result = S.Zero
        for h, g, k in partial:
            result += h.as_expr().subs(solution)/g.as_expr()**k
        partial = result
    terms = S.Zero
    for term in Add.make_args(partial):
        terms += factor(term)
    return common*(poly.as_expr() + terms)
