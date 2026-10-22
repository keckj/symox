"""Dense univariate polynomial helpers used by the DomainMatrix algorithms.

These are the reference ``sympy.polys.dense*`` functions that
``sympy.polys.matrices`` calls, ported over native domain elements. A dense
polynomial is a list of coefficients, highest degree first, with no leading
zeros.
"""
from collections import OrderedDict


def dup_strip(f):
    if not f or f[0]:
        return f
    i = 0
    for c in f:
        if c:
            break
        i += 1
    return f[i:]


def dup_convert(f, K0, K1):
    if K0 is not None and K0 == K1:
        return f
    return dup_strip([K1.convert(c, K0) if K0 is not None else K1.convert(c) for c in f])


def dup_add(f, g, K):
    if not f:
        return g
    if not g:
        return f
    df, dg = len(f) - 1, len(g) - 1
    if df == dg:
        return dup_strip([a + b for a, b in zip(f, g)])
    k = abs(df - dg)
    if df > dg:
        h, f = f[:k], f[k:]
    else:
        h, g = g[:k], g[k:]
    return h + [a + b for a, b in zip(f, g)]


def dup_mul(f, g, K):
    if not (f and g):
        return []
    df, dg = len(f) - 1, len(g) - 1
    h = []
    for i in range(0, df + dg + 1):
        coeff = K.zero
        for j in range(max(0, i - dg), min(df, i) + 1):
            coeff += f[j]*g[i - j]
        h.append(coeff)
    return dup_strip(h)


def dup_mul_ground(f, c, K):
    if not c or not f:
        return []
    return [cf * c for cf in f]


def dup_quo_ground(f, c, K):
    if not c:
        raise ZeroDivisionError('polynomial division')
    if not f:
        return f
    if K.is_Field:
        return [K.quo(cf, c) for cf in f]
    return [cf // c for cf in f]


def dup_content(f, K):
    if not f:
        return K.zero
    cont = K.zero
    if K.is_QQ:
        for c in f:
            cont = K.gcd(cont, c)
    else:
        for c in f:
            cont = K.gcd(cont, c)
            if K.is_one(cont):
                break
    return cont


def dup_primitive(f, K):
    if not f:
        return K.zero, f
    cont = dup_content(f, K)
    if K.is_one(cont):
        return cont, f
    return cont, dup_quo_ground(f, cont, K)


def dup_clear_denoms(f, K0, K1=None, convert=False):
    if K1 is None:
        K1 = K0.get_ring() if K0.has_assoc_Ring else K0
    common = K1.one
    for c in f:
        common = K1.lcm(common, K0.denom(c))
    if K1.is_one(common):
        if not convert:
            return common, f
        return common, dup_convert(f, K0, K1)
    # The reference uses quo rather than exquo so that inexact domains
    # discard the remainder.
    f = [K0.numer(c)*K1.quo(common, K0.denom(c)) for c in f]
    if not convert:
        return common, dup_convert(f, K1, K0)
    return common, f


def dup_transform(f, p, q, K):
    if not f:
        return []
    n = len(f) - 1
    h, Q = [f[0]], [[K.one]]
    for i in range(0, n):
        Q.append(dup_mul(Q[-1], q, K))
    for c, q in zip(f[1:], Q[1:]):
        h = dup_mul(h, p, K)
        q = dup_mul_ground(q, c, K)
        h = dup_add(h, q, K)
    return h


def dup_factor_list(f, K):
    """The reference ``dup_factor_list``: the content and the sorted
    irreducible factors with multiplicities, as dense lists over ``K``."""
    from ...core.symbol import Dummy
    from ..polytools import Poly
    if K.is_EX or K.is_EXRAW:
        # The reference dup_factor_list has no expression-domain case.
        from ..polyerrors import DomainError
        raise DomainError('factorization not supported over %s' % K)
    if not f:
        return K.zero, []
    x = Dummy('x')
    coeff, factors = Poly.from_list(list(f), x, domain=K).factor_list()
    return K.convert(coeff), [([K.convert(c) for c in g.rep.to_list()], k)
                              for g, k in factors]


def _order_key(factor):
    """The reference ordering of factor coefficients.

    Native symox elements share one Python type, so the orderings of the
    reference element types are reproduced here: finite field elements by
    their integer value, algebraic numbers by their coefficient lists and
    polynomials by ``(len(p), p.terms())``.
    """
    if isinstance(factor, (list, tuple)):
        return [_order_key(f) for f in factor]
    K = getattr(factor, 'parent', None)
    if K is None:
        return factor
    K = K()
    if K.is_ZZ or K.is_QQ:
        return factor
    if K.is_FF:
        return int(factor)
    if K.is_ZZ_I or K.is_QQ_I:
        return [_order_key(factor.y), _order_key(factor.x)]
    if K.is_EX or K.is_EXRAW:
        # The reference compares EXRAW expressions with ``<``, which raises
        # TypeError for symbolic factors; sort_key gives a total order.
        return K.to_sympy(factor).sort_key()
    if K.is_AlgebraicField:
        return _order_key(factor.to_list())
    if K.is_PolynomialRing:
        from ..polytools import Poly
        terms = Poly(K.to_sympy(factor), *K.symbols, domain=K.domain).terms()
        return (len(terms), [(m, _order_key(K.domain.from_sympy(c))) for m, c in terms])
    return factor


def _sort_factors(factors, **args):
    def order_if_multiple_key(factor):
        (f, n) = factor
        return (len(f), n, _order_key(f))

    def order_no_multiple_key(f):
        return (len(f), _order_key(f))

    if args.get('multiple', True):
        return sorted(factors, key=order_if_multiple_key)
    return sorted(factors, key=order_no_multiple_key)


def _strongly_connected_components(V, Gmap):
    """Non-recursive Tarjan algorithm, in the reference vertex order."""
    lowlink = {}
    indices = {}
    stack = OrderedDict()
    callstack = []
    components = []
    nomore = object()

    def start(v):
        index = len(stack)
        indices[v] = lowlink[v] = index
        stack[v] = None
        callstack.append((v, iter(Gmap[v])))

    def finish(v1):
        if lowlink[v1] == indices[v1]:
            component = [stack.popitem()[0]]
            while component[-1] != v1:
                component.append(stack.popitem()[0])
            components.append(component[::-1])
        v2, _ = callstack.pop()
        if callstack:
            v1, _ = callstack[-1]
            lowlink[v1] = min(lowlink[v1], lowlink[v2])

    for v in V:
        if v in indices:
            continue
        start(v)
        while callstack:
            v1, it1 = callstack[-1]
            v2 = next(it1, nomore)
            if v2 is nomore:
                finish(v1)
            elif v2 not in indices:
                start(v2)
            elif v2 in stack:
                lowlink[v1] = min(lowlink[v1], indices[v2])

    return components


def doctest_depends_on(**dependencies):
    """The reference doctest marker; it leaves the object unchanged."""
    def decorator(obj):
        obj._doctest_depends_on = dependencies
        return obj
    return decorator


def _native_elements(rows, domain):
    """Native symox domain elements of ``rows``.

    The reference stores the given objects as they are: its ground types
    accept Python ints and SymPy numbers in place of ``ZZ`` and ``QQ``
    elements. Native symox elements are a distinct type, so entries that are
    not elements of ``domain`` yet are converted here, exactly.
    """
    of_type, convert = domain.of_type, domain.convert
    try:
        if isinstance(rows, list):
            if all(of_type(e) for row in rows for e in row):
                return rows
            return [[e if of_type(e) else convert(e) for e in row] for row in rows]
        if isinstance(rows, dict):
            if all(of_type(e) for row in rows.values() for e in row.values()):
                return rows
            return {i: {j: e if of_type(e) else convert(e) for j, e in row.items()}
                    for i, row in rows.items()}
    except (TypeError, AttributeError):
        return rows
    return rows
