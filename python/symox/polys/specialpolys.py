"""Special polynomials (``sympy.polys.specialpolys``).

Cyclotomic coefficients come from the native integer routine; the other
constructions follow SymPy's definitions with exact symox arithmetic.
"""
import random as _random

from .._core import Poly

__all__ = [
    'swinnerton_dyer_poly', 'cyclotomic_poly', 'symmetric_poly',
    'random_poly', 'interpolating_poly',
]


def _pure(poly, gen):
    """A ``PurePoly`` in ``gen`` with the coefficients of the univariate ``poly``."""
    from .. import PurePoly
    try:
        domain = poly.domain
    except NotImplementedError:
        domain = poly.domain_descriptor
    return PurePoly._from_coeffs(poly.all_coeffs(), gen, domain, poly.as_expr())


def _named_result(coeffs, x, domain, polys):
    """Poly (``x`` given) or PurePoly in a Dummy from descending ``coeffs``."""
    from .. import Add, Dummy, sympify
    gen = Dummy('x') if x is None else sympify(x)
    slot = gen if gen.is_Symbol else Dummy('x')
    n = len(coeffs) - 1
    expr = Add(*[c*slot**(n - i) for i, c in enumerate(coeffs)])
    poly = Poly(expr, slot, domain=domain) if domain is not None else Poly(expr, slot)
    if slot is not gen:
        poly = poly._with_gens(gen)
    if x is None:
        poly = _pure(poly, gen)
    return poly if polys else poly.as_expr()


def swinnerton_dyer_poly(n, x=None, polys=False):
    """Swinnerton-Dyer polynomial of order ``n``: the minimal polynomial of
    ``sqrt(2) + sqrt(3) + ... + sqrt(p_n)``."""
    from .. import Dummy, nextprime, sympify
    if n <= 0:
        raise ValueError("Cannot generate Swinnerton-Dyer polynomial of order %s" % n)
    pure = x is None
    x = Dummy('x') if pure else sympify(x)
    if n > 3:
        p, primes = 2, [2]
        for _ in range(2, n + 1):
            p = nextprime(p)
            primes.append(p)
        ex = _sqrt_sum_norm(primes, x)
        if polys:
            return Poly(ex, x, domain='QQ')
        return ex
    if n == 1:
        ex = x**2 - 2
    elif n == 2:
        ex = x**4 - 10*x**2 + 1
    else:
        ex = x**8 - 40*x**6 + 352*x**4 - 960*x**2 + 576
    return _pure(Poly(ex, x), x) if polys else ex


def _sqrt_sum_norm(primes, x):
    """``prod(x - sum(+-sqrt(p)))``: eliminate one square root at a time.

    With ``f(x + t) = A(x, t**2) + t*B(x, t**2)`` the product of the two
    conjugates ``f(x + sqrt(p))*f(x - sqrt(p))`` is ``A**2 - p*B**2``.
    """
    from .. import Dummy, expand
    t = Dummy('t')
    f = Poly(x**2 - primes[0], x)
    for p in primes[1:]:
        shifted = Poly(expand(f.as_expr().subs(x, x + t)), t)
        even, odd = 0, 0
        for (k,), c in shifted.terms():
            if k % 2:
                odd += c*p**(k//2)
            else:
                even += c*p**(k//2)
        f = Poly(expand(even**2 - p*odd**2), x)
    return f.as_expr()


def cyclotomic_poly(n, x=None, polys=False):
    """The ``n``-th cyclotomic polynomial."""
    from .. import Integer, _core
    if n <= 0:
        raise ValueError("Cannot generate cyclotomic polynomial of order %s" % n)
    coeffs = [Integer(c) for c in reversed(_core._cyclotomic_coefficients(int(n)))]
    return _named_result(coeffs, x, 'ZZ', polys)


def _analyze_gens(gens):
    if len(gens) == 1 and isinstance(gens[0], (list, tuple)):
        return tuple(gens[0])
    return tuple(gens)


def symmetric_poly(n, *gens, polys=False):
    """Elementary symmetric polynomial of degree ``n`` in ``gens``."""
    from itertools import combinations
    from .. import Add, Mul, S, sympify
    gens = tuple(sympify(g) for g in _analyze_gens(gens))
    if n < 0 or n > len(gens) or not gens:
        raise ValueError("Cannot generate symmetric polynomial of order %s for %s"
                         % (n, '(%s)' % ', '.join(map(str, gens)) if len(gens) != 1 else '(%s,)' % gens[0]))
    if not n:
        poly = S.One
    else:
        poly = Add(*[Mul(*s) for s in combinations(gens, int(n))])
    return Poly(poly, *gens) if polys else poly


def random_poly(x, n, inf, sup, domain='ZZ', polys=False):
    """Univariate polynomial of degree ``n`` with random coefficients in
    ``[inf, sup]`` and a nonzero leading coefficient."""
    from .. import Integer
    coeffs = [_random.randint(inf, sup) for _ in range(n + 1)]
    while not coeffs[0]:
        coeffs[0] = _random.randint(inf, sup)
    coeffs = [Integer(c) for c in coeffs]
    from .. import Add, sympify
    x = sympify(x)
    expr = Add(*[c*x**(n - i) for i, c in enumerate(coeffs)])
    poly = Poly(expr, x, domain=domain)
    return poly if polys else poly.as_expr()


def interpolating_poly(n, x, X='x', Y='y'):
    """Lagrange interpolating polynomial for ``n`` points ``(X[i], Y[i])``."""
    from .. import Add, Mul, Tuple, symbols, sympify
    x = sympify(x)
    ok = getattr(x, 'free_symbols', None)
    if isinstance(X, str):
        X = symbols("%s:%s" % (X, n))
    elif ok and ok & Tuple(*X).free_symbols:
        ok = False
    if isinstance(Y, str):
        Y = symbols("%s:%s" % (Y, n))
    elif ok and ok & Tuple(*Y).free_symbols:
        ok = False
    if not ok:
        raise ValueError("\nExpecting symbol for x that does not appear in X or Y. Use\n"
                         "`interpolate(list(zip(X, Y)), x)` instead.")
    X = [sympify(v) for v in X]
    Y = [sympify(v) for v in Y]
    coeffs = []
    numert = Mul(*[x - X[i] for i in range(n)])
    for i in range(n):
        numer = numert/(x - X[i])
        denom = Mul(*[(X[i] - X[j]) for j in range(n) if i != j])
        coeffs.append(numer/denom)
    return Add(*[coeff*y for coeff, y in zip(coeffs, Y)])
