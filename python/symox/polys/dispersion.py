"""Dispersion of univariate polynomials (``sympy.polys.dispersion``).

The dispersion set of ``p`` and ``q`` is the set of non-negative integers
``a`` such that ``gcd(p(x), q(x + a))`` is not constant. Following SymPy, it
is computed by pairing the irreducible factors of ``p`` and ``q`` over the
rationals: only factors of equal degree and leading coefficient can differ
by a shift, and the shift is read off the subleading coefficients. Over
``ZZ``/``QQ`` the pairing runs natively (``symox::poly::dispersion``).
"""

from .._core import _dispersion_set

__all__ = ["dispersionset", "dispersion"]


def dispersionset(p, q=None, *gens, **args):
    """Set of integer shifts ``a >= 0`` for which ``p(x)`` and ``q(x + a)`` share a factor."""
    from .. import Poly, S
    same = q is None
    if same:
        q = p
    p = Poly(p, *gens, **args)
    q = Poly(q, *gens, **args)
    if not p.is_univariate or not q.is_univariate:
        raise ValueError("Polynomials need to be univariate")
    if not p.gen == q.gen:
        raise ValueError("Polynomials must have the same generator")
    gen = p.gen
    if p.degree() < 1 or q.degree() < 1:
        return {0}
    native = _dispersion_set(p, None if same else q)
    if native is not None:
        return set(native)
    fp = p.factor_list()
    fq = q.factor_list() if not same else fp
    J = set()
    for s, _ in fp[1]:
        for t, _ in fq[1]:
            m = s.degree()
            n = t.degree()
            if n != m:
                continue
            an = s.LC()
            bn = t.LC()
            if not (an - bn).is_zero:
                continue
            anm1 = s.coeff_monomial(gen**(m - 1))
            bnm1 = t.coeff_monomial(gen**(n - 1))
            alpha = (anm1 - bnm1)/S(n*bn)
            if not alpha.is_integer:
                continue
            if alpha < 0 or alpha in J:
                continue
            if n > 1 and not (s - t.shift(alpha)).is_zero:
                continue
            J.add(alpha)
    return J


def dispersion(p, q=None, *gens, **args):
    """Largest element of the dispersion set (``-oo`` when it is empty)."""
    from .. import S
    J = dispersionset(p, q, *gens, **args)
    if not J:
        return S.NegativeInfinity
    return max(J)
