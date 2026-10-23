"""Native minimal polynomials for the supported algebraic expressions."""
from ..._core import minimal_polynomial, minpoly

__all__ = ['minimal_polynomial', 'minpoly']


def _is_sum_surds(p):
    """SymPy's test for a sum of products of rationals and real square
    roots of rationals."""
    from ... import Add, Mul
    return all(f.is_Rational or f.is_Pow and
        f.base.is_Rational and (2*f.exp).is_Integer and f.is_extended_real
        for t in Add.make_args(p) for f in Mul.make_args(t))


def _separate_sq(p):
    """SymPy's helper: a polynomial without the square roots of ``p``
    (``p`` is a sum of ``c*sqrt(q)`` terms)."""
    from ... import S, Add, Mul
    from ...simplify.radsimp import _split_gcd
    from ...simplify.sqrtdenest import _mexpand
    from ...utilities.iterables import sift

    def is_sqrt(expr):
        return expr.is_Pow and expr.exp is S.Half
    a = []
    for y in p.args:
        if not y.is_Mul:
            if is_sqrt(y):
                a.append((S.One, y**2))
            elif y.is_Atom:
                a.append((y, S.One))
            elif y.is_Pow and y.exp.is_integer:
                a.append((y, S.One))
            else:
                raise NotImplementedError
        else:
            T, F = sift(y.args, is_sqrt, binary=True)
            a.append((Mul(*F), Mul(*T)**2))
    a.sort(key=lambda z: z[1])
    if a[-1][1] is S.One:
        return p
    surds = [z for y, z in a]
    for i in range(len(surds)):
        if surds[i] != 1:
            break
    g, b1, b2 = _split_gcd(*surds[i:])
    a1 = []
    a2 = []
    for y, z in a:
        if z in b1:
            a1.append(y*z**S.Half)
        else:
            a2.append(y*z**S.Half)
    p1 = Add(*a1)
    p2 = Add(*a2)
    return _mexpand(p1**2) - _mexpand(p2**2)


def _choose_factor(factors, x, v, prec=200):
    """SymPy's ``_choose_factor`` over QQ: the factor that vanishes at
    ``v``, found by evaluation at increasing precision."""
    from ... import S
    if isinstance(factors[0], tuple):
        factors = [f[0] for f in factors]
    if len(factors) == 1:
        return factors[0]
    illegal = (S.NaN, S.Infinity, S.NegativeInfinity, S.ComplexInfinity)
    prec1 = 10
    while prec1 <= prec:
        fe = [f.as_expr().xreplace({x: v}) for f in factors]
        if v.is_number:
            fe = [f.n(prec) for f in fe]
        candidates = [(abs(f.n(prec1)), i) for i, f in enumerate(fe)]
        if not any(i in illegal for i, _ in candidates):
            can = sorted(candidates)
            (a, ix), (b, _) = can[:2]
            if b > a*10**6:
                return factors[ix]
        prec1 *= 2
    raise NotImplementedError("multiple candidates for the minimal polynomial of %s" % v)


def _minimal_polynomial_sq(p, n, x):
    """SymPy's minimal polynomial of ``p**(1/n)`` for a sum of surds ``p``
    (None when ``p`` is not one)."""
    from ... import Poly, Rational, factor_list, sympify
    p = sympify(p)
    n = sympify(n)
    if not n.is_Integer or not n > 0 or not _is_sum_surds(p):
        return None
    pn = p**Rational(1, n)
    p -= x
    while 1:
        p1 = _separate_sq(p)
        if p1 is p:
            p = p1.subs({x: x**n})
            break
        else:
            p = p1
    if n == 1:
        p1 = Poly(p)
        if p.coeff(x**p1.degree(x)) < 0:
            p = -p
        p = p.primitive()[1]
        return p
    factors = factor_list(p)[1]
    return _choose_factor(factors, x, pn)
