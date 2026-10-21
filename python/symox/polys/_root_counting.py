"""Sturm counting of real roots over coefficient domains other than ZZ, QQ.

This follows the reference's ``dup_count_real_roots`` step by step: the
polynomial moves to the field of its domain, its Sturm sequence is built from
the square-free part by negated remainders, and sign variations are counted
at both bounds. Exact domains give exact counts; over ``RR`` the remainders
are computed in the same binary floating-point arithmetic as the reference,
including its removal of leading terms created by rounding.
"""
from fractions import Fraction

import mpmath
from mpmath.libmp import to_rational as _mpf_to_rational

from .polyerrors import DomainError as _DomainError, PolynomialDivisionFailed


class _Ops:
    """Field arithmetic on dense lists, leading coefficient first."""

    exact = True

    def __init__(self, zero, one, quo, from_int):
        self.zero, self.one, self.quo, self.from_int = zero, one, quo, from_int

    def strip(self, f):
        k = 0
        for c in f:
            if c:
                break
            k += 1
        return f[k:]

    def degree(self, f):
        return len(f) - 1 if f else float('-inf')

    def neg(self, f):
        return [-c for c in f]

    def sub(self, f, g):
        if not f:
            return self.neg(g)
        if not g:
            return f
        df, dg = self.degree(f), self.degree(g)
        if df == dg:
            return self.strip([a - b for a, b in zip(f, g)])
        k = abs(df - dg)
        if df > dg:
            h, f = f[:k], f[k:]
        else:
            h, g = self.neg(g[:k]), g[k:]
        return h + [a - b for a, b in zip(f, g)]

    def add_term(self, f, c, i):
        if not c:
            return f
        n = len(f)
        m = n - i - 1
        if i == n - 1:
            return self.strip([f[0] + c] + f[1:])
        if i >= n:
            return [c] + [self.zero]*(i - n) + f
        return f[:m] + [f[m] + c] + f[m + 1:]

    def mul_term(self, f, c, i):
        if not c or not f:
            return []
        return [cf*c for cf in f] + [self.zero]*i

    def div(self, f, g):
        df, dg = self.degree(f), self.degree(g)
        q, r, dr = [], f, df
        if not g:
            raise ZeroDivisionError("polynomial division")
        if df < dg:
            return q, r
        lc_g = g[0]
        while True:
            lc_r = r[0]
            c = self.quo(lc_r, lc_g)
            j = dr - dg
            q = self.add_term(q, c, j)
            h = self.mul_term(g, c, j)
            r = self.sub(r, h)
            _dr, dr = dr, self.degree(r)
            if dr < dg:
                break
            elif dr == _dr and not self.exact:
                # remove leading term created by rounding error
                r = self.strip(r[1:])
                dr = self.degree(r)
                if dr < dg:
                    break
            elif not (dr < _dr):
                raise PolynomialDivisionFailed(f, g, None)
        return q, r

    def monic(self, f):
        if not f:
            return f
        lc = f[0]
        if lc == self.one:
            return f
        return [self.quo(c, lc) for c in f]

    def diff(self, f):
        n = self.degree(f)
        if n < 1:
            return []
        deriv = []
        for coeff in f[:-1]:
            deriv.append(self.from_int(n)*coeff)
            n -= 1
        return self.strip(deriv)

    def gcd(self, f, g):
        """The monic greatest common divisor, unique over a field."""
        while g:
            f, g = g, self.div(f, g)[1]
        return self.monic(f)

    def eval(self, f, a):
        if not a:
            return f[-1] if f else self.zero
        result = self.zero
        for c in f:
            result *= a
            result += c
        return result


class _RealOps(_Ops):
    """The reference's ``RR``: binary floating point at the field precision."""

    exact = False

    def __init__(self, prec):
        ctx = mpmath.MPContext()
        ctx.prec = prec
        self.ctx = ctx
        self.max_denom = max(2**prec // 200, 99)
        super().__init__(ctx.mpf(0), ctx.mpf(1), lambda a, b: a/b, ctx.mpf)

    def to_rational(self, a):
        """``RR.to_rational`` with the default denominator limit."""
        p, q = _mpf_to_rational(a._mpf_)
        p, q = int(p), int(q)
        if q <= self.max_denom:
            return Fraction(p, q)
        p0, q0, p1, q1 = 0, 1, 1, 0
        n, d = p, q
        while True:
            k = n//d
            q2 = q0 + k*q1
            if q2 > self.max_denom:
                break
            p0, q0, p1, q1 = p1, q1, p0 + k*p1, q2
            n, d = d, n - k*d
        k = (self.max_denom - q0)//q1
        number = Fraction(p, q)
        bound1 = Fraction(p0 + k*p1, q0 + k*q1)
        bound2 = Fraction(p1, q1)
        if not bound2 or not bound1:
            return Fraction(p, q)
        if abs(bound2 - number) <= abs(bound1 - number):
            return bound2
        return bound1

    def from_rational(self, a):
        return self.ctx.mpf(a.numerator)/a.denominator

    def gcd(self, f, g):
        # Inexact gcds are taken over the exact domain, as in the reference.
        exact = _Ops(Fraction(0), Fraction(1), lambda a, b: a/b, Fraction)
        f = exact.strip([self.to_rational(c) for c in f])
        g = exact.strip([self.to_rational(c) for c in g])
        h = exact.gcd(f, g)
        return self.strip([self.from_rational(c) for c in h])

    def eval(self, f, a):
        if not a:
            return f[-1] if f else self.zero
        raise TypeError("unsupported operand type(s) for *=: 'mpf' and 'PythonMPQ'")


def _real(ctx, value):
    """The binary floating-point value of an ``RR`` element."""
    mpf = getattr(value, '_mpf_', None)
    if mpf is not None:
        return ctx.make_mpf(mpf)
    return ctx.mpf(int(value.p))/int(value.q)


def _comparable(value):
    flag = getattr(value, 'is_comparable', None)
    if flag is None:
        flag = value.as_expr().is_comparable
    return bool(flag)


def _sign_rule(K):
    """The reference's admissible domains for sign variation counting."""
    if K.is_RR:
        return True
    if K.is_AlgebraicField and _comparable(K.ext):
        return True
    if ((K.is_PolynomialRing or K.is_FractionField) and len(K.symbols) == 1 and
            (K.dom.is_ZZ or K.dom.is_QQ or K.is_AlgebraicField) and
            K.symbols[0].is_transcendental and K.symbols[0].is_comparable):
        return True
    return False


def count_real_roots(f, inf, sup):
    """Distinct real roots of the univariate ``f`` in ``[inf, sup]``.

    ``inf`` and ``sup`` are rational numbers or ``None`` for an infinite end.
    """
    if f.degree() <= 0:
        return 0
    K = f.domain
    R = K
    if not K.is_Field:
        K = K.get_field()
    if not _sign_rule(K):
        raise _DomainError("sign variation counting not supported over %s" % K)
    coeffs = f.rep.to_list()
    if K.is_RR:
        ops = _RealOps(getattr(K, 'prec', 53))
        rep = [_real(ops.ctx, K.to_sympy(c)) for c in coeffs]
        is_negative = lambda a: a < 0
        point = lambda a: a
    else:
        ops = _Ops(K.zero, K.one, K.quo, K.convert)
        rep = [K.convert(c, R) for c in coeffs] if R is not K else list(coeffs)
        to_sympy = K.to_sympy
        is_negative = lambda a: bool(a) and bool(to_sympy(a) < 0)
        point = K.convert
    rep = ops.strip(rep)

    # dup_sqf_part
    g = rep
    if is_negative(g[0]):
        g = ops.neg(g)
    g = ops.monic(ops.div(g, ops.gcd(g, ops.diff(g)))[0])

    sturm = [g, ops.diff(g)]
    while sturm[-1]:
        s = ops.div(sturm[-2], sturm[-1])[1]
        sturm.append(ops.neg(s))
    sturm = sturm[:-1]

    def variations(values):
        prev, k = ops.zero, 0
        for coeff in values:
            if is_negative(coeff*prev):
                k += 1
            if coeff:
                prev = coeff
        return k

    if inf is None:
        signs_inf = variations([s[0]*(-1)**ops.degree(s) for s in sturm])
    else:
        a = point(inf)
        signs_inf = variations([ops.eval(s, a) for s in sturm])
    if sup is None:
        signs_sup = variations([s[0] for s in sturm])
    else:
        b = point(sup)
        signs_sup = variations([ops.eval(s, b) for s in sturm])
    count = abs(signs_inf - signs_sup)
    if inf is not None and not ops.eval(rep, point(inf)):
        count += 1
    return count


def _rebuild(f, K, ops, rep):
    """A polynomial over ``K`` in the generator of ``f`` from field values."""
    from .. import Float, Poly
    if K.is_RR:
        prec = getattr(K, 'prec', 53)
        coeffs = [K.from_sympy(Float(c, precision=prec)) for c in rep]
    else:
        coeffs = list(rep)
    return Poly.from_list(coeffs or [K.zero], *f.gens, domain=K)


def sturm(f, auto=True):
    """The reference's ``Poly.sturm`` outside ZZ and QQ.

    The Sturm sequence of the square-free part over the field of ``f``: every
    remainder is exact over exact fields, and over ``RR`` it is computed in
    binary floating point like the reference.
    """
    if auto and not f.domain.is_Field:
        f = f.to_field()
    if len(f.gens) != 1:
        raise ValueError('univariate polynomial expected')
    K = f.domain
    if not K.is_Field:
        raise _DomainError("Cannot compute Sturm sequence over %s" % K)
    coeffs = f.rep.to_list()
    if K.is_RR:
        ops = _RealOps(getattr(K, 'prec', 53))
        rep = ops.strip([_real(ops.ctx, K.to_sympy(c)) for c in coeffs])
    else:
        ops = _Ops(K.zero, K.one, K.quo, K.convert)
        rep = ops.strip(list(coeffs))
    if K.is_FiniteField:
        # The radical, as the product of the square-free factors (one for
        # the zero polynomial): in positive characteristic f/gcd(f, f')
        # loses factors whose multiplicity the characteristic divides.
        g = [K.one]
        for factor, _ in f.sqf_list()[1]:
            g = _poly_mul(ops, g, ops.monic(list(factor.rep.to_list())))
    elif not rep:
        g = rep
    else:
        g = ops.monic(ops.div(rep, ops.gcd(rep, ops.diff(rep)))[0])
    sequence = [g, ops.diff(g)]
    while sequence[-1]:
        sequence.append(ops.neg(ops.div(sequence[-2], sequence[-1])[1]))
    return [_rebuild(f, K, ops, s) for s in sequence[:-1]]


def _poly_mul(ops, f, g):
    if not f or not g:
        return []
    h = [ops.zero]*(len(f) + len(g) - 1)
    for i, a in enumerate(f):
        for j, b in enumerate(g):
            h[i + j] = h[i + j] + a*b
    return h
