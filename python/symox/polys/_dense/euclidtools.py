"""Generic dense polynomial routines over an arbitrary ground domain.

Ported from the reference SymPy 1.14 ``sympy.polys.euclidtools`` (BSD license,
Copyright (c) SymPy Development Team) with docstrings removed. Symox uses
them only for the ring Z/nZ of ``GF(n)`` moduli outside the native prime
field path, so the reference code paths (and their errors on zero divisors)
are followed exactly.
"""
from .densearith import dup_sub_mul, dup_neg, dmp_neg, dmp_add, dmp_sub, dup_mul, dmp_mul, dmp_pow, dup_div, dmp_div, dup_rem, dup_quo, dmp_quo, dup_prem, dmp_prem, dup_mul_ground, dmp_mul_ground, dmp_mul_term, dup_quo_ground, dmp_quo_ground, dup_max_norm, dmp_max_norm
from .densebasic import dup_strip, dmp_raise, dmp_zero, dmp_one, dmp_ground, dmp_one_p, dmp_zero_p, dmp_zeros, dup_degree, dmp_degree, dmp_degree_in, dup_LC, dmp_LC, dmp_ground_LC, dmp_multi_deflate, dmp_inflate, dup_convert, dmp_convert, dmp_apply_pairs
from .densetools import dup_clear_denoms, dmp_clear_denoms, dup_diff, dmp_diff, dup_eval, dmp_eval, dmp_eval_in, dup_trunc, dmp_ground_trunc, dup_monic, dmp_ground_monic, dup_primitive, dmp_ground_primitive, dup_extract, dmp_ground_extract
from .galoistools import gf_int, gf_crt
from ._support import query
from ..polyerrors import MultivariatePolynomialError, HeuristicGCDFailed, HomomorphismFailed, NotInvertible, DomainError

def dup_half_gcdex(f, g, K):
    if not K.is_Field:
        raise DomainError('Cannot compute half extended GCD over %s' % K)
    (a, b) = ([K.one], [])
    while g:
        (q, r) = dup_div(f, g, K)
        (f, g) = (g, r)
        (a, b) = (b, dup_sub_mul(a, q, b, K))
    a = dup_quo_ground(a, dup_LC(f, K), K)
    f = dup_monic(f, K)
    return (a, f)

def dmp_half_gcdex(f, g, u, K):
    if not u:
        return dup_half_gcdex(f, g, K)
    else:
        raise MultivariatePolynomialError(f, g)

def dup_gcdex(f, g, K):
    (s, h) = dup_half_gcdex(f, g, K)
    F = dup_sub_mul(h, s, f, K)
    t = dup_quo(F, g, K)
    return (s, t, h)

def dmp_gcdex(f, g, u, K):
    if not u:
        return dup_gcdex(f, g, K)
    else:
        raise MultivariatePolynomialError(f, g)

def dup_invert(f, g, K):
    (s, h) = dup_half_gcdex(f, g, K)
    if h == [K.one]:
        return dup_rem(s, g, K)
    else:
        raise NotInvertible('zero divisor')

def dmp_invert(f, g, u, K):
    if not u:
        return dup_invert(f, g, K)
    else:
        raise MultivariatePolynomialError(f, g)

def dup_euclidean_prs(f, g, K):
    prs = [f, g]
    h = dup_rem(f, g, K)
    while h:
        prs.append(h)
        (f, g) = (g, h)
        h = dup_rem(f, g, K)
    return prs

def dmp_euclidean_prs(f, g, u, K):
    if not u:
        return dup_euclidean_prs(f, g, K)
    else:
        raise MultivariatePolynomialError(f, g)

def dup_primitive_prs(f, g, K):
    prs = [f, g]
    (_, h) = dup_primitive(dup_prem(f, g, K), K)
    while h:
        prs.append(h)
        (f, g) = (g, h)
        (_, h) = dup_primitive(dup_prem(f, g, K), K)
    return prs

def dmp_primitive_prs(f, g, u, K):
    if not u:
        return dup_primitive_prs(f, g, K)
    else:
        raise MultivariatePolynomialError(f, g)

def dup_inner_subresultants(f, g, K):
    n = dup_degree(f)
    m = dup_degree(g)
    if n < m:
        (f, g) = (g, f)
        (n, m) = (m, n)
    if not f:
        return ([], [])
    if not g:
        return ([f], [K.one])
    R = [f, g]
    d = n - m
    b = (-K.one) ** (d + 1)
    h = dup_prem(f, g, K)
    h = dup_mul_ground(h, b, K)
    lc = dup_LC(g, K)
    c = lc ** d
    S = [K.one, c]
    c = -c
    while h:
        k = dup_degree(h)
        R.append(h)
        (f, g, m, d) = (g, h, k, m - k)
        b = -lc * c ** d
        h = dup_prem(f, g, K)
        h = dup_quo_ground(h, b, K)
        lc = dup_LC(g, K)
        if d > 1:
            q = c ** (d - 1)
            c = K.quo((-lc) ** d, q)
        else:
            c = -lc
        S.append(-c)
    return (R, S)

def dup_subresultants(f, g, K):
    return dup_inner_subresultants(f, g, K)[0]

def dup_prs_resultant(f, g, K):
    if not f or not g:
        return (K.zero, [])
    (R, S) = dup_inner_subresultants(f, g, K)
    if dup_degree(R[-1]) > 0:
        return (K.zero, R)
    return (S[-1], R)

def dup_resultant(f, g, K, includePRS=False):
    if includePRS:
        return dup_prs_resultant(f, g, K)
    return dup_prs_resultant(f, g, K)[0]

def dmp_inner_subresultants(f, g, u, K):
    if not u:
        return dup_inner_subresultants(f, g, K)
    n = dmp_degree(f, u)
    m = dmp_degree(g, u)
    if n < m:
        (f, g) = (g, f)
        (n, m) = (m, n)
    if dmp_zero_p(f, u):
        return ([], [])
    v = u - 1
    if dmp_zero_p(g, u):
        return ([f], [dmp_ground(K.one, v)])
    R = [f, g]
    d = n - m
    b = dmp_pow(dmp_ground(-K.one, v), d + 1, v, K)
    h = dmp_prem(f, g, u, K)
    h = dmp_mul_term(h, b, 0, u, K)
    lc = dmp_LC(g, K)
    c = dmp_pow(lc, d, v, K)
    S = [dmp_ground(K.one, v), c]
    c = dmp_neg(c, v, K)
    while not dmp_zero_p(h, u):
        k = dmp_degree(h, u)
        R.append(h)
        (f, g, m, d) = (g, h, k, m - k)
        b = dmp_mul(dmp_neg(lc, v, K), dmp_pow(c, d, v, K), v, K)
        h = dmp_prem(f, g, u, K)
        h = [dmp_quo(ch, b, v, K) for ch in h]
        lc = dmp_LC(g, K)
        if d > 1:
            p = dmp_pow(dmp_neg(lc, v, K), d, v, K)
            q = dmp_pow(c, d - 1, v, K)
            c = dmp_quo(p, q, v, K)
        else:
            c = dmp_neg(lc, v, K)
        S.append(dmp_neg(c, v, K))
    return (R, S)

def dmp_subresultants(f, g, u, K):
    return dmp_inner_subresultants(f, g, u, K)[0]

def dmp_prs_resultant(f, g, u, K):
    if not u:
        return dup_prs_resultant(f, g, K)
    if dmp_zero_p(f, u) or dmp_zero_p(g, u):
        return (dmp_zero(u - 1), [])
    (R, S) = dmp_inner_subresultants(f, g, u, K)
    if dmp_degree(R[-1], u) > 0:
        return (dmp_zero(u - 1), R)
    return (S[-1], R)

def dmp_zz_modular_resultant(f, g, p, u, K):
    if not u:
        return gf_int(dup_prs_resultant(f, g, K)[0] % p, p)
    v = u - 1
    n = dmp_degree(f, u)
    m = dmp_degree(g, u)
    N = dmp_degree_in(f, 1, u)
    M = dmp_degree_in(g, 1, u)
    B = n * M + m * N
    (D, a) = ([K.one], -K.one)
    r = dmp_zero(v)
    while dup_degree(D) <= B:
        while True:
            a += K.one
            if a == p:
                raise HomomorphismFailed('no luck')
            F = dmp_eval_in(f, gf_int(a, p), 1, u, K)
            if dmp_degree(F, v) == n:
                G = dmp_eval_in(g, gf_int(a, p), 1, u, K)
                if dmp_degree(G, v) == m:
                    break
        R = dmp_zz_modular_resultant(F, G, p, v, K)
        e = dmp_eval(r, a, v, K)
        if not v:
            R = dup_strip([R])
            e = dup_strip([e])
        else:
            R = [R]
            e = [e]
        d = K.invert(dup_eval(D, a, K), p)
        d = dup_mul_ground(D, d, K)
        d = dmp_raise(d, v, 0, K)
        c = dmp_mul(d, dmp_sub(R, e, v, K), v, K)
        r = dmp_add(r, c, v, K)
        r = dmp_ground_trunc(r, p, v, K)
        D = dup_mul(D, [K.one, -a], K)
        D = dup_trunc(D, p, K)
    return r

def _collins_crt(r, R, P, p, K):
    return gf_int(gf_crt([r, R], [P, p], K), P * p)

def dmp_zz_collins_resultant(f, g, u, K):
    n = dmp_degree(f, u)
    m = dmp_degree(g, u)
    if n < 0 or m < 0:
        return dmp_zero(u - 1)
    A = dmp_max_norm(f, u, K)
    B = dmp_max_norm(g, u, K)
    a = dmp_ground_LC(f, u, K)
    b = dmp_ground_LC(g, u, K)
    v = u - 1
    B = K(2) * K.factorial(K(n + m)) * A ** m * B ** n
    (r, p, P) = (dmp_zero(v), K.one, K.one)
    from ...ntheory import nextprime
    while P <= B:
        p = K(nextprime(p))
        while not a % p or not b % p:
            p = K(nextprime(p))
        F = dmp_ground_trunc(f, p, u, K)
        G = dmp_ground_trunc(g, p, u, K)
        try:
            R = dmp_zz_modular_resultant(F, G, p, u, K)
        except HomomorphismFailed:
            continue
        if K.is_one(P):
            r = R
        else:
            r = dmp_apply_pairs(r, R, _collins_crt, (P, p, K), v, K)
        P *= p
    return r

def dmp_qq_collins_resultant(f, g, u, K0):
    n = dmp_degree(f, u)
    m = dmp_degree(g, u)
    if n < 0 or m < 0:
        return dmp_zero(u - 1)
    K1 = K0.get_ring()
    (cf, f) = dmp_clear_denoms(f, u, K0, K1)
    (cg, g) = dmp_clear_denoms(g, u, K0, K1)
    f = dmp_convert(f, u, K0, K1)
    g = dmp_convert(g, u, K0, K1)
    r = dmp_zz_collins_resultant(f, g, u, K1)
    r = dmp_convert(r, u - 1, K1, K0)
    c = K0.convert(cf ** m * cg ** n, K1)
    return dmp_quo_ground(r, c, u - 1, K0)

def dmp_resultant(f, g, u, K, includePRS=False):
    if not u:
        return dup_resultant(f, g, K, includePRS=includePRS)
    if includePRS:
        return dmp_prs_resultant(f, g, u, K)
    if K.is_Field:
        if K.is_QQ and query('USE_COLLINS_RESULTANT'):
            return dmp_qq_collins_resultant(f, g, u, K)
    elif K.is_ZZ and query('USE_COLLINS_RESULTANT'):
        return dmp_zz_collins_resultant(f, g, u, K)
    return dmp_prs_resultant(f, g, u, K)[0]

def dup_discriminant(f, K):
    d = dup_degree(f)
    if d <= 0:
        return K.zero
    else:
        s = (-1) ** (d * (d - 1) // 2)
        c = dup_LC(f, K)
        r = dup_resultant(f, dup_diff(f, 1, K), K)
        return K.quo(r, c * K(s))

def dmp_discriminant(f, u, K):
    if not u:
        return dup_discriminant(f, K)
    (d, v) = (dmp_degree(f, u), u - 1)
    if d <= 0:
        return dmp_zero(v)
    else:
        s = (-1) ** (d * (d - 1) // 2)
        c = dmp_LC(f, K)
        r = dmp_resultant(f, dmp_diff(f, 1, u, K), u, K)
        c = dmp_mul_ground(c, K(s), v, K)
        return dmp_quo(r, c, v, K)

def _dup_rr_trivial_gcd(f, g, K):
    if not (f or g):
        return ([], [], [])
    elif not f:
        if K.is_nonnegative(dup_LC(g, K)):
            return (g, [], [K.one])
        else:
            return (dup_neg(g, K), [], [-K.one])
    elif not g:
        if K.is_nonnegative(dup_LC(f, K)):
            return (f, [K.one], [])
        else:
            return (dup_neg(f, K), [-K.one], [])
    return None

def _dup_ff_trivial_gcd(f, g, K):
    if not (f or g):
        return ([], [], [])
    elif not f:
        return (dup_monic(g, K), [], [dup_LC(g, K)])
    elif not g:
        return (dup_monic(f, K), [dup_LC(f, K)], [])
    else:
        return None

def _dmp_rr_trivial_gcd(f, g, u, K):
    zero_f = dmp_zero_p(f, u)
    zero_g = dmp_zero_p(g, u)
    if_contain_one = dmp_one_p(f, u, K) or dmp_one_p(g, u, K)
    if zero_f and zero_g:
        return tuple(dmp_zeros(3, u, K))
    elif zero_f:
        if K.is_nonnegative(dmp_ground_LC(g, u, K)):
            return (g, dmp_zero(u), dmp_one(u, K))
        else:
            return (dmp_neg(g, u, K), dmp_zero(u), dmp_ground(-K.one, u))
    elif zero_g:
        if K.is_nonnegative(dmp_ground_LC(f, u, K)):
            return (f, dmp_one(u, K), dmp_zero(u))
        else:
            return (dmp_neg(f, u, K), dmp_ground(-K.one, u), dmp_zero(u))
    elif if_contain_one:
        return (dmp_one(u, K), f, g)
    elif query('USE_SIMPLIFY_GCD'):
        return _dmp_simplify_gcd(f, g, u, K)
    else:
        return None

def _dmp_ff_trivial_gcd(f, g, u, K):
    zero_f = dmp_zero_p(f, u)
    zero_g = dmp_zero_p(g, u)
    if zero_f and zero_g:
        return tuple(dmp_zeros(3, u, K))
    elif zero_f:
        return (dmp_ground_monic(g, u, K), dmp_zero(u), dmp_ground(dmp_ground_LC(g, u, K), u))
    elif zero_g:
        return (dmp_ground_monic(f, u, K), dmp_ground(dmp_ground_LC(f, u, K), u), dmp_zero(u))
    elif query('USE_SIMPLIFY_GCD'):
        return _dmp_simplify_gcd(f, g, u, K)
    else:
        return None

def _dmp_simplify_gcd(f, g, u, K):
    df = dmp_degree(f, u)
    dg = dmp_degree(g, u)
    if df > 0 and dg > 0:
        return None
    if not (df or dg):
        F = dmp_LC(f, K)
        G = dmp_LC(g, K)
    elif not df:
        F = dmp_LC(f, K)
        G = dmp_content(g, u, K)
    else:
        F = dmp_content(f, u, K)
        G = dmp_LC(g, K)
    v = u - 1
    h = dmp_gcd(F, G, v, K)
    cff = [dmp_quo(cf, h, v, K) for cf in f]
    cfg = [dmp_quo(cg, h, v, K) for cg in g]
    return ([h], cff, cfg)

def dup_rr_prs_gcd(f, g, K):
    result = _dup_rr_trivial_gcd(f, g, K)
    if result is not None:
        return result
    (fc, F) = dup_primitive(f, K)
    (gc, G) = dup_primitive(g, K)
    c = K.gcd(fc, gc)
    h = dup_subresultants(F, G, K)[-1]
    (_, h) = dup_primitive(h, K)
    c *= K.canonical_unit(dup_LC(h, K))
    h = dup_mul_ground(h, c, K)
    cff = dup_quo(f, h, K)
    cfg = dup_quo(g, h, K)
    return (h, cff, cfg)

def dup_ff_prs_gcd(f, g, K):
    result = _dup_ff_trivial_gcd(f, g, K)
    if result is not None:
        return result
    h = dup_subresultants(f, g, K)[-1]
    h = dup_monic(h, K)
    cff = dup_quo(f, h, K)
    cfg = dup_quo(g, h, K)
    return (h, cff, cfg)

def dmp_rr_prs_gcd(f, g, u, K):
    if not u:
        return dup_rr_prs_gcd(f, g, K)
    result = _dmp_rr_trivial_gcd(f, g, u, K)
    if result is not None:
        return result
    (fc, F) = dmp_primitive(f, u, K)
    (gc, G) = dmp_primitive(g, u, K)
    h = dmp_subresultants(F, G, u, K)[-1]
    (c, _, _) = dmp_rr_prs_gcd(fc, gc, u - 1, K)
    (_, h) = dmp_primitive(h, u, K)
    h = dmp_mul_term(h, c, 0, u, K)
    unit = K.canonical_unit(dmp_ground_LC(h, u, K))
    if unit != K.one:
        h = dmp_mul_ground(h, unit, u, K)
    cff = dmp_quo(f, h, u, K)
    cfg = dmp_quo(g, h, u, K)
    return (h, cff, cfg)

def dmp_ff_prs_gcd(f, g, u, K):
    if not u:
        return dup_ff_prs_gcd(f, g, K)
    result = _dmp_ff_trivial_gcd(f, g, u, K)
    if result is not None:
        return result
    (fc, F) = dmp_primitive(f, u, K)
    (gc, G) = dmp_primitive(g, u, K)
    h = dmp_subresultants(F, G, u, K)[-1]
    (c, _, _) = dmp_ff_prs_gcd(fc, gc, u - 1, K)
    (_, h) = dmp_primitive(h, u, K)
    h = dmp_mul_term(h, c, 0, u, K)
    h = dmp_ground_monic(h, u, K)
    cff = dmp_quo(f, h, u, K)
    cfg = dmp_quo(g, h, u, K)
    return (h, cff, cfg)
HEU_GCD_MAX = 6

def _dup_zz_gcd_interpolate(h, x, K):
    f = []
    while h:
        g = h % x
        if g > x // 2:
            g -= x
        f.insert(0, g)
        h = (h - g) // x
    return f

def dup_zz_heu_gcd(f, g, K):
    result = _dup_rr_trivial_gcd(f, g, K)
    if result is not None:
        return result
    df = dup_degree(f)
    dg = dup_degree(g)
    (gcd, f, g) = dup_extract(f, g, K)
    if df == 0 or dg == 0:
        return ([gcd], f, g)
    f_norm = dup_max_norm(f, K)
    g_norm = dup_max_norm(g, K)
    B = K(2 * min(f_norm, g_norm) + 29)
    x = max(min(B, 99 * K.sqrt(B)), 2 * min(f_norm // abs(dup_LC(f, K)), g_norm // abs(dup_LC(g, K))) + 4)
    for i in range(0, HEU_GCD_MAX):
        ff = dup_eval(f, x, K)
        gg = dup_eval(g, x, K)
        if ff and gg:
            h = K.gcd(ff, gg)
            cff = ff // h
            cfg = gg // h
            h = _dup_zz_gcd_interpolate(h, x, K)
            h = dup_primitive(h, K)[1]
            (cff_, r) = dup_div(f, h, K)
            if not r:
                (cfg_, r) = dup_div(g, h, K)
                if not r:
                    h = dup_mul_ground(h, gcd, K)
                    return (h, cff_, cfg_)
            cff = _dup_zz_gcd_interpolate(cff, x, K)
            (h, r) = dup_div(f, cff, K)
            if not r:
                (cfg_, r) = dup_div(g, h, K)
                if not r:
                    h = dup_mul_ground(h, gcd, K)
                    return (h, cff, cfg_)
            cfg = _dup_zz_gcd_interpolate(cfg, x, K)
            (h, r) = dup_div(g, cfg, K)
            if not r:
                (cff_, r) = dup_div(f, h, K)
                if not r:
                    h = dup_mul_ground(h, gcd, K)
                    return (h, cff_, cfg)
        x = 73794 * x * K.sqrt(K.sqrt(x)) // 27011
    raise HeuristicGCDFailed('no luck')

def _dmp_zz_gcd_interpolate(h, x, v, K):
    f = []
    while not dmp_zero_p(h, v):
        g = dmp_ground_trunc(h, x, v, K)
        f.insert(0, g)
        h = dmp_sub(h, g, v, K)
        h = dmp_quo_ground(h, x, v, K)
    if K.is_negative(dmp_ground_LC(f, v + 1, K)):
        return dmp_neg(f, v + 1, K)
    else:
        return f

def dmp_zz_heu_gcd(f, g, u, K):
    if not u:
        return dup_zz_heu_gcd(f, g, K)
    result = _dmp_rr_trivial_gcd(f, g, u, K)
    if result is not None:
        return result
    (gcd, f, g) = dmp_ground_extract(f, g, u, K)
    f_norm = dmp_max_norm(f, u, K)
    g_norm = dmp_max_norm(g, u, K)
    B = K(2 * min(f_norm, g_norm) + 29)
    x = max(min(B, 99 * K.sqrt(B)), 2 * min(f_norm // abs(dmp_ground_LC(f, u, K)), g_norm // abs(dmp_ground_LC(g, u, K))) + 4)
    for i in range(0, HEU_GCD_MAX):
        ff = dmp_eval(f, x, u, K)
        gg = dmp_eval(g, x, u, K)
        v = u - 1
        if not (dmp_zero_p(ff, v) or dmp_zero_p(gg, v)):
            (h, cff, cfg) = dmp_zz_heu_gcd(ff, gg, v, K)
            h = _dmp_zz_gcd_interpolate(h, x, v, K)
            h = dmp_ground_primitive(h, u, K)[1]
            (cff_, r) = dmp_div(f, h, u, K)
            if dmp_zero_p(r, u):
                (cfg_, r) = dmp_div(g, h, u, K)
                if dmp_zero_p(r, u):
                    h = dmp_mul_ground(h, gcd, u, K)
                    return (h, cff_, cfg_)
            cff = _dmp_zz_gcd_interpolate(cff, x, v, K)
            (h, r) = dmp_div(f, cff, u, K)
            if dmp_zero_p(r, u):
                (cfg_, r) = dmp_div(g, h, u, K)
                if dmp_zero_p(r, u):
                    h = dmp_mul_ground(h, gcd, u, K)
                    return (h, cff, cfg_)
            cfg = _dmp_zz_gcd_interpolate(cfg, x, v, K)
            (h, r) = dmp_div(g, cfg, u, K)
            if dmp_zero_p(r, u):
                (cff_, r) = dmp_div(f, h, u, K)
                if dmp_zero_p(r, u):
                    h = dmp_mul_ground(h, gcd, u, K)
                    return (h, cff_, cfg)
        x = 73794 * x * K.sqrt(K.sqrt(x)) // 27011
    raise HeuristicGCDFailed('no luck')

def dup_qq_heu_gcd(f, g, K0):
    result = _dup_ff_trivial_gcd(f, g, K0)
    if result is not None:
        return result
    K1 = K0.get_ring()
    (cf, f) = dup_clear_denoms(f, K0, K1)
    (cg, g) = dup_clear_denoms(g, K0, K1)
    f = dup_convert(f, K0, K1)
    g = dup_convert(g, K0, K1)
    (h, cff, cfg) = dup_zz_heu_gcd(f, g, K1)
    h = dup_convert(h, K1, K0)
    c = dup_LC(h, K0)
    h = dup_monic(h, K0)
    cff = dup_convert(cff, K1, K0)
    cfg = dup_convert(cfg, K1, K0)
    cff = dup_mul_ground(cff, K0.quo(c, cf), K0)
    cfg = dup_mul_ground(cfg, K0.quo(c, cg), K0)
    return (h, cff, cfg)

def dmp_qq_heu_gcd(f, g, u, K0):
    result = _dmp_ff_trivial_gcd(f, g, u, K0)
    if result is not None:
        return result
    K1 = K0.get_ring()
    (cf, f) = dmp_clear_denoms(f, u, K0, K1)
    (cg, g) = dmp_clear_denoms(g, u, K0, K1)
    f = dmp_convert(f, u, K0, K1)
    g = dmp_convert(g, u, K0, K1)
    (h, cff, cfg) = dmp_zz_heu_gcd(f, g, u, K1)
    h = dmp_convert(h, u, K1, K0)
    c = dmp_ground_LC(h, u, K0)
    h = dmp_ground_monic(h, u, K0)
    cff = dmp_convert(cff, u, K1, K0)
    cfg = dmp_convert(cfg, u, K1, K0)
    cff = dmp_mul_ground(cff, K0.quo(c, cf), u, K0)
    cfg = dmp_mul_ground(cfg, K0.quo(c, cg), u, K0)
    return (h, cff, cfg)

def dup_inner_gcd(f, g, K):
    if K.is_RR or K.is_CC:
        try:
            exact = K.get_exact()
        except DomainError:
            return ([K.one], f, g)
        f = dup_convert(f, K, exact)
        g = dup_convert(g, K, exact)
        (h, cff, cfg) = dup_inner_gcd(f, g, exact)
        h = dup_convert(h, exact, K)
        cff = dup_convert(cff, exact, K)
        cfg = dup_convert(cfg, exact, K)
        return (h, cff, cfg)
    elif K.is_Field:
        if K.is_QQ and query('USE_HEU_GCD'):
            try:
                return dup_qq_heu_gcd(f, g, K)
            except HeuristicGCDFailed:
                pass
        return dup_ff_prs_gcd(f, g, K)
    else:
        if K.is_ZZ and query('USE_HEU_GCD'):
            try:
                return dup_zz_heu_gcd(f, g, K)
            except HeuristicGCDFailed:
                pass
        return dup_rr_prs_gcd(f, g, K)

def _dmp_inner_gcd(f, g, u, K):
    if not K.is_Exact:
        try:
            exact = K.get_exact()
        except DomainError:
            return (dmp_one(u, K), f, g)
        f = dmp_convert(f, u, K, exact)
        g = dmp_convert(g, u, K, exact)
        (h, cff, cfg) = _dmp_inner_gcd(f, g, u, exact)
        h = dmp_convert(h, u, exact, K)
        cff = dmp_convert(cff, u, exact, K)
        cfg = dmp_convert(cfg, u, exact, K)
        return (h, cff, cfg)
    elif K.is_Field:
        if K.is_QQ and query('USE_HEU_GCD'):
            try:
                return dmp_qq_heu_gcd(f, g, u, K)
            except HeuristicGCDFailed:
                pass
        return dmp_ff_prs_gcd(f, g, u, K)
    else:
        if K.is_ZZ and query('USE_HEU_GCD'):
            try:
                return dmp_zz_heu_gcd(f, g, u, K)
            except HeuristicGCDFailed:
                pass
        return dmp_rr_prs_gcd(f, g, u, K)

def dmp_inner_gcd(f, g, u, K):
    if not u:
        return dup_inner_gcd(f, g, K)
    (J, (f, g)) = dmp_multi_deflate((f, g), u, K)
    (h, cff, cfg) = _dmp_inner_gcd(f, g, u, K)
    return (dmp_inflate(h, J, u, K), dmp_inflate(cff, J, u, K), dmp_inflate(cfg, J, u, K))

def dup_gcd(f, g, K):
    return dup_inner_gcd(f, g, K)[0]

def dmp_gcd(f, g, u, K):
    return dmp_inner_gcd(f, g, u, K)[0]

def dup_rr_lcm(f, g, K):
    if not f or not g:
        return dmp_zero(0)
    (fc, f) = dup_primitive(f, K)
    (gc, g) = dup_primitive(g, K)
    c = K.lcm(fc, gc)
    h = dup_quo(dup_mul(f, g, K), dup_gcd(f, g, K), K)
    u = K.canonical_unit(dup_LC(h, K))
    return dup_mul_ground(h, c * u, K)

def dup_ff_lcm(f, g, K):
    h = dup_quo(dup_mul(f, g, K), dup_gcd(f, g, K), K)
    return dup_monic(h, K)

def dup_lcm(f, g, K):
    if K.is_Field:
        return dup_ff_lcm(f, g, K)
    else:
        return dup_rr_lcm(f, g, K)

def dmp_rr_lcm(f, g, u, K):
    (fc, f) = dmp_ground_primitive(f, u, K)
    (gc, g) = dmp_ground_primitive(g, u, K)
    c = K.lcm(fc, gc)
    h = dmp_quo(dmp_mul(f, g, u, K), dmp_gcd(f, g, u, K), u, K)
    return dmp_mul_ground(h, c, u, K)

def dmp_ff_lcm(f, g, u, K):
    h = dmp_quo(dmp_mul(f, g, u, K), dmp_gcd(f, g, u, K), u, K)
    return dmp_ground_monic(h, u, K)

def dmp_lcm(f, g, u, K):
    if not u:
        return dup_lcm(f, g, K)
    if K.is_Field:
        return dmp_ff_lcm(f, g, u, K)
    else:
        return dmp_rr_lcm(f, g, u, K)

def dmp_content(f, u, K):
    (cont, v) = (dmp_LC(f, K), u - 1)
    if dmp_zero_p(f, u):
        return cont
    for c in f[1:]:
        cont = dmp_gcd(cont, c, v, K)
        if dmp_one_p(cont, v, K):
            break
    if K.is_negative(dmp_ground_LC(cont, v, K)):
        return dmp_neg(cont, v, K)
    else:
        return cont

def dmp_primitive(f, u, K):
    (cont, v) = (dmp_content(f, u, K), u - 1)
    if dmp_zero_p(f, u) or dmp_one_p(cont, v, K):
        return (cont, f)
    else:
        return (cont, [dmp_quo(c, cont, v, K) for c in f])

def dup_cancel(f, g, K, include=True):
    return dmp_cancel(f, g, 0, K, include=include)

def dmp_cancel(f, g, u, K, include=True):
    K0 = None
    if K.is_Field and K.has_assoc_Ring:
        (K0, K) = (K, K.get_ring())
        (cq, f) = dmp_clear_denoms(f, u, K0, K, convert=True)
        (cp, g) = dmp_clear_denoms(g, u, K0, K, convert=True)
    else:
        (cp, cq) = (K.one, K.one)
    (_, p, q) = dmp_inner_gcd(f, g, u, K)
    if K0 is not None:
        (_, cp, cq) = K.cofactors(cp, cq)
        p = dmp_convert(p, u, K, K0)
        q = dmp_convert(q, u, K, K0)
        K = K0
    p_neg = K.is_negative(dmp_ground_LC(p, u, K))
    q_neg = K.is_negative(dmp_ground_LC(q, u, K))
    if p_neg and q_neg:
        (p, q) = (dmp_neg(p, u, K), dmp_neg(q, u, K))
    elif p_neg:
        (cp, p) = (-cp, dmp_neg(p, u, K))
    elif q_neg:
        (cp, q) = (-cp, dmp_neg(q, u, K))
    if not include:
        return (cp, cq, p, q)
    p = dmp_mul_ground(p, cp, u, K)
    q = dmp_mul_ground(q, cq, u, K)
    return (p, q)
