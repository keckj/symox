"""Generic dense polynomial routines over an arbitrary ground domain.

Ported from the reference SymPy 1.14 ``sympy.polys.sqfreetools`` (BSD license,
Copyright (c) SymPy Development Team) with docstrings removed. Symox uses
them only for the ring Z/nZ of ``GF(n)`` moduli outside the native prime
field path, so the reference code paths (and their errors on zero divisors)
are followed exactly.
"""
from .densearith import dup_neg, dmp_neg, dup_sub, dmp_sub, dup_mul, dmp_mul, dup_quo, dmp_quo, dup_mul_ground, dmp_mul_ground
from .densebasic import dup_strip, dup_LC, dmp_ground_LC, dmp_zero_p, dmp_ground, dup_degree, dmp_degree, dmp_degree_in, dmp_degree_list, dmp_raise, dmp_inject, dup_convert
from .densetools import dup_diff, dmp_diff, dmp_diff_in, dup_shift, dmp_shift, dup_monic, dmp_ground_monic, dup_primitive, dmp_ground_primitive
from .euclidtools import dup_inner_gcd, dmp_inner_gcd, dup_gcd, dmp_gcd, dmp_resultant, dmp_primitive
from .galoistools import gf_sqf_list, gf_sqf_part
from ..polyerrors import MultivariatePolynomialError, DomainError

def _dup_check_degrees(f, result):
    deg = sum((k * dup_degree(fac) for (fac, k) in result))
    assert deg == dup_degree(f)

def _dmp_check_degrees(f, u, result):
    degs = [0] * (u + 1)
    for (fac, k) in result:
        degs_fac = dmp_degree_list(fac, u)
        degs = [d1 + k * d2 for (d1, d2) in zip(degs, degs_fac)]
    assert tuple(degs) == dmp_degree_list(f, u)

def dup_sqf_p(f, K):
    if not f:
        return True
    else:
        return not dup_degree(dup_gcd(f, dup_diff(f, 1, K), K))

def dmp_sqf_p(f, u, K):
    if dmp_zero_p(f, u):
        return True
    for i in range(u + 1):
        fp = dmp_diff_in(f, 1, i, u, K)
        if dmp_zero_p(fp, u):
            continue
        gcd = dmp_gcd(f, fp, u, K)
        if dmp_degree_in(gcd, i, u) != 0:
            return False
    return True

def dup_sqf_norm(f, K):
    if not K.is_Algebraic:
        raise DomainError('ground domain must be algebraic')
    (s, g) = (0, dmp_raise(K.mod.to_list(), 1, 0, K.dom))
    while True:
        (h, _) = dmp_inject(f, 0, K, front=True)
        r = dmp_resultant(g, h, 1, K.dom)
        if dup_sqf_p(r, K.dom):
            break
        else:
            (f, s) = (dup_shift(f, -K.unit, K), s + 1)
    return (s, f, r)

def _dmp_sqf_norm_shifts(f, u, K):
    n = u + 1
    s0 = [0] * n
    yield (s0, f)
    a = K.unit
    d = dmp_degree_list(f, u)
    var_indices = [i for (di, i) in sorted(zip(d, range(u + 1))) if di > 0]
    for i in var_indices:
        s1 = s0.copy()
        s1[i] = 1
        a1 = [-a * s1i for s1i in s1]
        f1 = dmp_shift(f, a1, u, K)
        yield (s1, f1)
    j = 0
    while True:
        j += 1
        sj = [j] * n
        aj = [-a * j] * n
        fj = dmp_shift(f, aj, u, K)
        yield (sj, fj)

def dmp_sqf_norm(f, u, K):
    if not u:
        (s, g, r) = dup_sqf_norm(f, K)
        return ([s], g, r)
    if not K.is_Algebraic:
        raise DomainError('ground domain must be algebraic')
    g = dmp_raise(K.mod.to_list(), u + 1, 0, K.dom)
    for (s, f) in _dmp_sqf_norm_shifts(f, u, K):
        (h, _) = dmp_inject(f, u, K, front=True)
        r = dmp_resultant(g, h, u + 1, K.dom)
        if dmp_sqf_p(r, u, K.dom):
            break
    return (s, f, r)

def dmp_norm(f, u, K):
    if not K.is_Algebraic:
        raise DomainError('ground domain must be algebraic')
    g = dmp_raise(K.mod.to_list(), u + 1, 0, K.dom)
    (h, _) = dmp_inject(f, u, K, front=True)
    return dmp_resultant(g, h, u + 1, K.dom)

def dup_gf_sqf_part(f, K):
    f = dup_convert(f, K, K.dom)
    g = gf_sqf_part(f, K.mod, K.dom)
    return dup_convert(g, K.dom, K)

def dmp_gf_sqf_part(f, u, K):
    raise NotImplementedError('multivariate polynomials over finite fields')

def dup_sqf_part(f, K):
    if K.is_FiniteField:
        return dup_gf_sqf_part(f, K)
    if not f:
        return f
    if K.is_negative(dup_LC(f, K)):
        f = dup_neg(f, K)
    gcd = dup_gcd(f, dup_diff(f, 1, K), K)
    sqf = dup_quo(f, gcd, K)
    if K.is_Field:
        return dup_monic(sqf, K)
    else:
        return dup_primitive(sqf, K)[1]

def dmp_sqf_part(f, u, K):
    if not u:
        return dup_sqf_part(f, K)
    if K.is_FiniteField:
        return dmp_gf_sqf_part(f, u, K)
    if dmp_zero_p(f, u):
        return f
    if K.is_negative(dmp_ground_LC(f, u, K)):
        f = dmp_neg(f, u, K)
    gcd = f
    for i in range(u + 1):
        gcd = dmp_gcd(gcd, dmp_diff_in(f, 1, i, u, K), u, K)
    sqf = dmp_quo(f, gcd, u, K)
    if K.is_Field:
        return dmp_ground_monic(sqf, u, K)
    else:
        return dmp_ground_primitive(sqf, u, K)[1]

def dup_gf_sqf_list(f, K, all=False):
    f_orig = f
    f = dup_convert(f, K, K.dom)
    (coeff, factors) = gf_sqf_list(f, K.mod, K.dom, all=all)
    for (i, (f, k)) in enumerate(factors):
        factors[i] = (dup_convert(f, K.dom, K), k)
    _dup_check_degrees(f_orig, factors)
    return (K.convert(coeff, K.dom), factors)

def dmp_gf_sqf_list(f, u, K, all=False):
    raise NotImplementedError('multivariate polynomials over finite fields')

def dup_sqf_list(f, K, all=False):
    if K.is_FiniteField:
        return dup_gf_sqf_list(f, K, all=all)
    f_orig = f
    if K.is_Field:
        coeff = dup_LC(f, K)
        f = dup_monic(f, K)
    else:
        (coeff, f) = dup_primitive(f, K)
        if K.is_negative(dup_LC(f, K)):
            f = dup_neg(f, K)
            coeff = -coeff
    if dup_degree(f) <= 0:
        return (coeff, [])
    (result, i) = ([], 1)
    h = dup_diff(f, 1, K)
    (g, p, q) = dup_inner_gcd(f, h, K)
    while True:
        d = dup_diff(p, 1, K)
        h = dup_sub(q, d, K)
        if not h:
            result.append((p, i))
            break
        (g, p, q) = dup_inner_gcd(p, h, K)
        if all or dup_degree(g) > 0:
            result.append((g, i))
        i += 1
    _dup_check_degrees(f_orig, result)
    return (coeff, result)

def dup_sqf_list_include(f, K, all=False):
    (coeff, factors) = dup_sqf_list(f, K, all=all)
    if factors and factors[0][1] == 1:
        g = dup_mul_ground(factors[0][0], coeff, K)
        return [(g, 1)] + factors[1:]
    else:
        g = dup_strip([coeff])
        return [(g, 1)] + factors

def dmp_sqf_list(f, u, K, all=False):
    if not u:
        return dup_sqf_list(f, K, all=all)
    if K.is_FiniteField:
        return dmp_gf_sqf_list(f, u, K, all=all)
    f_orig = f
    if K.is_Field:
        coeff = dmp_ground_LC(f, u, K)
        f = dmp_ground_monic(f, u, K)
    else:
        (coeff, f) = dmp_ground_primitive(f, u, K)
        if K.is_negative(dmp_ground_LC(f, u, K)):
            f = dmp_neg(f, u, K)
            coeff = -coeff
    deg = dmp_degree(f, u)
    if deg < 0:
        return (coeff, [])
    (content, f) = dmp_primitive(f, u, K)
    result = {}
    if deg != 0:
        h = dmp_diff(f, 1, u, K)
        (g, p, q) = dmp_inner_gcd(f, h, u, K)
        i = 1
        while True:
            d = dmp_diff(p, 1, u, K)
            h = dmp_sub(q, d, u, K)
            if dmp_zero_p(h, u):
                result[i] = p
                break
            (g, p, q) = dmp_inner_gcd(p, h, u, K)
            if all or dmp_degree(g, u) > 0:
                result[i] = g
            i += 1
    (coeff_content, result_content) = dmp_sqf_list(content, u - 1, K, all=all)
    coeff *= coeff_content
    for (fac, i) in result_content:
        fac = [fac]
        if i in result:
            result[i] = dmp_mul(result[i], fac, u, K)
        else:
            result[i] = fac
    result = [(result[i], i) for i in sorted(result)]
    _dmp_check_degrees(f_orig, u, result)
    return (coeff, result)

def dmp_sqf_list_include(f, u, K, all=False):
    if not u:
        return dup_sqf_list_include(f, K, all=all)
    (coeff, factors) = dmp_sqf_list(f, u, K, all=all)
    if factors and factors[0][1] == 1:
        g = dmp_mul_ground(factors[0][0], coeff, u, K)
        return [(g, 1)] + factors[1:]
    else:
        g = dmp_ground(coeff, u)
        return [(g, 1)] + factors

def dup_gff_list(f, K):
    if not f:
        raise ValueError("greatest factorial factorization doesn't exist for a zero polynomial")
    f = dup_monic(f, K)
    if not dup_degree(f):
        return []
    else:
        g = dup_gcd(f, dup_shift(f, K.one, K), K)
        H = dup_gff_list(g, K)
        for (i, (h, k)) in enumerate(H):
            g = dup_mul(g, dup_shift(h, -K(k), K), K)
            H[i] = (h, k + 1)
        f = dup_quo(f, g, K)
        if not dup_degree(f):
            return H
        else:
            return [(f, 1)] + H

def dmp_gff_list(f, u, K):
    if not u:
        return dup_gff_list(f, K)
    else:
        raise MultivariatePolynomialError(f)
