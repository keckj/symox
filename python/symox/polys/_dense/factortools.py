"""Generic dense polynomial routines over an arbitrary ground domain.

Ported from the reference SymPy 1.14 ``sympy.polys.factortools`` (BSD license,
Copyright (c) SymPy Development Team) with docstrings removed. Symox uses
them only for the ring Z/nZ of ``GF(n)`` moduli outside the native prime
field path, so the reference code paths (and their errors on zero divisors)
are followed exactly.
"""
from ._support import GROUND_TYPES
from ._support import _randint
from .galoistools import gf_from_int_poly, gf_to_int_poly, gf_lshift, gf_add_mul, gf_mul, gf_div, gf_rem, gf_gcdex, gf_sqf_p, gf_factor_sqf, gf_factor
from .densebasic import dup_LC, dmp_LC, dmp_ground_LC, dup_TC, dup_convert, dmp_convert, dup_degree, dmp_degree, dmp_degree_in, dmp_degree_list, dmp_from_dict, dmp_zero_p, dmp_one, dmp_nest, dmp_raise, dup_strip, dmp_ground, dup_inflate, dmp_exclude, dmp_include, dmp_inject, dmp_eject, dup_terms_gcd, dmp_terms_gcd
from .densearith import dup_neg, dmp_neg, dup_add, dmp_add, dup_sub, dmp_sub, dup_mul, dmp_mul, dup_sqr, dmp_pow, dup_div, dmp_div, dup_quo, dmp_quo, dmp_expand, dmp_add_mul, dup_sub_mul, dmp_sub_mul, dup_lshift, dup_max_norm, dmp_max_norm, dup_l1_norm, dup_mul_ground, dmp_mul_ground, dup_quo_ground, dmp_quo_ground
from .densetools import dup_clear_denoms, dmp_clear_denoms, dup_trunc, dmp_ground_trunc, dup_content, dup_monic, dmp_ground_monic, dup_primitive, dmp_ground_primitive, dmp_eval_tail, dmp_eval_in, dmp_diff_eval_in, dup_shift, dmp_shift, dup_mirror
from .euclidtools import dmp_primitive, dup_inner_gcd, dmp_inner_gcd
from .sqfreetools import dup_sqf_p, dup_sqf_norm, dmp_sqf_norm, dup_sqf_part, dmp_sqf_part, _dup_check_degrees, _dmp_check_degrees
from ._support import _sort_factors
from ._support import query
from ..polyerrors import ExtraneousFactors, DomainError, CoercionFailed, EvaluationFailed
from ...utilities.iterables import subsets
from math import ceil as _ceil, log as _log, log2 as _log2
if GROUND_TYPES == 'flint':
    from flint import fmpz_poly
else:
    fmpz_poly = None

def dup_trial_division(f, factors, K):
    result = []
    for factor in factors:
        k = 0
        while True:
            (q, r) = dup_div(f, factor, K)
            if not r:
                (f, k) = (q, k + 1)
            else:
                break
        if k == 0:
            raise RuntimeError('trial division failed')
        result.append((factor, k))
    return _sort_factors(result)

def dmp_trial_division(f, factors, u, K):
    result = []
    for factor in factors:
        k = 0
        while True:
            (q, r) = dmp_div(f, factor, u, K)
            if dmp_zero_p(r, u):
                (f, k) = (q, k + 1)
            else:
                break
        if k == 0:
            raise RuntimeError('trial division failed')
        result.append((factor, k))
    return _sort_factors(result)

def dup_zz_mignotte_bound(f, K):
    from ... import binomial
    d = dup_degree(f)
    delta = _ceil(d / 2)
    delta2 = _ceil(delta / 2)
    eucl_norm = K.sqrt(sum((cf ** 2 for cf in f)))
    t1 = binomial(delta - 1, delta2)
    t2 = binomial(delta - 1, delta2 - 1)
    lc = K.abs(dup_LC(f, K))
    bound = t1 * eucl_norm + t2 * lc
    bound += dup_max_norm(f, K)
    bound = _ceil(bound / 2) * 2
    return bound

def dmp_zz_mignotte_bound(f, u, K):
    a = dmp_max_norm(f, u, K)
    b = abs(dmp_ground_LC(f, u, K))
    n = sum(dmp_degree_list(f, u))
    return K.sqrt(K(n + 1)) * 2 ** n * a * b

def dup_zz_hensel_step(m, f, g, h, s, t, K):
    M = m ** 2
    e = dup_sub_mul(f, g, h, K)
    e = dup_trunc(e, M, K)
    (q, r) = dup_div(dup_mul(s, e, K), h, K)
    q = dup_trunc(q, M, K)
    r = dup_trunc(r, M, K)
    u = dup_add(dup_mul(t, e, K), dup_mul(q, g, K), K)
    G = dup_trunc(dup_add(g, u, K), M, K)
    H = dup_trunc(dup_add(h, r, K), M, K)
    u = dup_add(dup_mul(s, G, K), dup_mul(t, H, K), K)
    b = dup_trunc(dup_sub(u, [K.one], K), M, K)
    (c, d) = dup_div(dup_mul(s, b, K), H, K)
    c = dup_trunc(c, M, K)
    d = dup_trunc(d, M, K)
    u = dup_add(dup_mul(t, b, K), dup_mul(c, G, K), K)
    S = dup_trunc(dup_sub(s, d, K), M, K)
    T = dup_trunc(dup_sub(t, u, K), M, K)
    return (G, H, S, T)

def dup_zz_hensel_lift(p, f, f_list, l, K):
    r = len(f_list)
    lc = dup_LC(f, K)
    if r == 1:
        F = dup_mul_ground(f, K.gcdex(lc, p ** l)[0], K)
        return [dup_trunc(F, p ** l, K)]
    m = p
    k = r // 2
    d = int(_ceil(_log2(l)))
    g = gf_from_int_poly([lc], p)
    for f_i in f_list[:k]:
        g = gf_mul(g, gf_from_int_poly(f_i, p), p, K)
    h = gf_from_int_poly(f_list[k], p)
    for f_i in f_list[k + 1:]:
        h = gf_mul(h, gf_from_int_poly(f_i, p), p, K)
    (s, t, _) = gf_gcdex(g, h, p, K)
    g = gf_to_int_poly(g, p)
    h = gf_to_int_poly(h, p)
    s = gf_to_int_poly(s, p)
    t = gf_to_int_poly(t, p)
    for _ in range(1, d + 1):
        ((g, h, s, t), m) = (dup_zz_hensel_step(m, f, g, h, s, t, K), m ** 2)
    return dup_zz_hensel_lift(p, g, f_list[:k], l, K) + dup_zz_hensel_lift(p, h, f_list[k:], l, K)

def _test_pl(fc, q, pl):
    if q > pl // 2:
        q = q - pl
    if not q:
        return True
    return fc % q == 0

def dup_zz_zassenhaus(f, K):
    n = dup_degree(f)
    if n == 1:
        return [f]
    from ...ntheory import isprime
    fc = f[-1]
    A = dup_max_norm(f, K)
    b = dup_LC(f, K)
    B = int(abs(K.sqrt(K(n + 1)) * 2 ** n * A * b))
    C = int((n + 1) ** (2 * n) * A ** (2 * n - 1))
    gamma = int(_ceil(2 * _log2(C)))
    bound = int(2 * gamma * _log(gamma))
    a = []
    for px in range(3, bound + 1):
        if not isprime(px) or b % px == 0:
            continue
        px = K.convert(px)
        F = gf_from_int_poly(f, px)
        if not gf_sqf_p(F, px, K):
            continue
        fsqfx = gf_factor_sqf(F, px, K)[1]
        a.append((px, fsqfx))
        if len(fsqfx) < 15 or len(a) > 4:
            break
    (p, fsqf) = min(a, key=lambda x: len(x[1]))
    l = int(_ceil(_log(2 * B + 1, p)))
    modular = [gf_to_int_poly(ff, p) for ff in fsqf]
    g = dup_zz_hensel_lift(p, f, modular, l, K)
    sorted_T = range(len(g))
    T = set(sorted_T)
    (factors, s) = ([], 1)
    pl = p ** l
    while 2 * s <= len(T):
        for S in subsets(sorted_T, s):
            if b == 1:
                q = 1
                for i in S:
                    q = q * g[i][-1]
                q = q % pl
                if not _test_pl(fc, q, pl):
                    continue
            else:
                G = [b]
                for i in S:
                    G = dup_mul(G, g[i], K)
                G = dup_trunc(G, pl, K)
                G = dup_primitive(G, K)[1]
                q = G[-1]
                if q and fc % q != 0:
                    continue
            H = [b]
            S = set(S)
            T_S = T - S
            if b == 1:
                G = [b]
                for i in S:
                    G = dup_mul(G, g[i], K)
                G = dup_trunc(G, pl, K)
            for i in T_S:
                H = dup_mul(H, g[i], K)
            H = dup_trunc(H, pl, K)
            G_norm = dup_l1_norm(G, K)
            H_norm = dup_l1_norm(H, K)
            if G_norm * H_norm <= B:
                T = T_S
                sorted_T = [i for i in sorted_T if i not in S]
                G = dup_primitive(G, K)[1]
                f = dup_primitive(H, K)[1]
                factors.append(G)
                b = dup_LC(f, K)
                break
        else:
            s += 1
    return factors + [f]

def dup_zz_irreducible_p(f, K):
    lc = dup_LC(f, K)
    tc = dup_TC(f, K)
    e_fc = dup_content(f[1:], K)
    if e_fc:
        from ...ntheory import factorint
        e_ff = factorint(int(e_fc))
        for p in e_ff.keys():
            if lc % p and tc % p ** 2:
                return True

def dup_cyclotomic_p(f, K, irreducible=False):
    if K.is_QQ:
        try:
            (K0, K) = (K, K.get_ring())
            f = dup_convert(f, K0, K)
        except CoercionFailed:
            return False
    elif not K.is_ZZ:
        return False
    lc = dup_LC(f, K)
    tc = dup_TC(f, K)
    if lc != 1 or (tc != -1 and tc != 1):
        return False
    if not irreducible:
        (coeff, factors) = dup_factor_list(f, K)
        if coeff != K.one or factors != [(f, 1)]:
            return False
    n = dup_degree(f)
    (g, h) = ([], [])
    for i in range(n, -1, -2):
        g.insert(0, f[i])
    for i in range(n - 1, -1, -2):
        h.insert(0, f[i])
    g = dup_sqr(dup_strip(g), K)
    h = dup_sqr(dup_strip(h), K)
    F = dup_sub(g, dup_lshift(h, 1, K), K)
    if K.is_negative(dup_LC(F, K)):
        F = dup_neg(F, K)
    if F == f:
        return True
    g = dup_mirror(f, K)
    if K.is_negative(dup_LC(g, K)):
        g = dup_neg(g, K)
    if F == g and dup_cyclotomic_p(g, K):
        return True
    G = dup_sqf_part(F, K)
    if dup_sqr(G, K) == F and dup_cyclotomic_p(G, K):
        return True
    return False

def dup_zz_cyclotomic_poly(n, K):
    from ...ntheory import factorint
    h = [K.one, -K.one]
    for (p, k) in factorint(n).items():
        h = dup_quo(dup_inflate(h, p, K), h, K)
        h = dup_inflate(h, p ** (k - 1), K)
    return h

def _dup_cyclotomic_decompose(n, K):
    from ...ntheory import factorint
    H = [[K.one, -K.one]]
    for (p, k) in factorint(n).items():
        Q = [dup_quo(dup_inflate(h, p, K), h, K) for h in H]
        H.extend(Q)
        for i in range(1, k):
            Q = [dup_inflate(q, p, K) for q in Q]
            H.extend(Q)
    return H

def dup_zz_cyclotomic_factor(f, K):
    (lc_f, tc_f) = (dup_LC(f, K), dup_TC(f, K))
    if dup_degree(f) <= 0:
        return None
    if lc_f != 1 or tc_f not in [-1, 1]:
        return None
    if any((bool(cf) for cf in f[1:-1])):
        return None
    n = dup_degree(f)
    F = _dup_cyclotomic_decompose(n, K)
    if not K.is_one(tc_f):
        return F
    else:
        H = []
        for h in _dup_cyclotomic_decompose(2 * n, K):
            if h not in F:
                H.append(h)
        return H

def dup_zz_factor_sqf(f, K):
    (cont, g) = dup_primitive(f, K)
    n = dup_degree(g)
    if dup_LC(g, K) < 0:
        (cont, g) = (-cont, dup_neg(g, K))
    if n <= 0:
        return (cont, [])
    elif n == 1:
        return (cont, [g])
    if query('USE_IRREDUCIBLE_IN_FACTOR'):
        if dup_zz_irreducible_p(g, K):
            return (cont, [g])
    factors = None
    if query('USE_CYCLOTOMIC_FACTOR'):
        factors = dup_zz_cyclotomic_factor(g, K)
    if factors is None:
        factors = dup_zz_zassenhaus(g, K)
    return (cont, _sort_factors(factors, multiple=False))

def dup_zz_factor(f, K):
    if GROUND_TYPES == 'flint':
        f_flint = fmpz_poly(f[::-1])
        (cont, factors) = f_flint.factor()
        factors = [(fac.coeffs()[::-1], exp) for (fac, exp) in factors]
        return (cont, _sort_factors(factors))
    (cont, g) = dup_primitive(f, K)
    n = dup_degree(g)
    if dup_LC(g, K) < 0:
        (cont, g) = (-cont, dup_neg(g, K))
    if n <= 0:
        return (cont, [])
    elif n == 1:
        return (cont, [(g, 1)])
    if query('USE_IRREDUCIBLE_IN_FACTOR'):
        if dup_zz_irreducible_p(g, K):
            return (cont, [(g, 1)])
    g = dup_sqf_part(g, K)
    H = None
    if query('USE_CYCLOTOMIC_FACTOR'):
        H = dup_zz_cyclotomic_factor(g, K)
    if H is None:
        H = dup_zz_zassenhaus(g, K)
    factors = dup_trial_division(f, H, K)
    _dup_check_degrees(f, factors)
    return (cont, factors)

def dmp_zz_wang_non_divisors(E, cs, ct, K):
    result = [cs * ct]
    for q in E:
        q = abs(q)
        for r in reversed(result):
            while r != 1:
                r = K.gcd(r, q)
                q = q // r
            if K.is_one(q):
                return None
        result.append(q)
    return result[1:]

def dmp_zz_wang_test_points(f, T, ct, A, u, K):
    if not dmp_eval_tail(dmp_LC(f, K), A, u - 1, K):
        raise EvaluationFailed('no luck')
    g = dmp_eval_tail(f, A, u, K)
    if not dup_sqf_p(g, K):
        raise EvaluationFailed('no luck')
    (c, h) = dup_primitive(g, K)
    if K.is_negative(dup_LC(h, K)):
        (c, h) = (-c, dup_neg(h, K))
    v = u - 1
    E = [dmp_eval_tail(t, A, v, K) for (t, _) in T]
    D = dmp_zz_wang_non_divisors(E, c, ct, K)
    if D is not None:
        return (c, h, E)
    else:
        raise EvaluationFailed('no luck')

def dmp_zz_wang_lead_coeffs(f, T, cs, E, H, A, u, K):
    (C, J, v) = ([], [0] * len(E), u - 1)
    for h in H:
        c = dmp_one(v, K)
        d = dup_LC(h, K) * cs
        for i in reversed(range(len(E))):
            (k, e, (t, _)) = (0, E[i], T[i])
            while not d % e:
                (d, k) = (d // e, k + 1)
            if k != 0:
                (c, J[i]) = (dmp_mul(c, dmp_pow(t, k, v, K), v, K), 1)
        C.append(c)
    if not all(J):
        raise ExtraneousFactors
    (CC, HH) = ([], [])
    for (c, h) in zip(C, H):
        d = dmp_eval_tail(c, A, v, K)
        lc = dup_LC(h, K)
        if K.is_one(cs):
            cc = lc // d
        else:
            g = K.gcd(lc, d)
            (d, cc) = (d // g, lc // g)
            (h, cs) = (dup_mul_ground(h, d, K), cs // d)
        c = dmp_mul_ground(c, cc, v, K)
        CC.append(c)
        HH.append(h)
    if K.is_one(cs):
        return (f, HH, CC)
    (CCC, HHH) = ([], [])
    for (c, h) in zip(CC, HH):
        CCC.append(dmp_mul_ground(c, cs, v, K))
        HHH.append(dmp_mul_ground(h, cs, 0, K))
    f = dmp_mul_ground(f, cs ** (len(H) - 1), u, K)
    return (f, HHH, CCC)

def dup_zz_diophantine(F, m, p, K):
    if len(F) == 2:
        (a, b) = F
        f = gf_from_int_poly(a, p)
        g = gf_from_int_poly(b, p)
        (s, t, G) = gf_gcdex(g, f, p, K)
        s = gf_lshift(s, m, K)
        t = gf_lshift(t, m, K)
        (q, s) = gf_div(s, f, p, K)
        t = gf_add_mul(t, q, g, p, K)
        s = gf_to_int_poly(s, p)
        t = gf_to_int_poly(t, p)
        result = [s, t]
    else:
        G = [F[-1]]
        for f in reversed(F[1:-1]):
            G.insert(0, dup_mul(f, G[0], K))
        (S, T) = ([], [[1]])
        for (f, g) in zip(F, G):
            (t, s) = dmp_zz_diophantine([g, f], T[-1], [], 0, p, 1, K)
            T.append(t)
            S.append(s)
        (result, S) = ([], S + [T[-1]])
        for (s, f) in zip(S, F):
            s = gf_from_int_poly(s, p)
            f = gf_from_int_poly(f, p)
            r = gf_rem(gf_lshift(s, m, K), f, p, K)
            s = gf_to_int_poly(r, p)
            result.append(s)
    return result

def dmp_zz_diophantine(F, c, A, d, p, u, K):
    if not A:
        S = [[] for _ in F]
        n = dup_degree(c)
        for (i, coeff) in enumerate(c):
            if not coeff:
                continue
            T = dup_zz_diophantine(F, n - i, p, K)
            for (j, (s, t)) in enumerate(zip(S, T)):
                t = dup_mul_ground(t, coeff, K)
                S[j] = dup_trunc(dup_add(s, t, K), p, K)
    else:
        n = len(A)
        e = dmp_expand(F, u, K)
        (a, A) = (A[-1], A[:-1])
        (B, G) = ([], [])
        for f in F:
            B.append(dmp_quo(e, f, u, K))
            G.append(dmp_eval_in(f, a, n, u, K))
        C = dmp_eval_in(c, a, n, u, K)
        v = u - 1
        S = dmp_zz_diophantine(G, C, A, d, p, v, K)
        S = [dmp_raise(s, 1, v, K) for s in S]
        for (s, b) in zip(S, B):
            c = dmp_sub_mul(c, s, b, u, K)
        c = dmp_ground_trunc(c, p, u, K)
        m = dmp_nest([K.one, -a], n, K)
        M = dmp_one(n, K)
        for k in range(0, d):
            if dmp_zero_p(c, u):
                break
            M = dmp_mul(M, m, u, K)
            C = dmp_diff_eval_in(c, k + 1, a, n, u, K)
            if not dmp_zero_p(C, v):
                C = dmp_quo_ground(C, K.factorial(K(k) + 1), v, K)
                T = dmp_zz_diophantine(G, C, A, d, p, v, K)
                for (i, t) in enumerate(T):
                    T[i] = dmp_mul(dmp_raise(t, 1, v, K), M, u, K)
                for (i, (s, t)) in enumerate(zip(S, T)):
                    S[i] = dmp_add(s, t, u, K)
                for (t, b) in zip(T, B):
                    c = dmp_sub_mul(c, t, b, u, K)
                c = dmp_ground_trunc(c, p, u, K)
        S = [dmp_ground_trunc(s, p, u, K) for s in S]
    return S

def dmp_zz_wang_hensel_lifting(f, H, LC, A, p, u, K):
    (S, n, v) = ([f], len(A), u - 1)
    H = list(H)
    for (i, a) in enumerate(reversed(A[1:])):
        s = dmp_eval_in(S[0], a, n - i, u - i, K)
        S.insert(0, dmp_ground_trunc(s, p, v - i, K))
    d = max(dmp_degree_list(f, u)[1:])
    for (j, s, a) in zip(range(2, n + 2), S, A):
        (G, w) = (list(H), j - 1)
        (I, J) = (A[:j - 2], A[j - 1:])
        for (i, (h, lc)) in enumerate(zip(H, LC)):
            lc = dmp_ground_trunc(dmp_eval_tail(lc, J, v, K), p, w - 1, K)
            H[i] = [lc] + dmp_raise(h[1:], 1, w - 1, K)
        m = dmp_nest([K.one, -a], w, K)
        M = dmp_one(w, K)
        c = dmp_sub(s, dmp_expand(H, w, K), w, K)
        dj = dmp_degree_in(s, w, w)
        for k in range(0, dj):
            if dmp_zero_p(c, w):
                break
            M = dmp_mul(M, m, w, K)
            C = dmp_diff_eval_in(c, k + 1, a, w, w, K)
            if not dmp_zero_p(C, w - 1):
                C = dmp_quo_ground(C, K.factorial(K(k) + 1), w - 1, K)
                T = dmp_zz_diophantine(G, C, I, d, p, w - 1, K)
                for (i, (h, t)) in enumerate(zip(H, T)):
                    h = dmp_add_mul(h, dmp_raise(t, 1, w - 1, K), M, w, K)
                    H[i] = dmp_ground_trunc(h, p, w, K)
                h = dmp_sub(s, dmp_expand(H, w, K), w, K)
                c = dmp_ground_trunc(h, p, w, K)
    if dmp_expand(H, u, K) != f:
        raise ExtraneousFactors
    else:
        return H

def dmp_zz_wang(f, u, K, mod=None, seed=None):
    from ...ntheory import nextprime
    randint = _randint(seed)
    (ct, T) = dmp_zz_factor(dmp_LC(f, K), u - 1, K)
    b = dmp_zz_mignotte_bound(f, u, K)
    p = K(nextprime(b))
    if mod is None:
        if u == 1:
            mod = 2
        else:
            mod = 1
    (history, configs, A, r) = (set(), [], [K.zero] * u, None)
    try:
        (cs, s, E) = dmp_zz_wang_test_points(f, T, ct, A, u, K)
        (_, H) = dup_zz_factor_sqf(s, K)
        r = len(H)
        if r == 1:
            return [f]
        configs = [(s, cs, E, H, A)]
    except EvaluationFailed:
        pass
    eez_num_configs = query('EEZ_NUMBER_OF_CONFIGS')
    eez_num_tries = query('EEZ_NUMBER_OF_TRIES')
    eez_mod_step = query('EEZ_MODULUS_STEP')
    while len(configs) < eez_num_configs:
        for _ in range(eez_num_tries):
            A = [K(randint(-mod, mod)) for _ in range(u)]
            if tuple(A) not in history:
                history.add(tuple(A))
            else:
                continue
            try:
                (cs, s, E) = dmp_zz_wang_test_points(f, T, ct, A, u, K)
            except EvaluationFailed:
                continue
            (_, H) = dup_zz_factor_sqf(s, K)
            rr = len(H)
            if r is not None:
                if rr != r:
                    if rr < r:
                        (configs, r) = ([], rr)
                    else:
                        continue
            else:
                r = rr
            if r == 1:
                return [f]
            configs.append((s, cs, E, H, A))
            if len(configs) == eez_num_configs:
                break
        else:
            mod += eez_mod_step
    (s_norm, s_arg, i) = (None, 0, 0)
    for (s, _, _, _, _) in configs:
        _s_norm = dup_max_norm(s, K)
        if s_norm is not None:
            if _s_norm < s_norm:
                s_norm = _s_norm
                s_arg = i
        else:
            s_norm = _s_norm
        i += 1
    (_, cs, E, H, A) = configs[s_arg]
    orig_f = f
    try:
        (f, H, LC) = dmp_zz_wang_lead_coeffs(f, T, cs, E, H, A, u, K)
        factors = dmp_zz_wang_hensel_lifting(f, H, LC, A, p, u, K)
    except ExtraneousFactors:
        if query('EEZ_RESTART_IF_NEEDED'):
            return dmp_zz_wang(orig_f, u, K, mod + 1)
        else:
            raise ExtraneousFactors('we need to restart algorithm with better parameters')
    result = []
    for f in factors:
        (_, f) = dmp_ground_primitive(f, u, K)
        if K.is_negative(dmp_ground_LC(f, u, K)):
            f = dmp_neg(f, u, K)
        result.append(f)
    return result

def dmp_zz_factor(f, u, K):
    if not u:
        return dup_zz_factor(f, K)
    if dmp_zero_p(f, u):
        return (K.zero, [])
    (cont, g) = dmp_ground_primitive(f, u, K)
    if dmp_ground_LC(g, u, K) < 0:
        (cont, g) = (-cont, dmp_neg(g, u, K))
    if all((d <= 0 for d in dmp_degree_list(g, u))):
        return (cont, [])
    (G, g) = dmp_primitive(g, u, K)
    factors = []
    if dmp_degree(g, u) > 0:
        g = dmp_sqf_part(g, u, K)
        H = dmp_zz_wang(g, u, K)
        factors = dmp_trial_division(f, H, u, K)
    for (g, k) in dmp_zz_factor(G, u - 1, K)[1]:
        factors.insert(0, ([g], k))
    _dmp_check_degrees(f, u, factors)
    return (cont, _sort_factors(factors))

def dup_qq_i_factor(f, K0):
    K1 = K0.as_AlgebraicField()
    f = dup_convert(f, K0, K1)
    (coeff, factors) = dup_factor_list(f, K1)
    factors = [(dup_convert(fac, K1, K0), i) for (fac, i) in factors]
    coeff = K0.convert(coeff, K1)
    return (coeff, factors)

def dup_zz_i_factor(f, K0):
    K1 = K0.get_field()
    f = dup_convert(f, K0, K1)
    (coeff, factors) = dup_qq_i_factor(f, K1)
    new_factors = []
    for (fac, i) in factors:
        (fac_denom, fac_num) = dup_clear_denoms(fac, K1)
        fac_num_ZZ_I = dup_convert(fac_num, K1, K0)
        (content, fac_prim) = dmp_ground_primitive(fac_num_ZZ_I, 0, K0)
        coeff = coeff * content ** i // fac_denom ** i
        new_factors.append((fac_prim, i))
    factors = new_factors
    coeff = K0.convert(coeff, K1)
    return (coeff, factors)

def dmp_qq_i_factor(f, u, K0):
    K1 = K0.as_AlgebraicField()
    f = dmp_convert(f, u, K0, K1)
    (coeff, factors) = dmp_factor_list(f, u, K1)
    factors = [(dmp_convert(fac, u, K1, K0), i) for (fac, i) in factors]
    coeff = K0.convert(coeff, K1)
    return (coeff, factors)

def dmp_zz_i_factor(f, u, K0):
    K1 = K0.get_field()
    f = dmp_convert(f, u, K0, K1)
    (coeff, factors) = dmp_qq_i_factor(f, u, K1)
    new_factors = []
    for (fac, i) in factors:
        (fac_denom, fac_num) = dmp_clear_denoms(fac, u, K1)
        fac_num_ZZ_I = dmp_convert(fac_num, u, K1, K0)
        (content, fac_prim) = dmp_ground_primitive(fac_num_ZZ_I, u, K0)
        coeff = coeff * content ** i // fac_denom ** i
        new_factors.append((fac_prim, i))
    factors = new_factors
    coeff = K0.convert(coeff, K1)
    return (coeff, factors)

def dup_ext_factor(f, K):
    (n, lc) = (dup_degree(f), dup_LC(f, K))
    f = dup_monic(f, K)
    if n <= 0:
        return (lc, [])
    if n == 1:
        return (lc, [(f, 1)])
    (f, F) = (dup_sqf_part(f, K), f)
    (s, g, r) = dup_sqf_norm(f, K)
    factors = dup_factor_list_include(r, K.dom)
    if len(factors) == 1:
        return (lc, [(f, n // dup_degree(f))])
    H = s * K.unit
    for (i, (factor, _)) in enumerate(factors):
        h = dup_convert(factor, K.dom, K)
        (h, _, g) = dup_inner_gcd(h, g, K)
        h = dup_shift(h, H, K)
        factors[i] = h
    factors = dup_trial_division(F, factors, K)
    _dup_check_degrees(F, factors)
    return (lc, factors)

def dmp_ext_factor(f, u, K):
    if not u:
        return dup_ext_factor(f, K)
    lc = dmp_ground_LC(f, u, K)
    f = dmp_ground_monic(f, u, K)
    if all((d <= 0 for d in dmp_degree_list(f, u))):
        return (lc, [])
    (f, F) = (dmp_sqf_part(f, u, K), f)
    (s, g, r) = dmp_sqf_norm(f, u, K)
    factors = dmp_factor_list_include(r, u, K.dom)
    if len(factors) == 1:
        factors = [f]
    else:
        for (i, (factor, _)) in enumerate(factors):
            h = dmp_convert(factor, u, K.dom, K)
            (h, _, g) = dmp_inner_gcd(h, g, u, K)
            a = [si * K.unit for si in s]
            h = dmp_shift(h, a, u, K)
            factors[i] = h
    result = dmp_trial_division(F, factors, u, K)
    _dmp_check_degrees(F, u, result)
    return (lc, result)

def dup_gf_factor(f, K):
    f = dup_convert(f, K, K.dom)
    (coeff, factors) = gf_factor(f, K.mod, K.dom)
    for (i, (f, k)) in enumerate(factors):
        factors[i] = (dup_convert(f, K.dom, K), k)
    return (K.convert(coeff, K.dom), factors)

def dmp_gf_factor(f, u, K):
    raise NotImplementedError('multivariate polynomials over finite fields')

def dup_factor_list(f, K0):
    (j, f) = dup_terms_gcd(f, K0)
    (cont, f) = dup_primitive(f, K0)
    if K0.is_FiniteField:
        (coeff, factors) = dup_gf_factor(f, K0)
    elif K0.is_Algebraic:
        (coeff, factors) = dup_ext_factor(f, K0)
    elif K0.is_GaussianRing:
        (coeff, factors) = dup_zz_i_factor(f, K0)
    elif K0.is_GaussianField:
        (coeff, factors) = dup_qq_i_factor(f, K0)
    else:
        if not K0.is_Exact:
            (K0_inexact, K0) = (K0, K0.get_exact())
            f = dup_convert(f, K0_inexact, K0)
        else:
            K0_inexact = None
        if K0.is_Field:
            K = K0.get_ring()
            (denom, f) = dup_clear_denoms(f, K0, K)
            f = dup_convert(f, K0, K)
        else:
            K = K0
        if K.is_ZZ:
            (coeff, factors) = dup_zz_factor(f, K)
        elif K.is_Poly:
            (f, u) = dmp_inject(f, 0, K)
            (coeff, factors) = dmp_factor_list(f, u, K.dom)
            for (i, (f, k)) in enumerate(factors):
                factors[i] = (dmp_eject(f, u, K), k)
            coeff = K.convert(coeff, K.dom)
        else:
            raise DomainError('factorization not supported over %s' % K0)
        if K0.is_Field:
            for (i, (f, k)) in enumerate(factors):
                factors[i] = (dup_convert(f, K, K0), k)
            coeff = K0.convert(coeff, K)
            coeff = K0.quo(coeff, denom)
            if K0_inexact:
                for (i, (f, k)) in enumerate(factors):
                    max_norm = dup_max_norm(f, K0)
                    f = dup_quo_ground(f, max_norm, K0)
                    f = dup_convert(f, K0, K0_inexact)
                    factors[i] = (f, k)
                    coeff = K0.mul(coeff, K0.pow(max_norm, k))
                coeff = K0_inexact.convert(coeff, K0)
                K0 = K0_inexact
    if j:
        factors.insert(0, ([K0.one, K0.zero], j))
    return (coeff * cont, _sort_factors(factors))

def dup_factor_list_include(f, K):
    (coeff, factors) = dup_factor_list(f, K)
    if not factors:
        return [(dup_strip([coeff]), 1)]
    else:
        g = dup_mul_ground(factors[0][0], coeff, K)
        return [(g, factors[0][1])] + factors[1:]

def dmp_factor_list(f, u, K0):
    if not u:
        return dup_factor_list(f, K0)
    (J, f) = dmp_terms_gcd(f, u, K0)
    (cont, f) = dmp_ground_primitive(f, u, K0)
    if K0.is_FiniteField:
        (coeff, factors) = dmp_gf_factor(f, u, K0)
    elif K0.is_Algebraic:
        (coeff, factors) = dmp_ext_factor(f, u, K0)
    elif K0.is_GaussianRing:
        (coeff, factors) = dmp_zz_i_factor(f, u, K0)
    elif K0.is_GaussianField:
        (coeff, factors) = dmp_qq_i_factor(f, u, K0)
    else:
        if not K0.is_Exact:
            (K0_inexact, K0) = (K0, K0.get_exact())
            f = dmp_convert(f, u, K0_inexact, K0)
        else:
            K0_inexact = None
        if K0.is_Field:
            K = K0.get_ring()
            (denom, f) = dmp_clear_denoms(f, u, K0, K)
            f = dmp_convert(f, u, K0, K)
        else:
            K = K0
        if K.is_ZZ:
            (levels, f, v) = dmp_exclude(f, u, K)
            (coeff, factors) = dmp_zz_factor(f, v, K)
            for (i, (f, k)) in enumerate(factors):
                factors[i] = (dmp_include(f, levels, v, K), k)
        elif K.is_Poly:
            (f, v) = dmp_inject(f, u, K)
            (coeff, factors) = dmp_factor_list(f, v, K.dom)
            for (i, (f, k)) in enumerate(factors):
                factors[i] = (dmp_eject(f, v, K), k)
            coeff = K.convert(coeff, K.dom)
        else:
            raise DomainError('factorization not supported over %s' % K0)
        if K0.is_Field:
            for (i, (f, k)) in enumerate(factors):
                factors[i] = (dmp_convert(f, u, K, K0), k)
            coeff = K0.convert(coeff, K)
            coeff = K0.quo(coeff, denom)
            if K0_inexact:
                for (i, (f, k)) in enumerate(factors):
                    max_norm = dmp_max_norm(f, u, K0)
                    f = dmp_quo_ground(f, max_norm, u, K0)
                    f = dmp_convert(f, u, K0, K0_inexact)
                    factors[i] = (f, k)
                    coeff = K0.mul(coeff, K0.pow(max_norm, k))
                coeff = K0_inexact.convert(coeff, K0)
                K0 = K0_inexact
    for (i, j) in enumerate(reversed(J)):
        if not j:
            continue
        term = {(0,) * (u - i) + (1,) + (0,) * i: K0.one}
        factors.insert(0, (dmp_from_dict(term, u, K0), j))
    return (coeff * cont, _sort_factors(factors))

def dmp_factor_list_include(f, u, K):
    if not u:
        return dup_factor_list_include(f, K)
    (coeff, factors) = dmp_factor_list(f, u, K)
    if not factors:
        return [(dmp_ground(coeff, u), 1)]
    else:
        g = dmp_mul_ground(factors[0][0], coeff, u, K)
        return [(g, factors[0][1])] + factors[1:]

def dup_irreducible_p(f, K):
    return dmp_irreducible_p(f, 0, K)

def dmp_irreducible_p(f, u, K):
    (_, factors) = dmp_factor_list(f, u, K)
    if not factors:
        return True
    elif len(factors) > 1:
        return False
    else:
        (_, k) = factors[0]
        return k == 1
