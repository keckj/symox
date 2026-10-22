"""Generic dense polynomial routines over an arbitrary ground domain.

Ported from the reference SymPy 1.14 ``sympy.polys.galoistools`` (BSD license,
Copyright (c) SymPy Development Team) with docstrings removed. Symox uses
them only for the ring Z/nZ of ``GF(n)`` moduli outside the native prime
field path, so the reference code paths (and their errors on zero divisors)
are followed exactly.
"""
from math import ceil as _ceil, sqrt as _sqrt, prod
from ._support import uniform, _randint
from ._support import SYMPY_INTS, MPZ, invert
from ._support import query
from ..polyerrors import ExactQuotientFailed
from ._support import _sort_factors

def gf_crt(U, M, K=None):
    p = prod(M, start=K.one)
    v = K.zero
    for (u, m) in zip(U, M):
        e = p // m
        (s, _, _) = K.gcdex(e, m)
        v += e * (u * s % m)
    return v % p

def gf_crt1(M, K):
    (E, S) = ([], [])
    p = prod(M, start=K.one)
    for m in M:
        E.append(p // m)
        S.append(K.gcdex(E[-1], m)[0] % m)
    return (p, E, S)

def gf_crt2(U, M, p, E, S, K):
    v = K.zero
    for (u, m, e, s) in zip(U, M, E, S):
        v += e * (u * s % m)
    return v % p

def gf_int(a, p):
    if a <= p // 2:
        return a
    else:
        return a - p

def gf_degree(f):
    return len(f) - 1

def gf_LC(f, K):
    if not f:
        return K.zero
    else:
        return f[0]

def gf_TC(f, K):
    if not f:
        return K.zero
    else:
        return f[-1]

def gf_strip(f):
    if not f or f[0]:
        return f
    k = 0
    for coeff in f:
        if coeff:
            break
        else:
            k += 1
    return f[k:]

def gf_trunc(f, p):
    return gf_strip([a % p for a in f])

def gf_normal(f, p, K):
    return gf_trunc(list(map(K, f)), p)

def gf_from_dict(f, p, K):
    (n, h) = (max(f.keys()), [])
    if isinstance(n, SYMPY_INTS):
        for k in range(n, -1, -1):
            h.append(f.get(k, K.zero) % p)
    else:
        (n,) = n
        for k in range(n, -1, -1):
            h.append(f.get((k,), K.zero) % p)
    return gf_trunc(h, p)

def gf_to_dict(f, p, symmetric=True):
    (n, result) = (gf_degree(f), {})
    for k in range(0, n + 1):
        if symmetric:
            a = gf_int(f[n - k], p)
        else:
            a = f[n - k]
        if a:
            result[k] = a
    return result

def gf_from_int_poly(f, p):
    return gf_trunc(f, p)

def gf_to_int_poly(f, p, symmetric=True):
    if symmetric:
        return [gf_int(c, p) for c in f]
    else:
        return f

def gf_neg(f, p, K):
    return [-coeff % p for coeff in f]

def gf_add_ground(f, a, p, K):
    if not f:
        a = a % p
    else:
        a = (f[-1] + a) % p
        if len(f) > 1:
            return f[:-1] + [a]
    if not a:
        return []
    else:
        return [a]

def gf_sub_ground(f, a, p, K):
    if not f:
        a = -a % p
    else:
        a = (f[-1] - a) % p
        if len(f) > 1:
            return f[:-1] + [a]
    if not a:
        return []
    else:
        return [a]

def gf_mul_ground(f, a, p, K):
    if not a:
        return []
    else:
        return [a * b % p for b in f]

def gf_quo_ground(f, a, p, K):
    return gf_mul_ground(f, K.invert(a, p), p, K)

def gf_add(f, g, p, K):
    if not f:
        return g
    if not g:
        return f
    df = gf_degree(f)
    dg = gf_degree(g)
    if df == dg:
        return gf_strip([(a + b) % p for (a, b) in zip(f, g)])
    else:
        k = abs(df - dg)
        if df > dg:
            (h, f) = (f[:k], f[k:])
        else:
            (h, g) = (g[:k], g[k:])
        return h + [(a + b) % p for (a, b) in zip(f, g)]

def gf_sub(f, g, p, K):
    if not g:
        return f
    if not f:
        return gf_neg(g, p, K)
    df = gf_degree(f)
    dg = gf_degree(g)
    if df == dg:
        return gf_strip([(a - b) % p for (a, b) in zip(f, g)])
    else:
        k = abs(df - dg)
        if df > dg:
            (h, f) = (f[:k], f[k:])
        else:
            (h, g) = (gf_neg(g[:k], p, K), g[k:])
        return h + [(a - b) % p for (a, b) in zip(f, g)]

def gf_mul(f, g, p, K):
    df = gf_degree(f)
    dg = gf_degree(g)
    dh = df + dg
    h = [0] * (dh + 1)
    for i in range(0, dh + 1):
        coeff = K.zero
        for j in range(max(0, i - dg), min(i, df) + 1):
            coeff += f[j] * g[i - j]
        h[i] = coeff % p
    return gf_strip(h)

def gf_sqr(f, p, K):
    df = gf_degree(f)
    dh = 2 * df
    h = [0] * (dh + 1)
    for i in range(0, dh + 1):
        coeff = K.zero
        jmin = max(0, i - df)
        jmax = min(i, df)
        n = jmax - jmin + 1
        jmax = jmin + n // 2 - 1
        for j in range(jmin, jmax + 1):
            coeff += f[j] * f[i - j]
        coeff += coeff
        if n & 1:
            elem = f[jmax + 1]
            coeff += elem ** 2
        h[i] = coeff % p
    return gf_strip(h)

def gf_add_mul(f, g, h, p, K):
    return gf_add(f, gf_mul(g, h, p, K), p, K)

def gf_sub_mul(f, g, h, p, K):
    return gf_sub(f, gf_mul(g, h, p, K), p, K)

def gf_expand(F, p, K):
    if isinstance(F, tuple):
        (lc, F) = F
    else:
        lc = K.one
    g = [lc]
    for (f, k) in F:
        f = gf_pow(f, k, p, K)
        g = gf_mul(g, f, p, K)
    return g

def gf_div(f, g, p, K):
    df = gf_degree(f)
    dg = gf_degree(g)
    if not g:
        raise ZeroDivisionError('polynomial division')
    elif df < dg:
        return ([], f)
    inv = K.invert(g[0], p)
    (h, dq, dr) = (list(f), df - dg, dg - 1)
    for i in range(0, df + 1):
        coeff = h[i]
        for j in range(max(0, dg - i), min(df - i, dr) + 1):
            coeff -= h[i + j - dg] * g[dg - j]
        if i <= dq:
            coeff *= inv
        h[i] = coeff % p
    return (h[:dq + 1], gf_strip(h[dq + 1:]))

def gf_rem(f, g, p, K):
    return gf_div(f, g, p, K)[1]

def gf_quo(f, g, p, K):
    df = gf_degree(f)
    dg = gf_degree(g)
    if not g:
        raise ZeroDivisionError('polynomial division')
    elif df < dg:
        return []
    inv = K.invert(g[0], p)
    (h, dq, dr) = (f[:], df - dg, dg - 1)
    for i in range(0, dq + 1):
        coeff = h[i]
        for j in range(max(0, dg - i), min(df - i, dr) + 1):
            coeff -= h[i + j - dg] * g[dg - j]
        h[i] = coeff * inv % p
    return h[:dq + 1]

def gf_exquo(f, g, p, K):
    (q, r) = gf_div(f, g, p, K)
    if not r:
        return q
    else:
        raise ExactQuotientFailed(f, g)

def gf_lshift(f, n, K):
    if not f:
        return f
    else:
        return f + [K.zero] * n

def gf_rshift(f, n, K):
    if not n:
        return (f, [])
    else:
        return (f[:-n], f[-n:])

def gf_pow(f, n, p, K):
    if not n:
        return [K.one]
    elif n == 1:
        return f
    elif n == 2:
        return gf_sqr(f, p, K)
    h = [K.one]
    while True:
        if n & 1:
            h = gf_mul(h, f, p, K)
            n -= 1
        n >>= 1
        if not n:
            break
        f = gf_sqr(f, p, K)
    return h

def gf_frobenius_monomial_base(g, p, K):
    n = gf_degree(g)
    if n == 0:
        return []
    b = [0] * n
    b[0] = [1]
    if p < n:
        for i in range(1, n):
            mon = gf_lshift(b[i - 1], p, K)
            b[i] = gf_rem(mon, g, p, K)
    elif n > 1:
        b[1] = gf_pow_mod([K.one, K.zero], p, g, p, K)
        for i in range(2, n):
            b[i] = gf_mul(b[i - 1], b[1], p, K)
            b[i] = gf_rem(b[i], g, p, K)
    return b

def gf_frobenius_map(f, g, b, p, K):
    m = gf_degree(g)
    if gf_degree(f) >= m:
        f = gf_rem(f, g, p, K)
    if not f:
        return []
    n = gf_degree(f)
    sf = [f[-1]]
    for i in range(1, n + 1):
        v = gf_mul_ground(b[i], f[n - i], p, K)
        sf = gf_add(sf, v, p, K)
    return sf

def _gf_pow_pnm1d2(f, n, g, b, p, K):
    f = gf_rem(f, g, p, K)
    h = f
    r = f
    for i in range(1, n):
        h = gf_frobenius_map(h, g, b, p, K)
        r = gf_mul(r, h, p, K)
        r = gf_rem(r, g, p, K)
    res = gf_pow_mod(r, (p - 1) // 2, g, p, K)
    return res

def gf_pow_mod(f, n, g, p, K):
    if not n:
        return [K.one]
    elif n == 1:
        return gf_rem(f, g, p, K)
    elif n == 2:
        return gf_rem(gf_sqr(f, p, K), g, p, K)
    h = [K.one]
    while True:
        if n & 1:
            h = gf_mul(h, f, p, K)
            h = gf_rem(h, g, p, K)
            n -= 1
        n >>= 1
        if not n:
            break
        f = gf_sqr(f, p, K)
        f = gf_rem(f, g, p, K)
    return h

def gf_gcd(f, g, p, K):
    while g:
        (f, g) = (g, gf_rem(f, g, p, K))
    return gf_monic(f, p, K)[1]

def gf_lcm(f, g, p, K):
    if not f or not g:
        return []
    h = gf_quo(gf_mul(f, g, p, K), gf_gcd(f, g, p, K), p, K)
    return gf_monic(h, p, K)[1]

def gf_cofactors(f, g, p, K):
    if not f and (not g):
        return ([], [], [])
    h = gf_gcd(f, g, p, K)
    return (h, gf_quo(f, h, p, K), gf_quo(g, h, p, K))

def gf_gcdex(f, g, p, K):
    if not (f or g):
        return ([K.one], [], [])
    (p0, r0) = gf_monic(f, p, K)
    (p1, r1) = gf_monic(g, p, K)
    if not f:
        return ([], [K.invert(p1, p)], r1)
    if not g:
        return ([K.invert(p0, p)], [], r0)
    (s0, s1) = ([K.invert(p0, p)], [])
    (t0, t1) = ([], [K.invert(p1, p)])
    while True:
        (Q, R) = gf_div(r0, r1, p, K)
        if not R:
            break
        ((lc, r1), r0) = (gf_monic(R, p, K), r1)
        inv = K.invert(lc, p)
        s = gf_sub_mul(s0, s1, Q, p, K)
        t = gf_sub_mul(t0, t1, Q, p, K)
        (s1, s0) = (gf_mul_ground(s, inv, p, K), s1)
        (t1, t0) = (gf_mul_ground(t, inv, p, K), t1)
    return (s1, t1, r1)

def gf_monic(f, p, K):
    if not f:
        return (K.zero, [])
    else:
        lc = f[0]
        if K.is_one(lc):
            return (lc, list(f))
        else:
            return (lc, gf_quo_ground(f, lc, p, K))

def gf_diff(f, p, K):
    df = gf_degree(f)
    (h, n) = ([K.zero] * df, df)
    for coeff in f[:-1]:
        coeff *= K(n)
        coeff %= p
        if coeff:
            h[df - n] = coeff
        n -= 1
    return gf_strip(h)

def gf_eval(f, a, p, K):
    result = K.zero
    for c in f:
        result *= a
        result += c
        result %= p
    return result

def gf_multi_eval(f, A, p, K):
    return [gf_eval(f, a, p, K) for a in A]

def gf_compose(f, g, p, K):
    if len(g) <= 1:
        return gf_strip([gf_eval(f, gf_LC(g, K), p, K)])
    if not f:
        return []
    h = [f[0]]
    for c in f[1:]:
        h = gf_mul(h, g, p, K)
        h = gf_add_ground(h, c, p, K)
    return h

def gf_compose_mod(g, h, f, p, K):
    if not g:
        return []
    comp = [g[0]]
    for a in g[1:]:
        comp = gf_mul(comp, h, p, K)
        comp = gf_add_ground(comp, a, p, K)
        comp = gf_rem(comp, f, p, K)
    return comp

def gf_trace_map(a, b, c, n, f, p, K):
    u = gf_compose_mod(a, b, f, p, K)
    v = b
    if n & 1:
        U = gf_add(a, u, p, K)
        V = b
    else:
        U = a
        V = c
    n >>= 1
    while n:
        u = gf_add(u, gf_compose_mod(u, v, f, p, K), p, K)
        v = gf_compose_mod(v, v, f, p, K)
        if n & 1:
            U = gf_add(U, gf_compose_mod(u, V, f, p, K), p, K)
            V = gf_compose_mod(v, V, f, p, K)
        n >>= 1
    return (gf_compose_mod(a, V, f, p, K), U)

def _gf_trace_map(f, n, g, b, p, K):
    f = gf_rem(f, g, p, K)
    h = f
    r = f
    for i in range(1, n):
        h = gf_frobenius_map(h, g, b, p, K)
        r = gf_add(r, h, p, K)
        r = gf_rem(r, g, p, K)
    return r

def gf_random(n, p, K):
    pi = int(p)
    return [K.one] + [K(int(uniform(0, pi))) for i in range(0, n)]

def gf_irreducible(n, p, K):
    while True:
        f = gf_random(n, p, K)
        if gf_irreducible_p(f, p, K):
            return f

def gf_irred_p_ben_or(f, p, K):
    n = gf_degree(f)
    if n <= 1:
        return True
    (_, f) = gf_monic(f, p, K)
    if n < 5:
        H = h = gf_pow_mod([K.one, K.zero], p, f, p, K)
        for i in range(0, n // 2):
            g = gf_sub(h, [K.one, K.zero], p, K)
            if gf_gcd(f, g, p, K) == [K.one]:
                h = gf_compose_mod(h, H, f, p, K)
            else:
                return False
    else:
        b = gf_frobenius_monomial_base(f, p, K)
        H = h = gf_frobenius_map([K.one, K.zero], f, b, p, K)
        for i in range(0, n // 2):
            g = gf_sub(h, [K.one, K.zero], p, K)
            if gf_gcd(f, g, p, K) == [K.one]:
                h = gf_frobenius_map(h, f, b, p, K)
            else:
                return False
    return True

def gf_irred_p_rabin(f, p, K):
    n = gf_degree(f)
    if n <= 1:
        return True
    (_, f) = gf_monic(f, p, K)
    x = [K.one, K.zero]
    from ...ntheory import factorint
    indices = {n // d for d in factorint(n)}
    b = gf_frobenius_monomial_base(f, p, K)
    h = b[1]
    for i in range(1, n):
        if i in indices:
            g = gf_sub(h, x, p, K)
            if gf_gcd(f, g, p, K) != [K.one]:
                return False
        h = gf_frobenius_map(h, f, b, p, K)
    return h == x
_irred_methods = {'ben-or': gf_irred_p_ben_or, 'rabin': gf_irred_p_rabin}

def gf_irreducible_p(f, p, K):
    method = query('GF_IRRED_METHOD')
    if method is not None:
        irred = _irred_methods[method](f, p, K)
    else:
        irred = gf_irred_p_rabin(f, p, K)
    return irred

def gf_sqf_p(f, p, K):
    (_, f) = gf_monic(f, p, K)
    if not f:
        return True
    else:
        return gf_gcd(f, gf_diff(f, p, K), p, K) == [K.one]

def gf_sqf_part(f, p, K):
    (_, sqf) = gf_sqf_list(f, p, K)
    g = [K.one]
    for (f, _) in sqf:
        g = gf_mul(g, f, p, K)
    return g

def gf_sqf_list(f, p, K, all=False):
    (n, sqf, factors, r) = (1, False, [], int(p))
    (lc, f) = gf_monic(f, p, K)
    if gf_degree(f) < 1:
        return (lc, [])
    while True:
        F = gf_diff(f, p, K)
        if F != []:
            g = gf_gcd(f, F, p, K)
            h = gf_quo(f, g, p, K)
            i = 1
            while h != [K.one]:
                G = gf_gcd(g, h, p, K)
                H = gf_quo(h, G, p, K)
                if gf_degree(H) > 0:
                    factors.append((H, i * n))
                (g, h, i) = (gf_quo(g, G, p, K), G, i + 1)
            if g == [K.one]:
                sqf = True
            else:
                f = g
        if not sqf:
            d = gf_degree(f) // r
            for i in range(0, d + 1):
                f[i] = f[i * r]
            (f, n) = (f[:d + 1], n * r)
        else:
            break
    if all:
        raise ValueError("'all=True' is not supported yet")
    return (lc, factors)

def gf_Qmatrix(f, p, K):
    (n, r) = (gf_degree(f), int(p))
    q = [K.one] + [K.zero] * (n - 1)
    Q = [list(q)] + [[]] * (n - 1)
    for i in range(1, (n - 1) * r + 1):
        (qq, c) = ([-q[-1] * f[-1] % p], q[-1])
        for j in range(1, n):
            qq.append((q[j - 1] - c * f[-j - 1]) % p)
        if not i % r:
            Q[i // r] = list(qq)
        q = qq
    return Q

def gf_Qbasis(Q, p, K):
    (Q, n) = ([list(q) for q in Q], len(Q))
    for k in range(0, n):
        Q[k][k] = (Q[k][k] - K.one) % p
    for k in range(0, n):
        for i in range(k, n):
            if Q[k][i]:
                break
        else:
            continue
        inv = K.invert(Q[k][i], p)
        for j in range(0, n):
            Q[j][i] = Q[j][i] * inv % p
        for j in range(0, n):
            t = Q[j][k]
            Q[j][k] = Q[j][i]
            Q[j][i] = t
        for i in range(0, n):
            if i != k:
                q = Q[k][i]
                for j in range(0, n):
                    Q[j][i] = (Q[j][i] - Q[j][k] * q) % p
    for i in range(0, n):
        for j in range(0, n):
            if i == j:
                Q[i][j] = (K.one - Q[i][j]) % p
            else:
                Q[i][j] = -Q[i][j] % p
    basis = []
    for q in Q:
        if any(q):
            basis.append(q)
    return basis

def gf_berlekamp(f, p, K):
    Q = gf_Qmatrix(f, p, K)
    V = gf_Qbasis(Q, p, K)
    for (i, v) in enumerate(V):
        V[i] = gf_strip(list(reversed(v)))
    factors = [f]
    for k in range(1, len(V)):
        for f in list(factors):
            s = K.zero
            while s < p:
                g = gf_sub_ground(V[k], s, p, K)
                h = gf_gcd(f, g, p, K)
                if h != [K.one] and h != f:
                    factors.remove(f)
                    f = gf_quo(f, h, p, K)
                    factors.extend([f, h])
                if len(factors) == len(V):
                    return _sort_factors(factors, multiple=False)
                s += K.one
    return _sort_factors(factors, multiple=False)

def gf_ddf_zassenhaus(f, p, K):
    (i, g, factors) = (1, [K.one, K.zero], [])
    b = gf_frobenius_monomial_base(f, p, K)
    while 2 * i <= gf_degree(f):
        g = gf_frobenius_map(g, f, b, p, K)
        h = gf_gcd(f, gf_sub(g, [K.one, K.zero], p, K), p, K)
        if h != [K.one]:
            factors.append((h, i))
            f = gf_quo(f, h, p, K)
            g = gf_rem(g, f, p, K)
            b = gf_frobenius_monomial_base(f, p, K)
        i += 1
    if f != [K.one]:
        return factors + [(f, gf_degree(f))]
    else:
        return factors

def gf_edf_zassenhaus(f, n, p, K):
    factors = [f]
    if gf_degree(f) <= n:
        return factors
    N = gf_degree(f) // n
    if p != 2:
        b = gf_frobenius_monomial_base(f, p, K)
    t = [K.one, K.zero]
    while len(factors) < N:
        if p == 2:
            h = r = t
            for i in range(n - 1):
                r = gf_pow_mod(r, 2, f, p, K)
                h = gf_add(h, r, p, K)
            g = gf_gcd(f, h, p, K)
            t += [K.zero, K.zero]
        else:
            r = gf_random(2 * n - 1, p, K)
            h = _gf_pow_pnm1d2(r, n, f, b, p, K)
            g = gf_gcd(f, gf_sub_ground(h, K.one, p, K), p, K)
        if g != [K.one] and g != f:
            factors = gf_edf_zassenhaus(g, n, p, K) + gf_edf_zassenhaus(gf_quo(f, g, p, K), n, p, K)
    return _sort_factors(factors, multiple=False)

def gf_ddf_shoup(f, p, K):
    n = gf_degree(f)
    k = int(_ceil(_sqrt(n // 2)))
    b = gf_frobenius_monomial_base(f, p, K)
    h = gf_frobenius_map([K.one, K.zero], f, b, p, K)
    U = [[K.one, K.zero], h] + [K.zero] * (k - 1)
    for i in range(2, k + 1):
        U[i] = gf_frobenius_map(U[i - 1], f, b, p, K)
    (h, U) = (U[k], U[:k])
    V = [h] + [K.zero] * (k - 1)
    for i in range(1, k):
        V[i] = gf_compose_mod(V[i - 1], h, f, p, K)
    factors = []
    for (i, v) in enumerate(V):
        (h, j) = ([K.one], k - 1)
        for u in U:
            g = gf_sub(v, u, p, K)
            h = gf_mul(h, g, p, K)
            h = gf_rem(h, f, p, K)
        g = gf_gcd(f, h, p, K)
        f = gf_quo(f, g, p, K)
        for u in reversed(U):
            h = gf_sub(v, u, p, K)
            F = gf_gcd(g, h, p, K)
            if F != [K.one]:
                factors.append((F, k * (i + 1) - j))
            (g, j) = (gf_quo(g, F, p, K), j - 1)
    if f != [K.one]:
        factors.append((f, gf_degree(f)))
    return factors

def gf_edf_shoup(f, n, p, K):
    (N, q) = (gf_degree(f), int(p))
    if not N:
        return []
    if N <= n:
        return [f]
    (factors, x) = ([f], [K.one, K.zero])
    r = gf_random(N - 1, p, K)
    if p == 2:
        h = gf_pow_mod(x, q, f, p, K)
        H = gf_trace_map(r, h, x, n - 1, f, p, K)[1]
        h1 = gf_gcd(f, H, p, K)
        h2 = gf_quo(f, h1, p, K)
        factors = gf_edf_shoup(h1, n, p, K) + gf_edf_shoup(h2, n, p, K)
    else:
        b = gf_frobenius_monomial_base(f, p, K)
        H = _gf_trace_map(r, n, f, b, p, K)
        h = gf_pow_mod(H, (q - 1) // 2, f, p, K)
        h1 = gf_gcd(f, h, p, K)
        h2 = gf_gcd(f, gf_sub_ground(h, K.one, p, K), p, K)
        h3 = gf_quo(f, gf_mul(h1, h2, p, K), p, K)
        factors = gf_edf_shoup(h1, n, p, K) + gf_edf_shoup(h2, n, p, K) + gf_edf_shoup(h3, n, p, K)
    return _sort_factors(factors, multiple=False)

def gf_zassenhaus(f, p, K):
    factors = []
    for (factor, n) in gf_ddf_zassenhaus(f, p, K):
        factors += gf_edf_zassenhaus(factor, n, p, K)
    return _sort_factors(factors, multiple=False)

def gf_shoup(f, p, K):
    factors = []
    for (factor, n) in gf_ddf_shoup(f, p, K):
        factors += gf_edf_shoup(factor, n, p, K)
    return _sort_factors(factors, multiple=False)
_factor_methods = {'berlekamp': gf_berlekamp, 'zassenhaus': gf_zassenhaus, 'shoup': gf_shoup}

def gf_factor_sqf(f, p, K, method=None):
    (lc, f) = gf_monic(f, p, K)
    if gf_degree(f) < 1:
        return (lc, [])
    method = method or query('GF_FACTOR_METHOD')
    if method is not None:
        factors = _factor_methods[method](f, p, K)
    else:
        factors = gf_zassenhaus(f, p, K)
    return (lc, factors)

def gf_factor(f, p, K):
    (lc, f) = gf_monic(f, p, K)
    if gf_degree(f) < 1:
        return (lc, [])
    factors = []
    for (g, n) in gf_sqf_list(f, p, K)[1]:
        for h in gf_factor_sqf(g, p, K)[1]:
            factors.append((h, n))
    return (lc, _sort_factors(factors))

def gf_value(f, a):
    result = 0
    for c in f:
        result *= a
        result += c
    return result

def linear_congruence(a, b, m):
    from ... import gcdex
    if a % m == 0:
        if b % m == 0:
            return list(range(m))
        else:
            return []
    (r, _, g) = gcdex(a, m)
    if b % g != 0:
        return []
    return [(r * b // g + t * m // g) % m for t in range(g)]

def _raise_mod_power(x, s, p, f):
    from ._support import ZZ
    f_f = gf_diff(f, p, ZZ)
    alpha = gf_value(f_f, x)
    beta = -gf_value(f, x) // p ** s
    return linear_congruence(alpha, beta, p)

def _csolve_prime_las_vegas(f, p, seed=None):
    from ._support import ZZ
    from ...ntheory import sqrt_mod
    randint = _randint(seed)
    root = set()
    g = gf_pow_mod([1, 0], p - 1, f, p, ZZ)
    g = gf_sub_ground(g, 1, p, ZZ)
    factors = [gf_gcd(f, g, p, ZZ)]
    while factors:
        f = factors.pop()
        if len(f) <= 1:
            continue
        if len(f) == 2:
            root.add(-invert(f[0], p) * f[1] % p)
            continue
        if len(f) == 3:
            inv = invert(f[0], p)
            b = f[1] * inv % p
            b = (b + p * (b % 2)) // 2
            root.update(((r - b) % p for r in sqrt_mod(b ** 2 - f[2] * inv, p, all_roots=True)))
            continue
        while True:
            a = randint(0, p - 1)
            g = gf_pow_mod([1, a], (p - 1) // 2, f, p, ZZ)
            g = gf_sub_ground(g, 1, p, ZZ)
            g = gf_gcd(f, g, p, ZZ)
            if 1 < len(g) < len(f):
                factors.append(g)
                factors.append(gf_div(f, g, p, ZZ)[0])
                break
    return sorted(root)

def csolve_prime(f, p, e=1):
    from ._support import ZZ
    g = [MPZ(int(c)) for c in f]
    for i in range(len(g) - p):
        g[i + p - 1] += g[i]
        g[i] = 0
    g = gf_trunc(g, p)
    k = 0
    while k < len(g) and g[len(g) - k - 1] == 0:
        k += 1
    if k:
        g = g[:-k]
        root_zero = [0]
    else:
        root_zero = []
    if g == []:
        X1 = list(range(p))
    elif len(g) ** 2 < p:
        X1 = root_zero + _csolve_prime_las_vegas(g, p)
    else:
        X1 = root_zero + [i for i in range(p) if gf_eval(g, i, p, ZZ) == 0]
    if e == 1:
        return X1
    X = []
    S = list(zip(X1, [1] * len(X1)))
    while S:
        (x, s) = S.pop()
        if s == e:
            X.append(x)
        else:
            s1 = s + 1
            ps = p ** s
            S.extend([(x + v * ps, s1) for v in _raise_mod_power(x, s, p, f)])
    return sorted(X)

def gf_csolve(f, n):
    from ._support import ZZ
    from ...ntheory import factorint
    P = factorint(n)
    X = [csolve_prime(f, p, e) for (p, e) in P.items()]
    pools = list(map(tuple, X))
    perms = [[]]
    for pool in pools:
        perms = [x + [y] for x in perms for y in pool]
    dist_factors = [pow(p, e) for (p, e) in P.items()]
    return sorted([gf_crt(per, dist_factors, ZZ) for per in perms])
