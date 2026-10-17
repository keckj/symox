"""Chinese remaindering (``sympy.ntheory.modular``)."""

from math import gcd, prod

from ._intfunc import as_int, igcdex
from .primetest import isprime


def symmetric_residue(a, m):
    """The residue of a mod m in ``(-m/2, m/2]``."""
    if a <= m // 2:
        return a
    return a - m


def _gf_crt(U, M):
    p = prod(M)
    v = 0
    for u, m in zip(U, M):
        e = p // m
        s, _, _ = igcdex(e, m)
        v += e * (u * s % m)
    return v % p


def crt(m, v, symmetric=False, check=True):
    """Chinese remainder theorem: ``(x, lcm)`` with ``x = v[i] mod m[i]``."""
    if check:
        m = list(map(as_int, m))
        v = list(map(as_int, v))
    result = _gf_crt(v, m)
    mm = prod(m)
    if check:
        if not all(vi % mi == result % mi for vi, mi in zip(v, m)):
            result = solve_congruence(*list(zip(v, m)),
                                      check=False, symmetric=symmetric)
            if result is None:
                return result
            result, mm = result
    if symmetric:
        return int(symmetric_residue(result, mm)), int(mm)
    return int(result), int(mm)


def crt1(m):
    """First part of the CRT for repeated use with the same moduli."""
    p = prod(m)
    E, S = [], []
    for mi in m:
        e = p // mi
        s, _, _ = igcdex(e, mi)
        E.append(e)
        S.append(s % mi)
    return p, E, S


def crt2(m, v, mm, e, s, symmetric=False):
    """Second part of the CRT, using the output of ``crt1``."""
    result = 0
    for u, mi, ei, si in zip(v, m, e, s):
        result += ei * (u * si % mi)
    result %= mm
    if symmetric:
        return int(symmetric_residue(result, mm)), int(mm)
    return int(result), int(mm)


def solve_congruence(*remainder_modulus_pairs, **hint):
    """Solve a system of congruences ``x = r mod m``; ``None`` if none."""

    def combine(c1, c2):
        a1, m1 = c1
        a2, m2 = c2
        a, b, c = m1, a2 - a1, m2
        g = gcd(a, b, c)
        a, b, c = [i // g for i in [a, b, c]]
        if a != 1:
            inv_a, _, g = igcdex(a, c)
            if g != 1:
                return None
            b *= inv_a
        a, m = a1 + m1 * b, m1 * c
        return a, m

    rm = remainder_modulus_pairs
    symmetric = hint.get('symmetric', False)
    if hint.get('check', True):
        rm = [(as_int(r), as_int(m)) for r, m in rm]
        uniq = {}
        for r, m in rm:
            r %= m
            if m in uniq:
                if r != uniq[m]:
                    return None
                continue
            uniq[m] = r
        rm = [(r, m) for m, r in uniq.items()]
        del uniq
        if all(isprime(m) for r, m in rm):
            r, m = list(zip(*rm))
            return crt(m, r, symmetric=symmetric, check=False)
    rv = (0, 1)
    for rmi in rm:
        rv = combine(rv, rmi)
        if rv is None:
            break
        n, m = rv
        n = n % m
    else:
        if symmetric:
            return symmetric_residue(n, m), m
        return n, m
