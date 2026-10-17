"""Convolutions (``sympy.discrete.convolutions``)."""

from math import lcm

from ..ntheory._intfunc import as_int
from .transforms import (fft, ifft, ntt, intt, fwht, ifwht, mobius_transform,
                         inverse_mobius_transform, iterable)


def convolution(a, b, cycle=0, dps=None, prime=None, dyadic=None,
                subset=None):
    """Linear, cyclic, modular, dyadic or subset convolution of a and b."""
    from .. import Rational
    c = as_int(cycle)
    if c < 0:
        raise ValueError("The length for cyclic convolution "
                         "must be non-negative")
    dyadic = True if dyadic else None
    subset = True if subset else None
    if sum(x is not None for x in (prime, dps, dyadic, subset)) > 1:
        raise TypeError("Ambiguity in determining the type of convolution")
    if prime is not None:
        ls = convolution_ntt(a, b, prime=prime)
        return ls if not c else [sum(ls[i::c]) % prime for i in range(c)]
    if dyadic:
        ls = convolution_fwht(a, b)
    elif subset:
        ls = convolution_subset(a, b)
    else:
        def loop(a):
            dens = []
            for i in a:
                if isinstance(i, Rational) and i.q - 1:
                    dens.append(i.q)
                elif not isinstance(i, int):
                    return
            if dens:
                l = lcm(*dens)
                return [i * l if type(i) is int else i.p * (l // i.q)
                        for i in a], l
            return a, 1
        ls = None
        da = loop(a)
        if da is not None:
            db = loop(b)
            if db is not None:
                (ia, ma), (ib, mb) = da, db
                den = ma * mb
                ls = convolution_int(ia, ib)
                if den != 1:
                    ls = [Rational(i, den) for i in ls]
        if ls is None:
            ls = convolution_fft(a, b, dps)
    return ls if not c else [sum(ls[i::c]) for i in range(c)]


def convolution_fft(a, b, dps=None):
    """Linear convolution through the radix-2 FFT."""
    from .. import S, expand_mul
    a, b = a[:], b[:]
    n = m = len(a) + len(b) - 1
    if n > 0 and n & (n - 1):
        n = 2 ** n.bit_length()
    a += [S.Zero] * (n - len(a))
    b += [S.Zero] * (n - len(b))
    a, b = fft(a, dps), fft(b, dps)
    a = [expand_mul(x * y) for x, y in zip(a, b)]
    a = ifft(a, dps)[:m]
    return a


def convolution_ntt(a, b, prime):
    """Linear convolution modulo a prime through the NTT."""
    a, b, p = a[:], b[:], as_int(prime)
    n = m = len(a) + len(b) - 1
    if n > 0 and n & (n - 1):
        n = 2 ** n.bit_length()
    a += [0] * (n - len(a))
    b += [0] * (n - len(b))
    a, b = ntt(a, p), ntt(b, p)
    a = [x * y % p for x, y in zip(a, b)]
    a = intt(a, p)[:m]
    return a


def convolution_fwht(a, b):
    """Dyadic (XOR) convolution through the Walsh-Hadamard transform."""
    from .. import S, expand_mul
    if not a or not b:
        return []
    a, b = a[:], b[:]
    n = max(len(a), len(b))
    if n & (n - 1):
        n = 2 ** n.bit_length()
    a += [S.Zero] * (n - len(a))
    b += [S.Zero] * (n - len(b))
    a, b = fwht(a), fwht(b)
    a = [expand_mul(x * y) for x, y in zip(a, b)]
    a = ifwht(a)
    return a


def convolution_subset(a, b):
    """Subset convolution."""
    from .. import S, sympify, expand_mul
    if not a or not b:
        return []
    if not iterable(a) or not iterable(b):
        raise TypeError("Expected a sequence of coefficients for convolution")
    a = [sympify(arg) for arg in a]
    b = [sympify(arg) for arg in b]
    n = max(len(a), len(b))
    if n & (n - 1):
        n = 2 ** n.bit_length()
    a += [S.Zero] * (n - len(a))
    b += [S.Zero] * (n - len(b))
    c = [S.Zero] * n
    for mask in range(n):
        smask = mask
        while smask > 0:
            c[mask] += expand_mul(a[smask] * b[mask ^ smask])
            smask = (smask - 1) & mask
        c[mask] += expand_mul(a[smask] * b[mask ^ smask])
    return c


def covering_product(a, b):
    """Covering (OR) product through the subset Moebius transform."""
    from .. import S, expand_mul
    if not a or not b:
        return []
    a, b = a[:], b[:]
    n = max(len(a), len(b))
    if n & (n - 1):
        n = 2 ** n.bit_length()
    a += [S.Zero] * (n - len(a))
    b += [S.Zero] * (n - len(b))
    a, b = mobius_transform(a), mobius_transform(b)
    a = [expand_mul(x * y) for x, y in zip(a, b)]
    a = inverse_mobius_transform(a)
    return a


def intersecting_product(a, b):
    """Intersecting (AND) product through the superset Moebius transform."""
    from .. import S, expand_mul
    if not a or not b:
        return []
    a, b = a[:], b[:]
    n = max(len(a), len(b))
    if n & (n - 1):
        n = 2 ** n.bit_length()
    a += [S.Zero] * (n - len(a))
    b += [S.Zero] * (n - len(b))
    a, b = (mobius_transform(a, subset=False),
            mobius_transform(b, subset=False))
    a = [expand_mul(x * y) for x, y in zip(a, b)]
    a = inverse_mobius_transform(a, subset=False)
    return a


def convolution_int(a, b):
    """Exact integer convolution by Kronecker substitution."""
    B = max(abs(c) for c in a) * max(abs(c) for c in b) * \
        (1 + min(len(a) - 1, len(b) - 1))
    x, power = 1, 0
    while x <= (2 * B):
        x <<= 1
        power += 1

    def to_integer(poly):
        n, mul = 0, 0
        for c in reversed(poly):
            if c and not mul:
                mul = -1 if c < 0 else 1
            n <<= power
            n += mul * int(c)
        return mul, n

    (a_mul, a_packed), (b_mul, b_packed) = to_integer(a), to_integer(b)
    result = a_packed * b_packed
    mul = a_mul * b_mul
    mask, half, borrow, poly = x - 1, x >> 1, 0, []
    while result or borrow:
        coeff = (result & mask) + borrow
        result >>= power
        borrow = coeff >= half
        poly.append(mul * int(coeff if coeff < half else coeff - x))
    return poly or [0]
