"""Discrete transforms (``sympy.discrete.transforms``).

Radix-2 FFT, number theoretic transform, Walsh-Hadamard and Moebius
transforms, following SymPy's algorithms on symox expressions (or on
Python integers for the NTT).
"""

from .. import _core
from ..ntheory._intfunc import as_int


def iterable(i):
    """``sympy.utilities.iterables.iterable`` with its default exclusions."""
    if isinstance(i, (str, dict)):
        return False
    try:
        iter(i)
    except TypeError:
        return False
    return True


def _ints(seq):
    out = []
    for x in seq:
        if getattr(x, 'is_Integer', False) is True:
            out.append(int(x))
        else:
            return None
    return out


def _bitrev(i, b):
    return int(format(i, '0%db' % b)[::-1], 2) if b else 0


def _fourier_transform(seq, dps, inverse=False):
    from .. import S, Symbol, sympify, expand_mul, pi, I, sin, cos
    if not iterable(seq):
        raise TypeError("Expected a sequence of numeric coefficients "
                        "for Fourier Transform")
    a = [sympify(arg) for arg in seq]
    if any(x.atoms(Symbol) for x in a):
        raise ValueError("Expected non-symbolic coefficients")
    n = len(a)
    if n < 2:
        return a
    b = n.bit_length() - 1
    if n & (n - 1):
        b += 1
        n = 2 ** b
    a += [S.Zero] * (n - len(a))
    for i in range(1, n):
        j = _bitrev(i, b)
        if i < j:
            a[i], a[j] = a[j], a[i]
    ang = -2 * pi / n if inverse else 2 * pi / n
    if dps is not None:
        ang = ang.evalf(dps + 2)
    w = [cos(ang * i) + I * sin(ang * i) for i in range(n // 2)]
    h = 2
    while h <= n:
        hf, ut = h // 2, n // h
        for i in range(0, n, h):
            for j in range(hf):
                u, v = a[i + j], expand_mul(a[i + j + hf] * w[ut * j])
                a[i + j], a[i + j + hf] = u + v, u - v
        h *= 2
    if inverse:
        a = [(x / n).evalf(dps) for x in a] if dps is not None \
            else [x / n for x in a]
    return a


def fft(seq, dps=None):
    """Discrete Fourier transform (radix-2, zero padded)."""
    return _fourier_transform(seq, dps=dps)


def ifft(seq, dps=None):
    """Inverse discrete Fourier transform (radix-2, zero padded)."""
    return _fourier_transform(seq, dps=dps, inverse=True)


def _number_theoretic_transform(seq, prime, inverse=False):
    from ..ntheory import isprime
    if not iterable(seq):
        raise TypeError("Expected a sequence of integer coefficients "
                        "for Number Theoretic Transform")
    p = as_int(prime)
    if not isprime(p):
        raise ValueError("Expected prime modulus for "
                         "Number Theoretic Transform")
    a = [as_int(x) for x in seq]
    return [int(v) for v in _core.nt_ntt(a, p, bool(inverse))]


def ntt(seq, prime):
    """Number theoretic transform modulo a prime ``m*2**k + 1``."""
    return _number_theoretic_transform(seq, prime=prime)


def intt(seq, prime):
    """Inverse number theoretic transform."""
    return _number_theoretic_transform(seq, prime=prime, inverse=True)


def _walsh_hadamard_transform(seq, inverse=False):
    from .. import S, sympify
    if not iterable(seq):
        raise TypeError("Expected a sequence of coefficients "
                        "for Walsh Hadamard Transform")
    a = [sympify(arg) for arg in seq]
    n = len(a)
    if n < 2:
        return a
    ia = _ints(a)
    if ia is not None:
        from .. import Integer, Rational
        r = _core.nt_fwht(ia)
        if inverse:
            return [Rational(x, len(r)) for x in r]
        return [Integer(x) for x in r]
    if n & (n - 1):
        n = 2 ** n.bit_length()
    a += [S.Zero] * (n - len(a))
    h = 2
    while h <= n:
        hf = h // 2
        for i in range(0, n, h):
            for j in range(hf):
                u, v = a[i + j], a[i + j + hf]
                a[i + j], a[i + j + hf] = u + v, u - v
        h *= 2
    if inverse:
        a = [x / n for x in a]
    return a


def fwht(seq):
    """Walsh-Hadamard transform (Hadamard ordering)."""
    return _walsh_hadamard_transform(seq)


def ifwht(seq):
    """Inverse Walsh-Hadamard transform."""
    return _walsh_hadamard_transform(seq, inverse=True)


def _mobius_transform(seq, sgn, subset):
    from .. import S, sympify
    if not iterable(seq):
        raise TypeError("Expected a sequence of coefficients")
    a = [sympify(arg) for arg in seq]
    n = len(a)
    if n < 2:
        return a
    ia = _ints(a)
    if ia is not None:
        from .. import Integer
        return [Integer(x) for x in
                _core.nt_mobius_transform(ia, sgn, bool(subset))]
    if n & (n - 1):
        n = 2 ** n.bit_length()
    a += [S.Zero] * (n - len(a))
    i = 1
    while i < n:
        for j in range(n):
            if subset:
                if j & i:
                    a[j] += sgn * a[j ^ i]
            elif not j & i:
                a[j] += sgn * a[j ^ i]
        i *= 2
    return a


def mobius_transform(seq, subset=True):
    """Zeta (Moebius) transform over subset or superset lattices."""
    return _mobius_transform(seq, sgn=+1, subset=subset)


def inverse_mobius_transform(seq, subset=True):
    """Inverse zeta (Moebius) transform."""
    return _mobius_transform(seq, sgn=-1, subset=subset)
