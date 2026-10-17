"""Primality tests (``sympy.ntheory.primetest``)."""

import math

from .. import _core
from ._intfunc import as_int

MERSENNE_PRIME_EXPONENTS = tuple(_core.nt_mersenne_exponents())


def isprime(n):
    """Deterministic for n < 2**64, BPSW strong probable prime above."""
    n = as_int(n)
    return _core.isprime(n)


def is_mersenne_prime(n):
    """Whether n is a Mersenne prime ``2**p - 1``."""
    n = as_int(n)
    if n < 1:
        return False
    return _core.nt_is_mersenne_prime(n)


def is_square(n, prep=True):
    """Whether n is a perfect square."""
    if prep:
        n = as_int(n)
        if n < 0:
            return False
        if n in (0, 1):
            return True
    return math.isqrt(n) ** 2 == n


def is_gaussian_prime(num):
    """Whether ``num`` is a Gaussian prime."""
    from .. import sympify, re, im
    num = sympify(num)
    a, b = re(num), im(num)
    a = as_int(a, strict=False)
    b = as_int(b, strict=False)
    if a == 0:
        b = abs(b)
        return isprime(b) and b % 4 == 3
    elif b == 0:
        a = abs(a)
        return isprime(a) and a % 4 == 3
    return isprime(a ** 2 + b ** 2)


def is_fermat_pseudoprime(n, a):
    """Whether n is a Fermat probable prime to base a."""
    return _core.nt_prp("fermat", as_int(n), as_int(a))


def is_euler_pseudoprime(n, a):
    """Whether n is an Euler probable prime to base a."""
    return _core.nt_prp("euler", as_int(n), as_int(a))


def is_euler_jacobi_pseudoprime(n, a):
    """Whether n is an Euler-Jacobi probable prime to base a."""
    return _core.nt_prp("euler_jacobi", as_int(n), as_int(a))


def mr(n, bases):
    """Miller-Rabin strong pseudoprime test for the given bases."""
    return _core.nt_mr(as_int(n), [as_int(b) for b in bases])


def is_lucas_prp(n):
    """Standard Lucas probable prime test (Selfridge parameters)."""
    return _core.nt_lucas_prp("lucas", as_int(n))


def is_strong_lucas_prp(n):
    """Strong Lucas probable prime test (Selfridge parameters)."""
    return _core.nt_lucas_prp("strong", as_int(n))


def is_extra_strong_lucas_prp(n):
    """Extra strong Lucas probable prime test."""
    return _core.nt_lucas_prp("extra_strong", as_int(n))


def proth_test(n):
    """Primality test for Proth numbers ``k*2**m + 1`` with ``k < 2**m``."""
    return _core.nt_proth_test(as_int(n))
