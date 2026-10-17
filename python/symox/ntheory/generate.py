"""Prime generation (``sympy.ntheory.generate``).

The ``Sieve`` keeps SymPy's public state (``array`` lists of primes,
totients and Moebius values); sieving and primality run natively.
"""

from array import array as _array
from bisect import bisect, bisect_left
from itertools import count
import random as _random

from .. import _core
from ._intfunc import as_int, _int_ceiling

_NATIVE_SIEVE_LIMIT = 10 ** 12


def _as_int_ceiling(a):
    return as_int(_int_ceiling(a))


def _isprime(n):
    return _core.isprime(n)


def _native_primerange(a, b):
    """Primes in ``[a, b)`` for Python ints."""
    a = max(a, 2)
    if a >= b:
        return []
    if b <= _NATIVE_SIEVE_LIMIT:
        return _core.nt_primerange(a, b)
    out = []
    p = _core.nt_nextprime(a - 1, 1)
    while p < b:
        out.append(p)
        p = _core.nt_nextprime(p, 1)
    return out


class Sieve:
    """A list of prime numbers, implemented as a dynamically growing sieve."""

    def __init__(self, sieve_interval=1000000):
        self._n = 6
        self._list = _array('L', [2, 3, 5, 7, 11, 13])
        self._tlist = _array('L', [0, 1, 1, 2, 2, 4])
        self._mlist = _array('i', [0, 1, -1, -1, 0, -1])
        if sieve_interval <= 0:
            raise ValueError('sieve_interval should be a positive integer')
        self.sieve_interval = sieve_interval

    def __repr__(self):
        return ("<%s sieve (%i): %i, %i, %i, ... %i, %i\n"
                "%s sieve (%i): %i, %i, %i, ... %i, %i\n"
                "%s sieve (%i): %i, %i, %i, ... %i, %i>") % (
            'prime', len(self._list), self._list[0], self._list[1], self._list[2],
            self._list[-2], self._list[-1],
            'totient', len(self._tlist), self._tlist[0], self._tlist[1], self._tlist[2],
            self._tlist[-2], self._tlist[-1],
            'mobius', len(self._mlist), self._mlist[0], self._mlist[1], self._mlist[2],
            self._mlist[-2], self._mlist[-1])

    def _reset(self, prime=None, totient=None, mobius=None):
        if all(i is None for i in (prime, totient, mobius)):
            prime = totient = mobius = True
        if prime:
            self._list = self._list[:self._n]
        if totient:
            self._tlist = self._tlist[:self._n]
        if mobius:
            self._mlist = self._mlist[:self._n]

    def extend(self, n):
        """Grow the sieve to cover all primes <= n."""
        n = int(n)
        num = self._list[-1] + 1
        if n < num:
            return
        self._list += _array('L', _native_primerange(num, n + 1))

    def _primerange(self, a, b):
        yield from _native_primerange(a, b)

    def extend_to_no(self, i):
        """Extend to include the ith prime number."""
        i = as_int(i)
        while len(self._list) < i:
            self.extend(int(self._list[-1] * 1.5))

    def primerange(self, a, b=None):
        """Generate all prime numbers in the range [2, a) or [a, b)."""
        if b is None:
            b = _as_int_ceiling(a)
            a = 2
        else:
            a = max(2, _as_int_ceiling(a))
            b = _as_int_ceiling(b)
        if a >= b:
            return
        self.extend(b)
        yield from self._list[bisect_left(self._list, a):bisect_left(self._list, b)]

    def totientrange(self, a, b):
        """Generate all totient numbers for the range [a, b)."""
        a = max(1, _as_int_ceiling(a))
        b = _as_int_ceiling(b)
        n = len(self._tlist)
        if a >= b:
            return
        if b > n:
            self._tlist += _array('L', _core.nt_totient_range(n, b))
        for i in range(a, b):
            yield self._tlist[i]

    def mobiusrange(self, a, b):
        """Generate all mobius numbers for the range [a, b)."""
        a = max(1, _as_int_ceiling(a))
        b = _as_int_ceiling(b)
        n = len(self._mlist)
        if a >= b:
            return
        if b > n:
            self._mlist += _array('i', _core.nt_mobius_range(n, b))
        for i in range(a, b):
            yield self._mlist[i]

    def search(self, n):
        """Return the indices i, j of the primes that bound n."""
        test = _as_int_ceiling(n)
        n = as_int(n)
        if n < 2:
            raise ValueError('n should be >= 2 but got: %s' % n)
        if n > self._list[-1]:
            self.extend(n)
        b = bisect(self._list, n)
        if self._list[b - 1] == test:
            return b, b
        return b, b + 1

    def __contains__(self, n):
        try:
            n = as_int(n)
            assert n >= 2
        except (ValueError, AssertionError):
            return False
        if n % 2 == 0:
            return n == 2
        a, b = self.search(n)
        return a == b

    def __iter__(self):
        for n in count(1):
            yield self[n]

    def __getitem__(self, n):
        if isinstance(n, slice):
            self.extend_to_no(n.stop)
            start = n.start if n.start is not None else 0
            if start < 1:
                raise IndexError('Sieve indices start at 1.')
            return self._list[start - 1:n.stop - 1:n.step]
        if n < 1:
            raise IndexError('Sieve indices start at 1.')
        n = as_int(n)
        self.extend_to_no(n)
        return self._list[n - 1]


sieve = Sieve()


def prime(nth):
    """The nth prime (``prime(1) == 2``)."""
    n = as_int(nth)
    if n < 1:
        raise ValueError('nth must be a positive integer; prime(1) == 2')
    return int(_core.nt_prime(n))


def _primepi(n):
    if n < 2:
        return 0
    return int(_core.primepi(n))


def nextprime(n, ith=1):
    """The ith prime greater than n."""
    n = int(n)
    i = as_int(ith)
    if i <= 0:
        raise ValueError('ith should be positive')
    return int(_core.nt_nextprime(n, i))


def prevprime(n):
    """The largest prime smaller than n."""
    n = _as_int_ceiling(n)
    if n < 3:
        raise ValueError('no preceding primes')
    return int(_core.nt_prevprime(n))


def primerange(a, b=None):
    """Generate the primes in [2, a) or [a, b)."""
    if b is None:
        a, b = 2, a
    if a >= b:
        return
    lo = max(2, _int_ceiling(a))
    hi = _int_ceiling(b)
    yield from _native_primerange(lo, hi)


def randprime(a, b):
    """A random prime in [a, b) (uses Python's ``random`` module)."""
    if a >= b:
        return
    a, b = map(int, (a, b))
    n = _random.randint(a - 1, b)
    p = nextprime(n)
    if p >= b:
        p = prevprime(b)
    if p < a:
        raise ValueError('no primes exist in the specified range')
    return p


def primorial(n, nth=True):
    """Product of the first n primes (or of the primes <= n)."""
    if nth:
        n = as_int(n)
    else:
        n = int(n)
    if n < 1:
        raise ValueError('primorial argument must be >= 1')
    return int(_core.nt_primorial(n, bool(nth)))


def cycle_length(f, x0, nmax=None, values=False):
    """Brent's cycle detection: yields ``(lambda, mu)`` or the values."""
    nmax = int(nmax or 0)
    power = lam = 1
    tortoise, hare = x0, f(x0)
    i = 1
    if values:
        yield tortoise
    while tortoise != hare and (not nmax or i < nmax):
        i += 1
        if power == lam:
            tortoise = hare
            power *= 2
            lam = 0
        if values:
            yield hare
        hare = f(hare)
        lam += 1
    if nmax and i == nmax:
        if values:
            return
        yield nmax, None
        return
    if not values:
        mu = 0
        tortoise = hare = x0
        for i in range(lam):
            hare = f(hare)
        while tortoise != hare:
            tortoise = f(tortoise)
            hare = f(hare)
            mu += 1
        yield lam, mu


def composite(nth):
    """The nth composite number (``composite(1) == 4``)."""
    n = as_int(nth)
    if n < 1:
        raise ValueError('nth must be a positive integer; composite(1) == 4')
    return int(_core.nt_composite(n))


def compositepi(n):
    """Number of composite numbers <= n."""
    n = int(n)
    if n < 4:
        return 0
    return int(_core.nt_compositepi(n))


def primepi(n):
    """Deprecated in SymPy: the ``primepi`` function head."""
    from .. import primepi as _primepi_head
    return _primepi_head(n)
