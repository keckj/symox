"""Integer factorization (``sympy.ntheory.factor_``).

Factoring, divisor enumeration, Pollard's methods and perfect power
detection run natively; this module reproduces SymPy's argument handling
and result types (plain ``int`` values, dicts and lists).
"""

from collections import OrderedDict, defaultdict
from collections.abc import MutableMapping
import math

from .. import _core
from ._intfunc import as_int
from .primetest import isprime, is_mersenne_prime, MERSENNE_PRIME_EXPONENTS
from .generate import sieve, primerange, nextprime


def _sym():
    from .. import S, Mul, Pow, Rational, Integer
    return S, Mul, Pow, Rational, Integer


def _is_mul(n):
    return getattr(n, 'is_Mul', False) is True


def _remove(n, p):
    """``(n / p**m, m)`` with ``m`` maximal."""
    m = 0
    while n % p == 0:
        n //= p
        m += 1
    return n, m


def smoothness(n):
    """``(max prime factor, max prime power factor)`` of n."""
    if n == 1:
        return (1, 1)
    facs = factorint(n)
    return max(facs), max(m ** facs[m] for m in facs)


def smoothness_p(n, m=-1, power=0, visual=None):
    """Smoothness of ``p + m`` for the prime factors p of n."""
    if visual in (1, 0):
        visual = bool(visual)
    elif visual not in (True, False):
        visual = None
    if isinstance(n, str):
        if visual:
            return n
        d = {}
        for li in n.splitlines():
            k, v = [int(i) for i in
                    li.split('has')[0].split('=')[1].split('**')]
            d[k] = v
        if visual is not True and visual is not False:
            return d
        return smoothness_p(d, visual=False)
    elif not isinstance(n, tuple):
        facs = factorint(n, visual=False)
    k = -1 if power else 1
    if isinstance(n, tuple):
        rv = n
    else:
        rv = (m, sorted([(f, tuple([M] + list(smoothness(f + m))))
                         for f, M in list(facs.items())],
                        key=lambda x: (x[1][k], x[0])))
    if visual is False or (visual is not True) and (
            type(n) is int or _is_mul(n)):
        return rv
    lines = []
    for f, rest in rv[1]:
        dat = [f] + list(rest)
        dat.insert(2, m)
        lines.append('p**i=%i**%i has p%+i B=%i, B-pow=%i' % tuple(dat))
    return '\n'.join(lines)


def multiplicity(p, n):
    """The greatest integer m such that p**m divides n."""
    try:
        p, n = as_int(p), as_int(n)
    except ValueError:
        S, Mul, Pow, Rational, Integer = _sym()
        from .. import factorial

        def _is_rat(i):
            return isinstance(i, int) or getattr(i, 'is_Rational', False) is True
        if all(_is_rat(i) for i in (p, n)):
            p = Rational(p)
            n = Rational(n)
            if p.q == 1:
                if n.p == 1:
                    return -multiplicity(p.p, n.q)
                return multiplicity(p.p, n.p) - multiplicity(p.p, n.q)
            elif p.p == 1:
                return multiplicity(p.q, n.q)
            else:
                like = min(multiplicity(p.p, n.p), multiplicity(p.q, n.q))
                cross = min(multiplicity(p.q, n.p), multiplicity(p.p, n.q))
                return like - cross
        elif ((isinstance(p, int) or getattr(p, 'is_Integer', False) is True)
              and getattr(n, 'func', None) == factorial
              and getattr(n.args[0], 'is_Integer', False) is True and n.args[0] >= 0):
            return multiplicity_in_factorial(p, n.args[0])
        raise ValueError('expecting ints or fractions, got %s and %s' % (p, n))
    if n == 0:
        raise ValueError('no such integer exists: multiplicity of %s is not-defined' % n)
    return int(_core.nt_multiplicity(p, n))


def multiplicity_in_factorial(p, n):
    """The largest m such that p**m divides n!."""
    p, n = as_int(p), as_int(n)
    if p <= 0:
        raise ValueError('expecting positive integer got %s' % p)
    if n < 0:
        raise ValueError('expecting non-negative integer got %s' % n)
    return int(_core.nt_multiplicity_in_factorial(p, n))


def perfect_power(n, candidates=None, big=True, factor=True):
    """``(b, e)`` with ``b**e == n`` and ``e > 1``, else ``False``."""
    if getattr(n, 'is_Rational', False) is True and getattr(n, 'is_Integer', False) is not True:
        S, Mul, Pow, Rational, Integer = _sym()
        p, q = n.p, n.q
        if p == 1:
            qq = perfect_power(q, candidates, big, factor)
            return (S.One / qq[0], qq[1]) if qq is not False else False
        pp = perfect_power(p, factor=factor)
        if not pp:
            return False
        qq = perfect_power(q, factor=factor)
        if not qq:
            return False
        (num_base, num_exp), (den_base, den_exp) = pp, qq

        def compute_tuple(exponent):
            new_num = num_base ** (num_exp // exponent)
            new_den = den_base ** (den_exp // exponent)
            return Rational(new_num, new_den), exponent
        if candidates:
            valid = [i for i in candidates if num_exp % i == 0 and den_exp % i == 0]
            if not valid:
                return False
            e = max(valid) if big else min(valid)
            return compute_tuple(e)
        g = math.gcd(num_exp, den_exp)
        if g == 1:
            return False
        if big:
            return compute_tuple(g)
        e = next(p for p in primerange(2, g + 1) if g % p == 0)
        return compute_tuple(e)
    n = as_int(n)
    if candidates is not None:
        candidates = sorted({as_int(c) for c in candidates})
    r = _core.nt_perfect_power(n, candidates, bool(big))
    if r is None:
        return False
    return int(r[0]), int(r[1])


class FactorCache(MutableMapping):
    """Cache of a prime factor of composite numbers (SymPy compatible)."""

    def __init__(self, maxsize=None):
        self._cache = OrderedDict()
        self.maxsize = maxsize

    def __len__(self):
        return len(self._cache)

    def __contains__(self, n):
        return n in self._cache

    def __getitem__(self, n):
        factor = self.get(n)
        if factor is None:
            raise KeyError(f'{n} does not exist.')
        return factor

    def __setitem__(self, n, factor):
        if not (1 < factor <= n and n % factor == 0 and isprime(factor)):
            raise ValueError(f'{factor} is not a prime factor of {n}')
        self._cache[n] = max(self._cache.get(n, 0), factor)
        if self.maxsize is not None and len(self._cache) > self.maxsize:
            self._cache.popitem(False)

    def __delitem__(self, n):
        if n not in self._cache:
            raise KeyError(f'{n} does not exist.')
        del self._cache[n]

    def __iter__(self):
        return self._cache.__iter__()

    def cache_clear(self):
        self._cache = OrderedDict()

    @property
    def maxsize(self):
        return self._maxsize

    @maxsize.setter
    def maxsize(self, value):
        if value is not None and value <= 0:
            raise ValueError('maxsize must be None or a non-negative integer.')
        self._maxsize = value
        if value is not None:
            while len(self._cache) > value:
                self._cache.popitem(False)

    def get(self, n, default=None):
        if n <= sieve._list[-1]:
            from bisect import bisect_left
            if sieve._list[bisect_left(sieve._list, n)] == n:
                return n
        if n in self._cache:
            self._cache.move_to_end(n)
            return self._cache[n]
        if (factors := self.get_external(n)):
            self.add(n, factors)
            return self._cache[n]
        return default

    def add(self, n, factors):
        for p in sorted(factors, reverse=True):
            self[n] = p
            n, _ = _remove(n, p)

    def get_external(self, n):
        return None


factor_cache = FactorCache(maxsize=1000)


def pollard_rho(n, s=2, a=1, retries=5, seed=1234, max_steps=None, F=None):
    """Pollard's rho method; a nontrivial factor of n or ``None``."""
    n = int(n)
    if n < 5:
        raise ValueError('pollard_rho should receive n > 4')
    r = _core.nt_pollard_rho(n, int(s), int(a), int(retries), int(seed),
                             None if max_steps is None else int(max_steps), F or None)
    return None if r is None else int(r)


def pollard_pm1(n, B=10, a=2, retries=0, seed=1234):
    """Pollard's p-1 method; a nontrivial factor of n or ``None``."""
    n = int(n)
    if n < 4 or B < 3:
        raise ValueError('pollard_pm1 should receive n > 3 and B > 2')
    r = _core.nt_pollard_pm1(n, int(B), int(a), int(retries), int(seed))
    return None if r is None else int(r)


def _factorint_small(factors, n, limit, fail_max):
    """Trial division port of SymPy's ``_factorint_small`` (next_p == 2)."""
    def done(n, d):
        if d * d <= n:
            return n, d
        return n, 0

    def take(n, p):
        n, m = _remove(n, p)
        factors[p] = m
        return n
    limit2 = limit ** 2
    threshold2 = min(n, limit2)
    if not n & 1:
        n = take(n, 2)
        threshold2 = min(n, limit2)
    next_p = 3
    if threshold2 < 9:
        return done(n, next_p)
    if not n % 3:
        n = take(n, 3)
        threshold2 = min(n, limit2)
    next_p = 5
    if threshold2 < 25:
        return done(n, next_p)
    fails = 0
    while fails < fail_max:
        for step in (2, 4):
            if n % next_p:
                fails += 1
            else:
                n = take(n, next_p)
                fails = 0
                threshold2 = min(n, limit2)
            next_p += step
            if threshold2 < next_p ** 2:
                return done(n, next_p)
    return done(n, next_p)


def _factor_limited(n, limit):
    """SymPy's result for ``factorint(n, limit)`` when the trial phase
    decides it; ``None`` when the full factorization is needed."""
    factors = {}
    small = min(2 ** 15, limit)
    n, next_p = _factorint_small(factors, n, small, 600)
    if next_p == 0:
        if n > 1:
            factors[int(n)] = 1
        return factors
    if next_p > limit:
        if n == 1:
            return factors
        if n < next_p ** 2 or isprime(n):
            factors[int(n)] = 1
            return factors
        pp = _core.nt_perfect_power(n, None, True)
        if pp is not None:
            base, exp = int(pp[0]), int(pp[1])
            if base < next_p ** 2 or isprime(base):
                factors[base] = exp
            else:
                for b, e in factorint(base, limit).items():
                    factors[b] = exp * e
            return factors
        factors[int(n)] = 1
        return factors
    return None


def _factor_dict(n):
    """Native complete factorization of an int ``n > 0`` as a dict."""
    return {int(p): int(e) for p, e in _core.nt_factor_list(n)}


def factorint(n, limit=None, use_trial=True, use_rho=True, use_pm1=True,
              use_ecm=True, verbose=False, visual=None, multiple=False):
    """Prime factorization of n as ``{p: e}`` (or a list / visual Mul)."""
    S, Mul, Pow, Rational, Integer = _sym()
    if multiple:
        fac = factorint(n, limit=limit, use_trial=use_trial, use_rho=use_rho,
                        use_pm1=use_pm1, verbose=verbose, visual=False, multiple=False)
        return sum(([p] * fac[p] if fac[p] > 0 else [S.One / p] * -fac[p]
                    for p in sorted(fac)), [])
    factordict = {}
    if visual and not (_is_mul(n) or isinstance(n, dict)):
        factordict = factorint(n, limit=limit, use_trial=use_trial, use_rho=use_rho,
                               use_pm1=use_pm1, verbose=verbose, visual=False)
    elif _is_mul(n):
        factordict = {int(k): int(v) for k, v in n.as_powers_dict().items()}
    elif isinstance(n, dict):
        factordict = n
    if factordict and (_is_mul(n) or isinstance(n, dict)):
        for key in list(factordict.keys()):
            if isprime(key):
                continue
            e = factordict.pop(key)
            d = factorint(key, limit=limit, use_trial=use_trial, use_rho=use_rho,
                          use_pm1=use_pm1, verbose=verbose, visual=False)
            for k, v in d.items():
                if k in factordict:
                    factordict[k] += v * e
                else:
                    factordict[k] = v * e
    if visual or (type(n) is dict and visual is not True and visual is not False):
        if factordict == {}:
            return S.One
        if -1 in factordict:
            factordict.pop(-1)
            args = [S.NegativeOne]
        else:
            args = []
        args.extend([Pow(*i, evaluate=False) for i in sorted(factordict.items())])
        return Mul(*args, evaluate=False)
    elif _is_mul(n) or isinstance(n, dict):
        return factordict
    n = as_int(n)
    if limit:
        limit = int(limit)
    if n < 0:
        factors = factorint(-n, limit=limit, use_trial=use_trial, use_rho=use_rho,
                            use_pm1=use_pm1, verbose=verbose, visual=False)
        factors[-1] = 1
        return factors
    if limit and limit < 2:
        if n == 1:
            return {}
        return {n: 1}
    elif n < 10:
        return [{0: 1}, {}, {2: 1}, {3: 1}, {2: 2}, {5: 1}, {2: 1, 3: 1},
                {7: 1}, {2: 3}, {3: 2}][n]
    if limit:
        r = _factor_limited(n, limit)
        if r is not None:
            return r
    return _factor_dict(n)


def factorrat(rat, limit=None, use_trial=True, use_rho=True, use_pm1=True,
              verbose=False, visual=None, multiple=False):
    """Factorization of a Rational as ``{p: e}`` with negative exponents."""
    S, Mul, Pow, Rational, Integer = _sym()
    if multiple:
        fac = factorrat(rat, limit=limit, use_trial=use_trial, use_rho=use_rho,
                        use_pm1=use_pm1, verbose=verbose, visual=False, multiple=False)
        return sum(([p] * fac[p] if fac[p] > 0 else [S.One / p] * -fac[p]
                    for p, _ in sorted(fac.items(),
                                       key=lambda elem: elem[0] if elem[1] > 0 else 1 / elem[0])), [])
    f = factorint(rat.p, limit=limit, use_trial=use_trial, use_rho=use_rho,
                  use_pm1=use_pm1, verbose=verbose).copy()
    f = defaultdict(int, f)
    for p, e in factorint(rat.q, limit=limit, use_trial=use_trial, use_rho=use_rho,
                          use_pm1=use_pm1, verbose=verbose).items():
        f[p] += -e
    if len(f) > 1 and 1 in f:
        del f[1]
    if not visual:
        return dict(f)
    if -1 in f:
        f.pop(-1)
        args = [S.NegativeOne]
    else:
        args = []
    args.extend([Pow(*i, evaluate=False) for i in sorted(f.items())])
    return Mul(*args, evaluate=False)


def primefactors(n, limit=None, verbose=False, **kwargs):
    """Sorted list of the prime factors of n (without -1, 0 and 1)."""
    n = int(n)
    kwargs.update({'visual': None, 'multiple': False, 'limit': limit, 'verbose': verbose})
    factors = sorted(factorint(n=n, **kwargs).keys())
    s = [f for f in factors[:-1] if f not in [-1, 0, 1]]
    if factors and isprime(factors[-1]):
        s += [factors[-1]]
    return s


def _divisors(n, proper=False):
    return iter([int(d) for d in _core.nt_divisors(n, bool(proper), False)])


def divisors(n, generator=False, proper=False):
    """Divisors of n, sorted (or a generator in SymPy's order)."""
    n = as_int(abs(n))
    if generator:
        return _divisors(n, proper)
    return [int(d) for d in _core.nt_divisors(n, bool(proper), True)]


def divisor_count(n, modulus=1, proper=False):
    """Number of divisors of n (divisible by ``modulus``)."""
    from .. import Mul
    if not modulus:
        return 0
    elif modulus != 1:
        n, r = divmod(n, modulus)
        if r:
            return 0
    if n == 0:
        return 0
    n = Mul(*[v + 1 for k, v in factorint(n).items() if k > 1])
    if n and proper:
        n -= 1
    return n


def proper_divisors(n, generator=False):
    """Divisors of n other than n itself."""
    return divisors(n, generator=generator, proper=True)


def proper_divisor_count(n, modulus=1):
    """Number of proper divisors of n."""
    return divisor_count(n, modulus=modulus, proper=True)


def _udivisors(n):
    if n <= 1:
        if n == 1:
            yield 1
        return
    factorpows = [p ** e for p, e in factorint(n).items()]
    for i in range(2 ** len(factorpows)):
        d = 1
        for k in range(i.bit_length()):
            if i & 1:
                d *= factorpows[k]
            i >>= 1
        yield d


def udivisors(n, generator=False):
    """Unitary divisors of n."""
    rv = _udivisors(as_int(abs(n)))
    return rv if generator else sorted(rv)


def udivisor_count(n):
    """Number of unitary divisors of n."""
    if n == 0:
        return 0
    return 2 ** len([p for p in factorint(n) if p > 1])


def _antidivisors(n):
    if n <= 2:
        return
    for d in _divisors(n):
        y = 2 * d
        if n > y and n % y:
            yield y
    for d in _divisors(2 * n - 1):
        if n > d >= 2 and n % d:
            yield d
    for d in _divisors(2 * n + 1):
        if n > d >= 2 and n % d:
            yield d


def antidivisors(n, generator=False):
    """Antidivisors of n."""
    rv = _antidivisors(as_int(abs(n)))
    return rv if generator else sorted(rv)


def antidivisor_count(n):
    """Number of antidivisors of n."""
    n = as_int(abs(n))
    if n <= 2:
        return 0
    return divisor_count(2 * n - 1) + divisor_count(2 * n + 1) + \
        divisor_count(n) - divisor_count(n, 2) - 5


def _divisor_sigma(n, k=1):
    return int(_core.nt_divisor_sigma(n, k))


def core(n, t=2):
    """t-th power free part of n."""
    n = as_int(n)
    t = as_int(t)
    if n <= 0:
        raise ValueError('n must be a positive integer')
    elif t <= 1:
        raise ValueError('t must be >= 2')
    y = 1
    for p, e in factorint(n).items():
        y *= p ** (e % t)
    return y


def mersenne_prime_exponent(nth):
    """The nth Mersenne prime exponent."""
    n = as_int(nth)
    if n < 1:
        raise ValueError('nth must be a positive integer; mersenne_prime_exponent(1) == 2')
    if n > 51:
        raise ValueError('There are only 51 perfect numbers; nth must be less than or equal to 51')
    return MERSENNE_PRIME_EXPONENTS[n - 1]


def is_perfect(n):
    """Whether n is a perfect number."""
    n = as_int(n)
    if n < 1:
        return False
    return _core.nt_is_perfect(n)


def abundance(n):
    """``sigma(n) - 2*n``."""
    return int(_core.nt_abundance(as_int(n)))


def is_abundant(n):
    """Whether n is an abundant number."""
    n = as_int(n)
    if is_perfect(n):
        return False
    return n % 6 == 0 or bool(abundance(n) > 0)


def is_deficient(n):
    """Whether n is a deficient number."""
    n = as_int(n)
    if is_perfect(n):
        return False
    return bool(abundance(n) < 0)


def is_amicable(m, n):
    """Whether m and n are an amicable pair."""
    return m != n and m + n == _divisor_sigma(as_int(m)) == _divisor_sigma(as_int(n))


def is_carmichael(n):
    """Whether n is a Carmichael number (``0`` for even ``n >= 561``)."""
    if n < 561:
        return False
    r = _core.nt_is_carmichael(as_int(n))
    if r is None:
        return n % 2 if isinstance(n, int) else _sym()[0].Zero
    return r


def find_carmichael_numbers_in_range(x, y):
    """Carmichael numbers in [x, y)."""
    if 0 <= x <= y:
        if x % 2 == 0:
            return [i for i in range(x + 1, y, 2) if is_carmichael(i)]
        return [i for i in range(x, y, 2) if is_carmichael(i)]
    raise ValueError('The provided range is not valid. x and y must be non-negative integers and x <= y')


def find_first_n_carmichaels(n):
    """The first n Carmichael numbers."""
    i = 561
    carmichaels = []
    while len(carmichaels) < n:
        if is_carmichael(i):
            carmichaels.append(i)
        i += 2
    return carmichaels


def dra(n, b):
    """Additive digital root of n in base b."""
    num = abs(as_int(n))
    b = as_int(b)
    if b <= 1:
        raise ValueError('Base should be an integer greater than 1')
    if num == 0:
        return 0
    return 1 + (num - 1) % (b - 1)


def drm(n, b):
    """Multiplicative digital root of n in base b."""
    n = abs(as_int(n))
    b = as_int(b)
    if b <= 1:
        raise ValueError('Base should be an integer greater than 1')
    while n > b:
        mul = 1
        while n > 1:
            n, r = divmod(n, b)
            if r == 0:
                return 0
            mul *= r
        n = mul
    return n


def _head(name):
    def f(*args):
        import symox
        return getattr(symox, name)(*args)
    f.__name__ = name
    f.__doc__ = "Deprecated in SymPy: the ``%s`` function head." % name
    return f


totient = _head('totient')
reduced_totient = _head('reduced_totient')
divisor_sigma = _head('divisor_sigma')
udivisor_sigma = _head('udivisor_sigma')
primenu = _head('primenu')
primeomega = _head('primeomega')
