"""Helpers the ported dense routines import from the reference package."""
import random as _random

from ...polys.polyconfig import query

GROUND_TYPES = 'python'
SYMPY_INTS = (int,)
MPZ = int

_rng = _random.Random()
uniform = _rng.uniform


def invert(x, m):
    try:
        return pow(x, -1, m)
    except ValueError:
        raise ZeroDivisionError("invert() no inverse exists")


def _randint(seed=None):
    if seed is None:
        return _rng.randint
    if isinstance(seed, int):
        _rng.seed(seed)
        return _rng.randint
    seed = list(seed)
    seed.reverse()

    def give(a, b, seq=seed):
        x = seq.pop()
        if not a <= x <= b:
            raise ValueError('given sequence value %s is not in range [%s, %s]' % (x, a, b))
        return x
    return give


def _order_key(factor):
    if isinstance(factor, list):
        return [_order_key(f) for f in factor]
    try:
        return int(factor)
    except TypeError:
        return factor


def _sort_factors(factors, **args):
    if args.get('multiple', True):
        return sorted(factors, key=lambda fn: (len(fn[0]), fn[1], _order_key(fn[0])))
    return sorted(factors, key=lambda f: (len(f), _order_key(f)))


class _IntegerRing:
    """The reference ``ZZ`` over Python ints (the ground of ``GF(n)``)."""
    dtype = tp = int
    zero, one = 0, 1
    is_ZZ = is_IntegerRing = is_Ring = is_PID = is_Exact = is_Numerical = True
    is_Field = is_QQ = is_FiniteField = is_FF = False
    has_assoc_Field = True

    def __call__(self, a):
        return int(a)

    def __str__(self):
        return 'ZZ'

    __repr__ = __str__

    def __eq__(self, other):
        return isinstance(other, _IntegerRing) or (getattr(other, 'is_ZZ', False) is True)

    def __hash__(self):
        return hash('ZZ')

    def convert(self, a, base=None):
        return int(a)

    def from_ZZ(self, a, base=None):
        return int(a)

    def of_type(self, a):
        return type(a) is int

    def is_one(self, a):
        return a == 1

    def is_zero(self, a):
        return not a

    def is_negative(self, a):
        return a < 0

    def is_nonnegative(self, a):
        return a >= 0

    def canonical_unit(self, a):
        return -1 if a < 0 else 1

    def gcd(self, a, b):
        from math import gcd
        return gcd(a, b)

    def lcm(self, a, b):
        from math import lcm
        return lcm(a, b)

    def gcdex(self, a, b):
        from ..._modular_ring import _gcdex
        return _gcdex(a, b)

    def invert(self, a, b):
        s, t, h = self.gcdex(a, b)
        if h == 1:
            return s % b
        from ..polyerrors import NotInvertible
        raise NotInvertible("zero divisor")

    def quo(self, a, b):
        return a // b

    def exquo(self, a, b):
        if a % b:
            from ..polyerrors import ExactQuotientFailed
            raise ExactQuotientFailed(a, b, self)
        return a // b

    def rem(self, a, b):
        return a % b

    def div(self, a, b):
        return divmod(a, b)

    def to_sympy(self, a):
        from ... import Integer
        return Integer(a)

    def sqrt(self, a):
        from math import isqrt
        return isqrt(a)


ZZ = _IntegerRing()


def __getattr__(name):
    if name == 'QQ':
        from ..._core import QQ
        return QQ
    raise AttributeError(name)
