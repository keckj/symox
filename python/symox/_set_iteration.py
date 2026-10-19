"""Exact enumerations of the countable numeric sets."""
from itertools import count
from math import gcd


def _integers(self):
    from . import Integer
    yield Integer(0)
    for value in count(1):
        yield Integer(value)
        yield Integer(-value)


def _naturals(self):
    from . import Integer
    for value in count(0 if self.class_name == 'Naturals0' else 1):
        yield Integer(value)


def _rationals(self):
    from . import Integer, Rational
    yield Integer(0)
    yield Integer(1)
    yield Integer(-1)
    # Coprime pairs enumerate each positive rational below/above one once.
    for denominator in count(2):
        for numerator in range(1, denominator):
            if gcd(numerator, denominator) == 1:
                yield Rational(numerator, denominator)
                yield Rational(denominator, numerator)
                yield Rational(-numerator, denominator)
                yield Rational(-denominator, numerator)


def install(namespace):
    for name, iterator in [('Integers', _integers), ('Naturals', _naturals),
                           ('Naturals0', _naturals), ('Rationals', _rationals)]:
        type(namespace[name]).__iter__ = iterator
