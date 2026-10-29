"""Big-integer gcd and Rational normalisation (native Lehmer gcd).

The digests are sha256(str(value))[:16] and lengths of the SymPy 1.14
results (pure Python ground types) on smaller instances of the audit
witnesses (audits/performance/2026-10-29-biggcd). The randomized checks
are certified by Python's exact ``fractions.Fraction`` and ``math.gcd``.
"""
import functools
import hashlib
import math
import random
import sys
from fractions import Fraction

import pytest

from symox import (Add, Integer, Poly, Rational, Symbol, gamma, gcd, harmonic,
                   lcm)

x = Symbol('x')
G = 13**3000
A = G * 17**4000 + G
B = G * 19**3800 * 3 + 5 * G


def digest(value):
    old = sys.get_int_max_str_digits()
    sys.set_int_max_str_digits(0)
    try:
        text = str(value)
    finally:
        sys.set_int_max_str_digits(old)
    return len(text), hashlib.sha256(text.encode()).hexdigest()[:16]


SYMPY_DIGESTS = [
    ('rat_pow2', lambda: Rational(3**20001, 2**10000), (12555, '1408acf78eeb2008')),
    ('rat_pow2_shift', lambda: Rational(3**20001 * 2**7, -2**10000), (12554, '92aafd888ba69e47')),
    ('gamma_half', lambda: gamma(Rational(20001, 2)), (41689, 'ea6dea97eaf43f2c')),
    ('gcd_int', lambda: gcd(Integer(A), Integer(B)), (3343, '40a108e4250ea8b3')),
    ('lcm_int', lambda: lcm(Integer(A), Integer(B)), (13124, '52fc42e59337c5f0')),
    ('rat_AB', lambda: Rational(A, B), (9783, '22295814df909533')),
    ('div_AB', lambda: Integer(A) / Integer(-B), (9784, '1aad255e746c01bc')),
    ('primitive', lambda: Poly(A*x**2 + B*x + G, x).primitive(), (13164, '7805eb267517d1a2')),
    ('add_200', lambda: Add(*[Rational(1, 10**40 + k) for k in range(200)]), (15325, 'e83eccce05b2997a')),
    ('reduce_200', lambda: functools.reduce(lambda a, b: a + b, [Rational(k, 10**40 + k) for k in range(1, 201)]),
     (15313, '1247a25e759b659c')),
    ('harmonic_2000', lambda: harmonic(2000), (1734, '0d32fa67c85c913d')),
    ('harmonic_300_3', lambda: harmonic(300, 3), (779, 'fa21283480ccf71b')),
    ('mul_chain', lambda: functools.reduce(lambda a, b: a * b, [Rational(10**30 + k, 10**30 - k) for k in range(1, 101)]),
     (5623, 'd7d9fbb1e3cb1681')),
]


@pytest.mark.parametrize('name,build,expected', SYMPY_DIGESTS, ids=[c[0] for c in SYMPY_DIGESTS])
def test_matches_sympy_digest(name, build, expected):
    assert digest(build()) == expected


def test_power_of_two_denominator_is_exact():
    a = 3**20001
    r = Rational(a * 2**37, 2**10000)
    assert (r.p, r.q) == (a, 2**(10000 - 37))
    r = Rational(-(2**10000), 2**9000 * 5)
    assert (r.p, r.q) == (-(2**1000), 5)
    assert Rational(2**5000, 2**5000) == 1 and Rational(2**5000, 2**5000).is_Integer


def test_signs_zero_and_integer_results():
    big = 10**60 + 7
    assert Rational(0, -big) == 0 and Rational(0, -big).is_zero
    r = Rational(big * 6, -big * 4)
    assert (r.p, r.q) == (-3, 2)
    assert Rational(1, 3 * big) - Rational(1, 3 * big) == 0
    p = Rational(big, 3) * Rational(3, big)
    assert p == 1 and p.is_Integer
    q = Rational(big**2, 7) / Rational(big, 7)
    assert (q.p, q.q) == (big, 1)
    s = Rational(1, 2 * big) + Rational(1, 2 * big)
    assert (s.p, s.q) == (1, big)


def _rand_big(rng, digits):
    return rng.randrange(10**(digits - 1), 10**digits)


def test_random_rational_arithmetic_matches_fraction():
    rng = random.Random(20261029)
    common = [1, 2**64 + 13, 3**90, _rand_big(rng, 60)]
    for case in range(120):
        c = common[case % len(common)]
        n1 = (-1)**case * _rand_big(rng, 1 + case % 70) * c
        d1 = _rand_big(rng, 1 + (case * 7) % 50) * (c if case % 3 else 1) << (case % 90)
        n2 = _rand_big(rng, 1 + (case * 3) % 80) * (c if case % 2 else 1)
        d2 = -_rand_big(rng, 1 + (case * 5) % 60) * c << (case % 40)
        for value, ref in [
            (Rational(n1, d1), Fraction(n1, d1)),
            (Rational(n1, d1) + Rational(n2, d2), Fraction(n1, d1) + Fraction(n2, d2)),
            (Rational(n1, d1) - Rational(n2, d2), Fraction(n1, d1) - Fraction(n2, d2)),
            (Rational(n1, d1) * Rational(n2, d2), Fraction(n1, d1) * Fraction(n2, d2)),
            (Rational(n1, d1) / Rational(n2, d2), Fraction(n1, d1) / Fraction(n2, d2)),
            (Integer(n1) / Integer(d2), Fraction(n1, d2)),
        ]:
            assert (int(value.p), int(value.q)) == (ref.numerator, ref.denominator)


def test_random_integer_gcd_lcm_match_math():
    rng = random.Random(7)
    for case in range(80):
        g = _rand_big(rng, 1 + case % 300) << (case % 70)
        a = _rand_big(rng, 1 + (case * 11) % 400) * g * (-1)**case
        b = _rand_big(rng, 1 + (case * 13) % 350) * g
        if case % 9 == 0:
            b = 2**(case * 17)
        assert int(gcd(Integer(a), Integer(b))) == math.gcd(a, b)
        assert int(lcm(Integer(a), Integer(b))) == math.lcm(a, b)
        poly = Poly(a * x**2 + b * x + g, x)
        content, prim = poly.primitive()
        assert int(content) == math.gcd(a, b, g)
