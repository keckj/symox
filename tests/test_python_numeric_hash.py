"""Cross-type equality must obey Python's set and dictionary hash contract."""
from fractions import Fraction
import random
import sys

import pytest
import symox as sx


@pytest.mark.parametrize('value', [0, 1, -1, 2, -2, -(2**63), 2**63, 10**200, -10**201, sys.hash_info.modulus, sys.hash_info.modulus+1])
def test_integer_dictionary_keys(value):
    number = sx.Integer(value)
    assert number == value
    assert hash(number) == hash(value)
    assert len({number, value}) == 1
    assert {number: 'symbolic'}[value] == 'symbolic'
    assert {value: 'python'}[number] == 'python'


def test_rational_hashes_match_python_fractions():
    rng = random.Random(42017)
    p = sys.hash_info.modulus
    values = [Fraction(1, 2), Fraction(-1, 2), Fraction(1, p), Fraction(-1, p), Fraction(p, p*p), Fraction(0)]
    values += [Fraction(rng.randrange(-10**100, 10**100), rng.randrange(1, 10**100)) for _ in range(100)]
    for value in values:
        number = sx.Rational(value.numerator, value.denominator)
        assert number == value
        assert hash(number) == hash(value)
        assert len({number, value}) == 1


@pytest.mark.parametrize('value', [0.0, -0.0, 0.5, -0.5, 1.25, 1e300, 5e-324, float('inf'), float('-inf'), True, False])
def test_float_boolean_and_infinity_keys(value):
    number = sx.sympify(value)
    assert number == value
    assert hash(number) == hash(value)
    assert len({number, value}) == 1


def test_multiprecision_hash_uses_the_exact_binary_value():
    value = Fraction(1, 2**400)
    number = sx.Float('0.5', 150)**400
    assert hash(number) == hash(value)
    assert hash(sx.Float('-0.5', 150)) == hash(Fraction(-1, 2))
