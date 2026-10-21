"""factor_list and sqf_list of expressions with symbolic or fractional
multiplicities, and the relational ordering their sort relies on.

The expected values were produced by the reference SymPy 1.14. The
multiplicities keep their expressions, the multiplicities of a component of
exponent one are Python ints, and the factor sort raises the TypeError of
an undecided comparison where the reference raises it. The square-free
grouping order of two symbolic multiplicities depends on the hash seed in
the reference and is not tested.
"""
import pytest
from symox import *

# the star import brings the library test runner, which is not a test here
del test

a, z, p = symbols('a z p')
s = Symbol('s')
n = Symbol('n', positive=True)
x = Symbol('x', real=True)
w = Symbol('w', finite=True)

VALUES = [
    ('factor_list(sqrt(z + 1)*(z - 1)**2)', '(1, [(z + 1, 1/2), (z - 1, 2)])', ['Half', 'Integer']),
    ('sqf_list(sqrt(z + 1)*(z - 1)**2)', '(1, [(z + 1, 1/2), (z - 1, 2)])', ['Half', 'Integer']),
    ('factor_list(z*(z**2 - 1)**s)', '(1, [(z, 1), (z**2 - 1, s)])', ['int', 'Symbol']),
    ('sqf_list(2**s*(z + 1)**s)', '(1, [(2*z + 2, s)])', ['Symbol']),
    ('factor_list(sqrt(z + 1)*z)', '(1, [(z + 1, 1/2), (z, 1)])', ['Half', 'int']),
    ('sqf_list(sqrt(z + 1)*z)', '(1, [(z + 1, 1/2), (z, 1)])', ['Half', 'int']),
    ('factor_list((z + 1)**(S(3)/2)*z**2)', '(1, [(z + 1, 3/2), (z, 2)])', ['Rational', 'Integer']),
    ('sqf_list((z + 1)**(S(3)/2)*z**2)', '(1, [(z + 1, 3/2), (z, 2)])', ['Rational', 'Integer']),
    ('factor_list((z + 1)**s*(z - 1)**s)', '(1, [(z - 1, s), (z + 1, s)])', ['Symbol', 'Symbol']),
    ('sqf_list((z + 1)**s*(z - 1)**s)', '(1, [(z**2 - 1, s)])', ['Symbol']),
    ('Lt(s, s + 1)', 's < s + 1', None),
    ('Lt(x, x + 1)', 'True', None),
    ('Lt(w, w + 1)', 'w < w + 1', None),
    ('Gt(s + 1, s)', 's + 1 > s', None),
    ('Le(n, n + 2)', 'True', None),
    ('Lt(1, s)', '1 < s', None),
    ('Eq(s, s + 1)', 'False', None),
]

ERRORS = [
    ('factor_list(a**s*(z + 1))', 'TypeError', 'cannot determine truth value of Relational: s > 1'),
    ('sqf_list(a**s*(z + 1))', 'TypeError', 'cannot determine truth value of Relational: s < 1'),
    ('factor_list(a**s*(z + 1)**2)', 'TypeError', 'cannot determine truth value of Relational: 2 < s'),
    ('sqf_list(a**s*(z + 1)**2)', 'TypeError', 'cannot determine truth value of Relational: s < 2'),
    ('factor_list((z + 1)**s*(z - 1))', 'TypeError', 'cannot determine truth value of Relational: s > 1'),
    ('sqf_list((z + 1)**s*(z - 1))', 'TypeError', 'cannot determine truth value of Relational: s < 1'),
    ('sqf_list(z*(z**2 - 1)**s)', 'TypeError', 'cannot determine truth value of Relational: s < 1'),
    ('factor_list(p**s*(z + 1)**3)', 'TypeError', 'cannot determine truth value of Relational: 3 < s'),
    ('sqf_list(p**s*(z + 1)**3)', 'TypeError', 'cannot determine truth value of Relational: s < 3'),
    ('factor_list(2**s*(z + 1)**s)', 'PolificationFailed', 'Cannot construct a polynomial from 2'),
    ('factor_list(a**(s + 1)*(z + 1)**s)', 'TypeError', 'cannot determine truth value of Relational: s < s + 1'),
    ('factor_list(sqrt(z)*sqrt(z + 1)/(z - 1))', 'PolynomialError', 'a polynomial expected, got sqrt(z)*sqrt(z + 1)/(z - 1)'),
    ('sqf_list(sqrt(z)*sqrt(z + 1)/(z - 1))', 'PolynomialError', 'a polynomial expected, got sqrt(z)*sqrt(z + 1)/(z - 1)'),
    ('factor_list(2**s)', 'PolificationFailed', 'Cannot construct a polynomial from 2'),
    ('sqf_list(2**s)', 'PolificationFailed', 'Cannot construct a polynomial from 2'),
    ('factor_list((z + 1)**n*(z - 1))', 'TypeError', 'cannot determine truth value of Relational: n > 1'),
    ('sqf_list((z + 1)**n*(z - 1))', 'TypeError', 'cannot determine truth value of Relational: n < 1'),
    ('factor_list(a**n*(z + 1))', 'TypeError', 'cannot determine truth value of Relational: n > 1'),
    ('sqf_list(a**n*(z + 1))', 'TypeError', 'cannot determine truth value of Relational: n < 1'),
    ('factor_list((-z - 1)**s*(z - 1))', 'TypeError', 'cannot determine truth value of Relational: s > 1'),
    ('sqf_list((-z - 1)**s*(z - 1))', 'TypeError', 'cannot determine truth value of Relational: s < 1'),
]


@pytest.mark.parametrize('call, expected, kinds', VALUES)
def test_symbolic_multiplicities_match_the_reference(call, expected, kinds):
    result = eval(call)
    assert str(result) == expected
    if kinds is not None:
        assert [type(k) is int for _, k in result[1]] == [k == 'int' for k in kinds]


@pytest.mark.parametrize('call, error, message', ERRORS)
def test_symbolic_multiplicity_errors_match_the_reference(call, error, message):
    with pytest.raises(Exception) as info:
        eval(call)
    assert type(info.value).__name__ == error
    assert str(info.value) == message
