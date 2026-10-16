"""Held scalar views keep raw tails while canonical comparison stays deterministic."""
import random

import pytest
import sympy as sp
import symox as sx

from test_factor_inspection import expression
from test_sortkey_compat import normalized


SOURCES = [
    'Mul(-1,x,Rational(1,2),evaluate=False)',
    'Mul(2,x,x,evaluate=False)', 'Mul(-2,x,y,evaluate=False)',
    'Mul(Float(-2.5),x,x,evaluate=False)',
    'Mul(x,2,y,evaluate=False)', 'Mul(-oo,x,evaluate=False)',
    'Mul(x,UnevaluatedExpr(y),z,evaluate=False)',
    'Add(y,Add(y,-2,-1,evaluate=False),evaluate=False)',
    'Add(Mul(-1,2,zoo,evaluate=False),Add(-oo,oo,x,evaluate=False),evaluate=False)',
    'Add(oo,Mul(-1,sin(x),-1,evaluate=False),evaluate=False)',
    'Pow(Rational(1,2),x,evaluate=False)', 'Pow(E,x,evaluate=False)',
    'oo', '-oo', 'zoo', 'nan', 'Float(-2.5)', 'I',
]


@pytest.mark.parametrize('source', SOURCES)
def test_scalar_key_preserves_held_structure_and_exact_nonfinite_atoms(source):
    actual, expected = expression(sx, source), expression(sp, source)
    assert normalized(sx.default_sort_key(actual)) == sp.default_sort_key(expected)


@pytest.mark.parametrize('source', SOURCES)
@pytest.mark.parametrize('rational', [False, True])
def test_leading_coefficient_keeps_raw_product_tail(source, rational):
    actual, expected = expression(sx, source), expression(sp, source)
    before = actual.args
    left = actual.as_coeff_Mul(rational=rational)
    right = expected.as_coeff_Mul(rational=rational)
    assert normalized(sx.default_sort_key(left[0])) == sp.default_sort_key(right[0])
    assert normalized(sx.default_sort_key(left[1])) == sp.default_sort_key(right[1])
    assert actual.args == before


@pytest.mark.parametrize('source', SOURCES[:7])
def test_ordered_factors_split_negative_coefficient_and_keep_noncommutative_tail(source):
    actual, expected = expression(sx, source), expression(sp, source)
    assert [normalized(sx.default_sort_key(a)) for a in actual.as_ordered_factors()] == [sp.default_sort_key(a) for a in expected.as_ordered_factors()]


def test_scalar_key_keeps_shared_python_tuples_without_expanding_dag():
    node = sx.Symbol('x')
    pair = sx.Function('pair')
    for _ in range(80):
        node = pair(node, node)
    key = sx.default_sort_key(node)
    seen = set()
    pending = [key]
    while pending:
        item = pending.pop()
        if not isinstance(item, tuple) or id(item) in seen:
            continue
        seen.add(id(item))
        pending.extend(item)
    assert len(seen) < 650
    # Repeated function arguments must refer to the very same key tuple.
    while key[0][2] == 'pair':
        children = key[1][1]
        assert children[0] is children[1]
        key = children[0]


def test_seeded_held_signs_match_oracle_without_canonicalizing_input():
    randomizer = random.Random(756913)
    def generate(depth):
        if not depth or randomizer.random() < .45:
            return randomizer.choice(['x', 'y', 'z', '1', '-1', '2', '-2', 'Rational(1,2)', 'I', 'sin(x)'])
        children = [generate(depth-1) for _ in range(randomizer.choice([2, 3]))]
        head = 'Add' if randomizer.random() < .65 else 'Mul'
        return f'{head}({",".join(children)},evaluate=False)'
    for _ in range(400):
        source = generate(2)
        actual, expected = expression(sx, source), expression(sp, source)
        actual, expected = sx.sympify(actual), sp.sympify(expected)
        assert actual.could_extract_minus_sign() == expected.could_extract_minus_sign(), source


@pytest.mark.xfail(strict=True, reason='Canonical Add reconstruction still changes the raw tail of a held product')
@pytest.mark.parametrize('source', [
    'Add(-2,Mul(y,-2,I,evaluate=False),evaluate=False)',
    'Add(y,Mul(-1,x,Rational(1,2),evaluate=False),evaluate=False)',
    'Add(Mul(-1,x,I,evaluate=False),Mul(x,y,-2,evaluate=False),evaluate=False)',
])
def test_remaining_held_sum_reconstruction_sign_cases(source):
    assert expression(sx, source).could_extract_minus_sign() == expression(sp, source).could_extract_minus_sign()


@pytest.mark.xfail(strict=True, reason='Held Add may reuse a canonical parent whose product children have different public argument order')
@pytest.mark.parametrize('source', [
    'Add(Rational(1,2),Mul(2,1,evaluate=False),evaluate=False)',
    'Add(x,Mul(2,2,Rational(1,2),evaluate=False),Mul(x,x,-1,evaluate=False),evaluate=False)',
])
def test_remaining_held_sum_reconstruction_key_cases(source):
    assert normalized(sx.default_sort_key(expression(sx, source))) == sp.default_sort_key(expression(sp, source))
