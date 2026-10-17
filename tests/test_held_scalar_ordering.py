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


@pytest.mark.parametrize('seed,extra', [
    (756913, ()), (756912, ('oo', '-oo', 'zoo', 'nan')),
])
def test_seeded_held_signs_and_keys_match_oracle_without_canonicalizing_input(seed, extra):
    randomizer = random.Random(seed)
    def generate(depth):
        if not depth or randomizer.random() < .45:
            return randomizer.choice(['x', 'y', 'z', '1', '-1', '2', '-2', 'Rational(1,2)', 'I', 'sin(x)'] + list(extra))
        children = [generate(depth-1) for _ in range(randomizer.choice([2, 3]))]
        head = 'Add' if randomizer.random() < .65 else 'Mul'
        return f'{head}({",".join(children)},evaluate=False)'
    for _ in range(1000):
        source = generate(2)
        actual, expected = expression(sx, source), expression(sp, source)
        actual, expected = sx.sympify(actual), sp.sympify(expected)
        try:
            expected_sign = expected.could_extract_minus_sign()
        except (TypeError, ValueError):
            # SymPy can reject its own NaN comparisons while choosing a sign.
            continue
        assert actual.could_extract_minus_sign() == expected_sign, source
        try:
            expected_key = sp.default_sort_key(expected)
        except (TypeError, ValueError):
            # A rejected oracle key is not a reference ordering.
            continue
        assert normalized(sx.default_sort_key(actual)) == expected_key, source


@pytest.mark.parametrize('source', [
    'Add(-2,Mul(y,-2,I,evaluate=False),evaluate=False)',
    'Add(y,Mul(-1,x,Rational(1,2),evaluate=False),evaluate=False)',
    'Add(Mul(-1,x,I,evaluate=False),Mul(x,y,-2,evaluate=False),evaluate=False)',
])
def test_held_sum_reconstruction_preserves_product_signs(source):
    assert expression(sx, source).could_extract_minus_sign() == expression(sp, source).could_extract_minus_sign()


@pytest.mark.parametrize('source', [
    'Add(Rational(1,2),Mul(2,1,evaluate=False),evaluate=False)',
    'Add(x,Mul(2,2,Rational(1,2),evaluate=False),Mul(x,x,-1,evaluate=False),evaluate=False)',
])
def test_held_sum_reconstruction_preserves_product_presentation(source):
    assert normalized(sx.default_sort_key(expression(sx, source))) == sp.default_sort_key(expression(sp, source))


@pytest.mark.parametrize('source', SOURCES + [
    'Mul(x,oo,evaluate=False)', 'Mul(x,-oo,evaluate=False)',
    'Mul(x,zoo,evaluate=False)', 'Mul(x,nan,evaluate=False)',
    'Add(Add(1,1,evaluate=False),2,evaluate=False)',
    'Add(x,Add(y,2,evaluate=False),evaluate=False)',
])
def test_negation_of_held_arithmetic_matches_oracle_without_changing_input(source):
    actual, expected = expression(sx, source), expression(sp, source)
    before = actual.args
    assert normalized(sx.default_sort_key(-actual)) == sp.default_sort_key(-expected)
    assert actual.args == before


@pytest.mark.parametrize('source', SOURCES + [
    'Mul(x,oo,evaluate=False)', 'Mul(x,-oo,evaluate=False)',
    'Add(Mul(x,oo,evaluate=False),y,evaluate=False)',
])
def test_public_term_order_uses_the_same_keys_as_default_sort_key(source):
    actual, expected = expression(sx, source), expression(sp, source)
    assert [normalized(sx.default_sort_key(term)) for term in actual.as_ordered_terms()] == [sp.default_sort_key(term) for term in expected.as_ordered_terms()]


def test_public_arithmetic_preserves_dag_sharing_for_held_sum_normalization():
    one = sx.Integer(1)
    x = sx.Symbol('x')
    value = x
    for _ in range(80):
        value = sx._core._restore_arithmetic('Add', sx.Integer(0), [(value, one), (value, one)], False)
    assert value + 1 == 2**80*x+1
    assert -value == -(2**80*x)
    assert 2*value == 2**81*x
    assert sx.Add(sx.Add(1,1,evaluate=False),2) == 4


@pytest.mark.parametrize('source', [
    'Add(-2,Add(I,2,oo,evaluate=False),evaluate=False)',
    'Add(-1,Add(z,-oo,evaluate=False),evaluate=False)',
    'Add(Add(I,-oo,I,evaluate=False),Mul(-1,1,evaluate=False),evaluate=False)',
    'Add(Add(Rational(1,2),zoo,y,evaluate=False),-oo,evaluate=False)',
])
def test_nested_nonfinite_sum_signs_reclassify_flattened_infinities(source):
    assert expression(sx, source).could_extract_minus_sign() == expression(sp, source).could_extract_minus_sign()


@pytest.mark.parametrize('source', [
    '(oo+x)-oo', '2*(oo+x)', '-(oo+x)', '(oo+I)+2',
    'oo+Symbol("r",real=True)', 'zoo+I',
])
def test_nonfinite_arithmetic_uses_exact_finiteness_facts(source):
    assert normalized(sx.default_sort_key(expression(sx, source))) == sp.default_sort_key(expression(sp, source))


@pytest.mark.parametrize('source', [
    'Add(Add(Rational(1,2),oo,evaluate=False),Mul(oo,-2,evaluate=False),evaluate=False)',
    'Add(nan,I,Add(-1,Rational(1,2),-oo,evaluate=False),evaluate=False)',
    'Add(Add(x,I,2,evaluate=False),Mul(-2,I,evaluate=False),Add(-1,-oo,1,evaluate=False),evaluate=False)',
])
def test_nonfinite_numeric_coefficients_stop_comparison_at_unequal_nan(source):
    actual, expected = expression(sx, source), expression(sp, source)
    assert normalized(sx.default_sort_key(actual)) == sp.default_sort_key(expected)



def test_held_numeric_power_constructor_does_not_evaluate_large_integer():
    import subprocess
    import sys
    subprocess.run([sys.executable, "-c", """
import symox as s
power = s.Pow(2, 10**9, evaluate=False)
assert power.args == (2, 10**9)
assert 0*power == 0
assert s.Pow(0, 1, evaluate=False)*power == 0
"""], check=True, timeout=10)
