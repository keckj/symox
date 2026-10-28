"""Ordinary nonfinite arithmetic: reference agreement and explicit undefined sums.

A zero times an Add with a literal infinite public term/factor is NaN.
This rule stops at functions, powers, and nested sums. SymPy1.14 has an
arity-dependent exception for some n-ary noncommutative constructors; those
cases deliberately follow the same undefined-product rule in every order.
"""
import itertools
import pytest
import sympy as sp
import symox as sx


def _infinite_sum(module, sign, commutative):
    b = module.Symbol('b', commutative=commutative)
    infinity = {'positive': module.oo, 'negative': -module.oo, 'complex': module.zoo}[sign]
    return b, module.Add(infinity, b, evaluate=False)


@pytest.mark.parametrize('commutative', [True, False])
@pytest.mark.parametrize('sign', ['positive', 'negative', 'complex'])
@pytest.mark.parametrize('reverse', [False, True])
def test_zero_times_literal_infinite_sum_matches_reference(commutative, sign, reverse):
    for module in (sp, sx):
        _, term = _infinite_sum(module, sign, commutative)
        args = (term, module.S.Zero) if reverse else (module.S.Zero, term)
        assert module.Mul(*args) == module.nan


@pytest.mark.parametrize('sign', ['positive', 'negative', 'complex'])
def test_zero_product_is_undefined_in_every_factor_order_and_grouping(sign):
    _, term = _infinite_sum(sx, sign, False)
    # A literal infinity cannot produce a finite sum: cancellation with an
    # opposing infinity is itself undefined. Thus none of these products is0.
    for factors in itertools.permutations([sx.S.Zero, sx.Integer(2), term]):
        assert sx.Mul(*factors) == sx.nan
    assert sx.Mul(0, term, term) == sx.nan
    b, c = sx.symbols('b c', commutative=False)
    held = sx.Mul(c, term, evaluate=False)
    assert sx.Mul(0, held) == sx.nan
    assert sx.Mul(held, 0) == sx.nan


@pytest.mark.parametrize('sign', ['positive', 'negative', 'complex'])
def test_infinite_factor_in_a_direct_sum_term_is_also_undefined(sign):
    for module in (sp, sx):
        b, c = module.symbols('b c', commutative=False)
        infinity = {'positive': module.oo, 'negative': -module.oo, 'complex': module.zoo}[sign]
        inner = module.Mul(infinity, b, evaluate=False)
        held = module.Add(inner, c, evaluate=False)
        assert module.Mul(0, held) == module.nan


@pytest.mark.parametrize('boundary', ['function', 'power', 'nested_sum', 'ordinary_symbol'])
def test_literal_scan_preserves_reference_boundaries(boundary):
    for module in (sp, sx):
        b, c = module.symbols('b c', commutative=False)
        held = module.Add(module.oo, b, evaluate=False)
        term = {'function': lambda: module.exp(held, evaluate=False),
                'power': lambda: module.Pow(held, 2, evaluate=False),
                'nested_sum': lambda: module.Add(held, c, evaluate=False),
                'ordinary_symbol': lambda: b}[boundary]()
        assert module.Mul(0, term) == 0


@pytest.mark.parametrize('value', [sx.oo, -sx.oo, sx.zoo, sx.nan, 0, 1, -1, sx.I])
def test_partial_simultaneous_and_reversed_substitution_cannot_leave_a_false_zero(value):
    a = sx.Symbol('a')
    b, term = _infinite_sum(sx, 'positive', False)
    product = a*term
    first = product.subs(a, 0).subs(b, value)
    second = product.subs(b, value).subs(a, 0)
    simultaneous = product.subs({a: 0, b: value}, simultaneous=True)
    assert first == second == simultaneous == sx.nan


@pytest.mark.parametrize('assumptions', [{}, {'real': True}, {'finite': True}, {'complex': True},
                                         {'extended_real': True}, {'finite': False}])
@pytest.mark.parametrize('sign', [1, -1])
def test_addition_does_not_assume_an_unknown_scalar_finite(assumptions, sign):
    actual = sx.Symbol('b', **assumptions)+sign*sx.oo
    reference = sp.Symbol('b', **assumptions)+sign*sp.oo
    assert sx.to_sympy(actual) == reference
    for native, expected in [(sx.oo, sp.oo), (-sx.oo, -sp.oo), (sx.zoo, sp.zoo),
                             (sx.nan, sp.nan), (sx.I, sp.I), (sx.S.One, sp.S.One)]:
        assert sx.to_sympy(actual.subs(sx.Symbol('b', **assumptions), native)) == reference.subs(sp.Symbol('b', **assumptions), expected)


def test_nonzero_product_retains_high_precision_and_ordered_factors():
    b, c = sx.symbols('b c', commutative=False)
    text = '1.12345678901234567890123456789012345678901234567890123456789012345678901234567890'
    coefficient = sx.Float(text, 90)
    result = sx.Mul(coefficient, b, c)
    assert result.args == (coefficient, b, c)
    assert sx.to_sympy(result.args[0])._prec >= 280


def test_an_extended_real_exponential_can_be_zero_at_a_valid_infinite_argument():
    b = sx.Symbol('b', extended_real=True)
    product = sx.oo*sx.exp(b)
    # -oo satisfies extended_real; exp(-oo)=0 makes the product undefined.
    # SymPy1.14 absorbs exp(b) too early here, so this is an exact domain
    # regression rather than an assertion matching that reference defect.
    assert product.subs(b, -sx.oo) == sx.oo*sx.exp(-sx.oo) == sx.nan
