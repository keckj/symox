"""Certified algebraic zero times zoo after radical rationalization."""
import pytest
import symox as S
import sympy as R


def original_example(library):
    sqrt = library.sqrt
    return (-240*sqrt(2)*sqrt(sqrt(5) + 5)*sqrt(8*sqrt(5) + 40) -
        360*sqrt(2)*sqrt(-8*sqrt(5) + 40)*sqrt(-sqrt(5) + 5) -
        120*sqrt(10)*sqrt(-8*sqrt(5) + 40)*sqrt(-sqrt(5) + 5) +
        120*sqrt(2)*sqrt(-sqrt(5) + 5)*sqrt(8*sqrt(5) + 40) +
        120*sqrt(2)*sqrt(-8*sqrt(5) + 40)*sqrt(sqrt(5) + 5) +
        120*sqrt(10)*sqrt(-sqrt(5) + 5)*sqrt(8*sqrt(5) + 40) +
        120*sqrt(10)*sqrt(-8*sqrt(5) + 40)*sqrt(sqrt(5) + 5))/(-36000 -
        7200*sqrt(5) + (12*sqrt(10)*sqrt(sqrt(5) + 5) +
        24*sqrt(10)*sqrt(-sqrt(5) + 5))**2)


@pytest.mark.parametrize('symbolic', [False, True])
@pytest.mark.parametrize('max_terms', [0, 1, 4])
def test_original_radical_zero_over_zero(symbolic, max_terms):
    value = original_example(S)
    expected = original_example(R)
    n, d = R.fraction(expected)
    assert R.simplify(n) == 0 and R.simplify(d) == 0
    assert R.radsimp(expected, symbolic=symbolic, max_terms=max_terms) is R.nan
    assert S.radsimp(value, symbolic=symbolic, max_terms=max_terms) is S.nan


@pytest.mark.parametrize('wrapper', [lambda e: e, lambda e: 1+e, lambda e: S.sin(e)])
def test_nested_undefined_product_is_preserved(wrapper):
    assert S.radsimp(wrapper(original_example(S))) is S.nan


def test_unknown_symbolic_zero_is_not_inferred():
    x = S.Symbol('x')
    value = S.radsimp(S.zoo*x)
    assert value is not S.nan
    assert value.free_symbols == {x}
    assert value.subs(x, 1) is S.zoo
    assert value.subs(x, 0) is S.nan


@pytest.mark.parametrize('power', [20, 80, 200])
def test_near_zero_algebraic_numerator_is_never_proved_zero(power):
    r = S.sqrt(5)
    numerator = S.sqrt(5-r)*S.sqrt(40-8*r)-2*S.sqrt(2)*(5-r)
    assert numerator.is_zero is True
    tiny = numerator+S.Rational(1, 10**power)
    assert S.radsimp(S.Mul(S.zoo, tiny, evaluate=False)) is not S.nan


def test_smaller_exact_radical_zero_product():
    r = S.sqrt(5)
    zero = S.sqrt(5-r)*S.sqrt(40-8*r)-2*S.sqrt(2)*(5-r)
    assert S.radsimp(S.Mul(S.zoo, zero, evaluate=False)) is S.nan
