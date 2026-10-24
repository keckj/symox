"""Resource limits: an exhausted budget never proves the absence of a result.

Expected values are hard-coded from SymPy 1.14.  Each case used to hit a
symox budget that returned a wrong or partial answer (an input reported
irreducible, an unfactored product, a series term reported nonzero, an
inexact integer part, ...).
"""
import random

import mpmath
import pytest

from symox import *
del test


x, y = symbols('x y')


def _random_poly(n, bound):
    return sum(random.randint(-bound, bound)*x**i for i in range(n)) + x**n


def test_trinomial_irreducible_without_work_exhaustion():
    assert factor_list(x**199 + x + 1) == (1, [(x**199 + x + 1, 1)])


def test_products_of_random_dense_factors_are_split():
    random.seed(3)
    p, q = _random_poly(60, 1000), _random_poly(60, 1000)
    assert sorted(Poly(f, x).degree() for f, _ in factor_list(expand(p*q))[1]) == [60, 60]
    random.seed(4)
    ps = [_random_poly(12, 10**60) for _ in range(4)]
    assert sorted(Poly(f, x).degree() for f, _ in factor_list(expand(Mul(*ps)))[1]) == [12]*4


def test_ten_variable_factorization():
    v = symbols('v0:10')
    a, b = sum(v) + 1, v[0] - v[1] + v[9]
    coeff, factors = factor_list(expand(a*b))
    assert coeff == 1 and set(factors) == {(a, 1), (b, 1)}
    coeff, factors = sqf_list(expand(a**2*b))
    assert coeff == 1 and set(factors) == {(a, 2), (b, 1)}


def test_bivariate_evaluation_points_are_searched_until_found():
    # Every y in -15..15 is a root of the leading coefficient: SymPy keeps
    # widening the evaluation range; symox used to report F irreducible.
    lead = Mul(*[y - k for k in range(-15, 16)])
    product = (lead*x + 1)*(x + y)
    expanded = expand(product)
    assert Poly(expanded, x, y).is_irreducible is False
    assert len(factor_list(expanded)[1]) == 2
    assert factor(expanded) == factor(product)


def test_large_hidden_zero_in_series_is_detected():
    a, b, c, d = symbols('a b c d')
    big = (a + b + c + d)**8
    g = big*sin(x) - expand(big)*sin(x) + x**2
    assert g.as_leading_term(x) == x**2
    assert limit(g/x**2, x, 0) == 1


def test_deep_array_kind():
    e = x
    for _ in range(70):
        e = (e + 1)**2
    assert str(Array([e]).kind) == 'ArrayKind(NumberKind)'


def test_decompogen_never_returns_a_truncated_chain():
    assert decompogen(x**4000 + x**2, x) == [x**2000 + x, x**2]
    expected = [x**2050 + x, x**2]
    try:
        result = decompogen(x**4100 + x**2, x)
    except NotImplementedError as error:
        assert 'decompogen' in str(error)
    else:
        assert result == expected


def test_eigenvalues_with_large_constant_terms():
    assert Matrix([[10**13, 0], [0, 1]]).eigenvals() == {10**13: 1, 1: 1}
    root = sqrt(99999999999940000000000017)/2
    assert Matrix([[10**13, 1], [2, 3]]).eigenvals() == {
        Rational(10000000000003, 2) - root: 1, Rational(10000000000003, 2) + root: 1}
    assert Matrix([[3*10**12 + 1, 0, 0], [0, 7, 1], [0, 1, 5]]).eigenvals() == {
        3000000000001: 1, 6 - sqrt(2): 1, 6 + sqrt(2): 1}


def test_checkodesol_high_order():
    f = Function('f')
    assert checkodesol(f(x).diff(x, 300) - f(x), Eq(f(x), exp(x))) == (True, 0)


@pytest.mark.parametrize('k', [300, 3000])
def test_mod_integer_part_is_exact(k):
    m = Mod(exp(k), 7)
    shifts = [arg for arg in m.args if arg.is_Integer]
    assert m.func is Add and len(shifts) == 1 and m - shifts[0] == exp(k)
    with mpmath.workdps(k + 100):
        exact = int(mpmath.floor(mpmath.exp(k)/7))*7
    assert int(-shifts[0]) == exact
    if k == 300:
        assert 0 <= m.evalf(30) < 7


def test_mod_of_unevaluated_huge_integer_is_never_inexact():
    # SymPy evaluates fibonacci(1000001) (residue 501 mod 1000); symox keeps
    # it unevaluated above its exact-evaluation cap, so Mod may stay
    # unevaluated, but it must never subtract an inexact multiple.
    m = Mod(fibonacci(1000001), 1000)
    assert m == 501 or m.func is Mod
