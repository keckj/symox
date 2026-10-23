"""Regression tests for the simp-rest2 parity work (expected values from SymPy 1.14)."""

from symox import *

del test

x, y, c = symbols("x y c")


def _real(lhs, rhs, domain=S.Reals):
    return solveset(Eq(lhs, rhs), x, domain)


def test_negative_base_real_inverter_matches_sympy():
    assert _real((-2) ** x, 4) == FiniteSet(2)
    assert _real((-5) ** (x + 1), 625) == FiniteSet(3)
    assert _real((-3) ** (x - 3), -(3**39)) == FiniteSet(42)
    assert _real((-2) ** x, -8) == FiniteSet(3)
    assert _real((-2) ** x, -2) == FiniteSet(1)
    assert _real((-2) ** (2 * x), 16) == FiniteSet(2)
    assert _real((-2) ** (x**2), 16) == FiniteSet(-2, 2)
    assert _real((-2) ** x + 1, 5) == FiniteSet(2)
    assert _real(3 * (-2) ** x, 12) == FiniteSet(2)
    assert _real((-2) ** (2 * x + 1), -8) == FiniteSet(1)
    assert _real((-2) ** (x / 3), 4) == FiniteSet(6)
    assert _real((-2) ** x, 4, Interval(0, 10)) == FiniteSet(2)
    assert _real((-3) ** x, 9, Interval(-5, 5)) == FiniteSet(2)
    assert _real((-3) ** x, 9, S.Naturals) == FiniteSet(2)
    assert _real((-2) ** x - 1, 0) == FiniteSet(0)
    assert solveset(exp((-2) ** x) - exp(4), x, S.Reals) == FiniteSet(2)
    assert solveset(log((-2) ** x) - log(4), x, S.Reals) == FiniteSet(2)


def test_negative_base_no_real_solution():
    # (-2)**x is real only on the integers, where it is (-1)**x*2**x.
    for lhs, rhs in [((-2) ** x, 2), ((-2) ** x, 8), ((-2) ** x, 3), ((-2) ** x, 0),
                     ((-4) ** x, 2), ((-8) ** x, 2), ((-8) ** x, -2), ((-2) ** x, -4),
                     ((-2) ** (x / 2), 2)]:
        assert _real(lhs, rhs) == S.EmptySet
    assert _real((-2) ** x, 4, Interval(3, 10)) == S.EmptySet


def test_negative_base_undecided_stays_conditionset():
    assert _real((-2) ** x, 4, S.Complexes) == ConditionSet(x, Eq((-2) ** x - 4, 0), S.Complexes)
    assert _real((-2) ** x * 2**x, 16) == ConditionSet(x, Eq((-2) ** x * 2**x - 16, 0), S.Reals)
    assert _real((-2) ** x + x, 0) == ConditionSet(x, Eq((-2) ** x + x, 0), S.Reals)
    assert isinstance(_real((-2) ** x, y), ConditionSet)
    assert isinstance(_real((-2) ** x, sqrt(2) * I), ConditionSet)
    assert isinstance(_real((-c) ** x, 4), ConditionSet)
    cp = Symbol("cp", positive=True)
    assert isinstance(_real((-cp) ** x, 4), ConditionSet)
    r = Symbol("r", real=True)
    # SymPy inverts only a positive or negative base; the principal
    # logarithm misses x = 2 for c = -2.
    for base in [c, r, r**2, I, 1 + I, -r**2 - 1]:
        lhs = base**x
        assert _real(lhs, 4) == ConditionSet(x, Eq(lhs - 4, 0), S.Reals)
    assert _real(c**x * (x - 1), 0) == ConditionSet(x, Eq(c**x * (x - 1), 0), S.Reals)
    assert _real((x - 1) * (c**x - 4), 0) == ConditionSet(x, Eq((c**x - 4) * (x - 1), 0), S.Reals)
    assert str(_real(c * exp(x), 1)) == "Intersection({log(1/c)}, Reals)"
    assert isinstance(_real((-2) ** x, Float(4)), ConditionSet)
    assert isinstance(_real((-2) ** (2 * x) + (-2) ** x, 20), ConditionSet)
    # SymPy returns EmptySet, but (-2)**x == 2**x holds on every even integer.
    assert isinstance(_real((-2) ** x, 2**x), ConditionSet)


def test_negative_base_intentional_corrections():
    # SymPy raises ValueError from integer_log for these; each value checks.
    assert _real((-2) ** x, Rational(1, 4)) == FiniteSet(-2)
    assert _real(Rational(-1, 2) ** x, 4) == FiniteSet(-2)
    assert _real((-2) ** x, Rational(-1, 8)) == FiniteSet(-3)
    assert _real(Rational(-2, 3) ** x, Rational(81, 16)) == FiniteSet(-4)
    assert _real((-E) ** x, exp(2)) == FiniteSet(2)
    assert _real((-pi) ** x, pi**2) == FiniteSet(2)
    assert _real((-sqrt(2)) ** x, 4) == FiniteSet(4)
    assert _real((-sqrt(2)) ** x, -2 * sqrt(2)) == FiniteSet(3)
    assert _real((-2) ** x, sqrt(2)) == S.EmptySet
    n = Dummy("n")
    assert _real((-1) ** x, 1).dummy_eq(ImageSet(Lambda(n, 2 * n), S.Integers))
    assert _real((-1) ** x, -1).dummy_eq(ImageSet(Lambda(n, 2 * n + 1), S.Integers))
    assert _real((-1) ** (3 * x + 1), -1).dummy_eq(ImageSet(Lambda(n, 2 * n / 3), S.Integers))
    assert _real((-1) ** x, 1, Interval(0, 7)) == FiniteSet(0, 2, 4, 6)
    assert _real((-1) ** x, 2) == S.EmptySet
    assert _real((-2) ** x * (x - 1), 0) == FiniteSet(1)


def test_negative_base_sin_exponent():
    n = Dummy("n")
    assert _real((-2) ** sin(x), -2).dummy_eq(ImageSet(Lambda(n, 2 * n * pi + pi / 2), S.Integers))


def test_solve_negative_base_integer_guess():
    xr = Symbol("xr", real=True)
    assert solve((-2) ** x - 4, x) == [2, log(4) / (log(2) + I * pi)]
    assert solve((-2) ** x + 8, x) == [3, (log(8) + I * pi) / (log(2) + I * pi)]
    # SymPy also keeps a non-real second root for the real symbol
    assert solve((-2) ** xr - 4, xr) == [2]


def test_nthroot():
    from symox.simplify.simplify import nthroot
    assert nthroot(90 + 34 * sqrt(7), 3) == sqrt(7) + 3
    assert nthroot(41 + 29 * sqrt(2), 5) == 1 + sqrt(2)
    assert nthroot(-41 - 29 * sqrt(2), 5) == -1 - sqrt(2)
    q = 1 + sqrt(2) - 2 * sqrt(3) + sqrt(6) + sqrt(7)
    assert nthroot(expand_multinomial(q**3), 3) == q
    assert nthroot(expand_multinomial(q**5), 5, 8) == q


def test_sqrtdenest_numeric_helpers():
    from symox.simplify.sqrtdenest import _sqrt_numeric_denest, _subsets
    assert _subsets(1) == [[1]]
    assert _subsets(2) == [[1, 0], [0, 1], [1, 1]]
    assert _subsets(3) == [[1, 0, 0], [0, 1, 0], [1, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]]
    # issue 19914
    a = Integer(-8)
    b = Integer(-1)
    r = Integer(63)
    d2 = a * a - b * b * r
    assert _sqrt_numeric_denest(a, b, r, d2) == sqrt(14) * I / 2 + 3 * sqrt(2) * I / 2
    assert sqrtdenest(sqrt(-8 - sqrt(63))) == sqrt(14) * I / 2 + 3 * sqrt(2) * I / 2


def test_minimal_polynomial_sq():
    from symox.polys.numberfields.minpoly import _minimal_polynomial_sq, _is_sum_surds, _separate_sq
    assert _minimal_polynomial_sq(1 + sqrt(2) + sqrt(3), 1, x) == x**4 - 4 * x**3 - 4 * x**2 + 16 * x - 8
    assert _is_sum_surds(1 + sqrt(2) + sqrt(3))
    assert not _is_sum_surds(1 + sqrt(2) + pi)
    p = -x + sqrt(2) + sqrt(3) + sqrt(7)
    assert _separate_sq(p) == -x**2 + 2 * sqrt(3) * x + 2 * sqrt(7) * x - 2 * sqrt(21) - 8
