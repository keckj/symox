"""Regression tests for the simp-rest3 parity work (expected values from SymPy 1.14)."""

import pytest

from symox import *

del test

x, y = symbols("x y")
r = Symbol("r", real=True)
xr = Symbol("xr", real=True)


def test_imaginary_unit_power_is_finite_for_real_exponent():
    assert (I**xr).is_finite is True
    assert (I**x).is_finite is None
    assert solveset(I**xr * (xr - 2), xr, S.Reals) == FiniteSet(2)


def test_even_power_of_real_base_folds_to_abs():
    assert (r**2) ** x == Abs(r) ** (2 * x)
    assert (r**4) ** x == Abs(r) ** (4 * x)
    assert (r**2) ** Rational(1, 3) == Abs(r) ** Rational(2, 3)
    assert str((x**2) ** y) == "(x**2)**y"


def test_nthroot_of_near_integer_radical_sum():
    from symox.simplify.simplify import nthroot

    for k in (20, 30):
        q = 1 + sqrt(2) + sqrt(3) + S.One / 10**k
        p = expand_multinomial(q**5)
        assert nthroot(p, 5) == q


def test_simplify_keeps_factored_denominator():
    e = simplify((x + 5) / (16 * x**2 + 28 * x + 6))
    assert str(e) == "(x + 5)/(2*(8*x**2 + 14*x + 3))"


def test_nc_simplify_matrix_products():
    from symox.simplify.simplify import nc_simplify

    A, B, C = [MatrixSymbol(n, 3, 3) for n in "ABC"]
    assert nc_simplify(A * B * A * B * A * B * C * (A * B) ** 3 * C) == ((A * B) ** 3 * C) ** 2
    assert nc_simplify(A * B * (A * B) ** -1) == Identity(3)
    assert nc_simplify((A * B) ** 2 * A * B * C) == (A * B) ** 2 * A * B * C
    assert nc_simplify(A * B * A * B) == (A * B) ** 2
    assert nc_simplify(A * B * A * B + C) == (A * B) ** 2 + C
    assert nc_simplify(A * (B * A) ** -1 * A) == B**-1 * A


def test_trigsimp_groebner():
    from symox.simplify.trigsimp import trigsimp_groebner

    c, s = cos(x), sin(x)
    ex = (4 * s * c + 12 * s + 5 * c**3 + 21 * c**2 + 23 * c + 15) / (
        -s * c**2 + 2 * s * c + 15 * s + 7 * c**3 + 31 * c**2 + 37 * c + 21
    )
    resnum = 5 * s - 5 * c + 1
    resdenom = 8 * s - 6 * c
    results = [resnum / resdenom, (-resnum) / (-resdenom)]
    assert trigsimp_groebner(ex) in results
    assert trigsimp_groebner(s / c, hints=[tan]) == tan(x)
    assert trigsimp_groebner(c * s) == c * s
    assert trigsimp((-s + 1) / c + c / (-s + 1), method="groebner") == 2 / c
    assert trigsimp((-s + 1) / c + c / (-s + 1), method="groebner", polynomial=True) == 2 / c
    assert trigsimp_groebner(ex, hints=[2]) in results
    assert trigsimp_groebner(sin(I * x) / cos(I * x), hints=[tanh]) == I * tanh(x)
    assert trigsimp_groebner(
        (tanh(x) + tanh(y)) / (1 + tanh(x) * tanh(y)), hints=[(tanh, x, y)]
    ) == tanh(x + y)
    assert trigsimp(log(2 * sin(x) ** 2 + 2 * cos(x) ** 2), method="groebner") == log(2)
    assert trigsimp(sin(x) ** 2 + cos(x) ** 2 + tan(x), method="combined") == tan(x) + 1
    assert trigsimp_groebner(sinh(2 * x) / sinh(x)) == 2 * cosh(x)
    assert trigsimp_groebner(sin(x) ** 2 - cos(x) ** 2, hints=[2]) == -cos(2 * x)
    assert trigsimp(Matrix.eye(1), method="groebner") == Matrix.eye(1)


def test_powsimp_polar():
    p, q = symbols("p q", polar=True)
    assert powsimp(p**x * q**x) == (p * q) ** x
    assert powsimp(p**x * q**x, combine="base") == (p * q) ** x
    assert powsimp(exp_polar(x) * exp_polar(y)) == exp_polar(x + y)
    assert powsimp(exp_polar(x) * exp_polar(y) * p**x * p**y) == (p * exp_polar(1)) ** (x + y)
    assert powsimp(exp_polar(x) * exp_polar(y) * p**x * p**y, combine="exp") == exp_polar(
        x + y
    ) * p ** (x + y)
    assert powsimp(
        exp_polar(x) * exp_polar(y) * exp_polar(2) * sin(x) + sin(y) + p**x * p**y
    ) == p ** (x + y) + sin(x) * exp_polar(2 + x + y) + sin(y)
    assert powsimp(sin(exp_polar(x) * exp_polar(y))) == sin(exp_polar(x) * exp_polar(y))
    assert powsimp(sin(exp_polar(x) * exp_polar(y)), deep=True) == sin(exp_polar(x + y))
    assert powsimp(exp_polar(x) * exp_polar(-x)) == 1
    assert powsimp(p**x * (1 / p) ** y) == p ** (x - y)
    assert powsimp(exp_polar(1) * exp_polar(x)) == exp_polar(x + 1)


def test_expand_splits_integral_and_sum_with_finite_limits():
    assert expand(Integral(x + y, (y, 1, 3))) == Integral(x, (y, 1, 3)) + Integral(y, (y, 1, 3))
    assert expand(Integral(x + y, (y, 1, oo))) == Integral(x + y, (y, 1, oo))
    assert expand(Sum(x + y, (y, 1, 3))) == Sum(x, (y, 1, 3)) + Sum(y, (y, 1, 3))


def test_groebner_trigsimp_keeps_rational_angle_capability():
    # SymPy 1.14 raises ValueError("1/2 is not an integer") in
    # analyse_gens here; the identity is 0, which symox keeps deciding.
    beta = Symbol("beta", real=True)
    e = -6 * sin(beta / 2) ** 8 + 12 * sin(beta / 2) ** 6 - 6 * sin(beta / 2) ** 4 + 3 * sin(beta) ** 4 / 8
    assert trigsimp(e, method="groebner") == 0
