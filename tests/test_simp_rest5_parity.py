"""Wave 11 simplification leftovers (expected values from SymPy 1.14)."""
import pytest

from symox import *

del test

x, y = symbols("x y")
i = symbols("i", imaginary=True)


def test_simplify_reference_sequence_forms():
    assert str(simplify((x + 5) / (16 * x**2 + 28 * x + 6))) == "(x + 5)/(2*(8*x**2 + 14*x + 3))"
    assert simplify(x / (x * y) + y / (x * y)) == (x + y) / (x * y)
    assert str(simplify(exp(-Rational(1, 2)) + exp(-Rational(3, 2)))) == "(1 + E)*exp(-3/2)"


def test_simplify_piecewise_form():
    assert str(simplify(Piecewise((x**2 - y**2, x < 1), (y * (x + y) / x, True)))) == (
        "Piecewise((x**2 - y**2, x < 1), (y*(x + y)/x, True))"
    )


def test_issue_17141():
    assert simplify(x ** (1 / acos(I))) == x ** (2 / (pi - 2 * I * log(1 + sqrt(2))))
    assert simplify(acos(I) ** 2) == (pi - 2 * I * log(1 + sqrt(2))) ** 2 / 4
    assert simplify(acos(-I) ** 2 * acos(I) ** 2) == (
        log(1 + sqrt(2)) ** 4 + pi**2 * log(1 + sqrt(2)) ** 2 / 2 + pi**4 / 16
    )
    assert simplify(2 ** acos(I) ** 2) == 2 ** ((pi - 2 * I * log(1 + sqrt(2))) ** 2 / 4)
    p = 2 ** acos(I + 1) ** 2
    assert simplify(p) == p


def _facts(e):
    return [getattr(e, "is_" + name) for name in ("real", "positive", "nonnegative", "imaginary", "zero", "nonzero")]


def test_imaginary_powers_of_imaginary_bases():
    # Pow._eval_is_extended_real/positive: log(I) is imaginary. SymPy's
    # positive test starts with `base == exp`, so I**I keeps an unknown sign.
    assert _facts(I**I) == [True, None, None, False, False, True]
    for e in ((-I) ** I, I ** (2 * I), I**i, (-I) ** i, exp(I) ** i):
        assert _facts(e) == [True, True, True, False, False, True], e
    assert (I**i).is_nonnegative is True
    assert (exp(I) ** i).is_nonnegative is True
    assert (I**I).rewrite(exp) == exp(-pi / 2)
    assert Abs(I**I) == exp(-pi / 2)


def test_products_of_imaginary_factors():
    assert _facts(I * i) == [True, None, None, False, False, True]
    assert (2 * I * i).is_real is True
    assert (2 * I * i).is_positive is None
    assert (I * i * x).is_real is None
    assert (I * i * symbols("r", real=True)).is_real is True


def test_held_two_argument_log():
    e = log(exp(x), 2, evaluate=False)
    assert str(e) == "log(exp(x), 2)"
    assert e.args == (exp(x), 2)
    assert srepr(e) == "log(exp(Symbol('x')), Integer(2))"
    assert e.subs(x, y) == log(exp(y)) / log(2)
    assert e.subs(x, 2) == 2 / log(2)
    assert e.doit() == log(exp(x)) / log(2)
    assert expand(e) == log(exp(x)) / log(2)
    assert e.func(*e.args) == log(exp(x)) / log(2)
    assert simplify(e) == log(exp(x)) / log(2)
    assert simplify(e, inverse=True) == x / log(2)


def test_matrix_doit_inv_expand():
    A = MatrixSymbol("A", 3, 3)
    B = MatrixSymbol("B", 3, 3)
    assert str(Inverse(A * B).doit(inv_expand=False)) == "(A*B)**(-1)"
    assert str(Inverse(A * B).doit()) == "B**(-1)*A**(-1)"
    assert str((A**-1 * B).doit(inv_expand=False)) == "A**(-1)*B"
    assert str((A * B).I.doit(inv_expand=False)) == "B**(-1)*A**(-1)"
