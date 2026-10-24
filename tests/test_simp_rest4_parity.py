"""Regression tests for the simp-rest4 parity work (expected values from SymPy 1.14)."""

from symox import *

del test

x, y = symbols("x y")
r = Symbol("r", real=True)
p = Symbol("p", positive=True)
ip = Symbol("ip", integer=True, positive=True)

A = Float("1.0000000000000000000001", 30)
B = Float("2.00000000000000000000000000003", 40)

FA = "Float('1.00000000000000000000009999999998', precision=103)"
FB = "Float('2.000000000000000000000000000029999999999991', precision=136)"
FA136 = "Float('1.000000000000000000000099999999981995874774', precision=136)"
F3 = "Float('3.000000000000000000000100000029981995874764', precision=136)"


def test_issue_11004_cancel_keeps_high_precision():
    assert srepr(cancel(A * x / x)) == FA
    assert srepr(cancel((A * x**2 - A) / (x - 1))) == f"Add(Mul({FA}, Symbol('x')), {FA})"
    assert srepr(simplify((A * x**2 - A) / (x - 1))) == f"Add(Mul({FA}, Symbol('x')), {FA})"
    assert srepr(cancel((x**2 - B) / (x + A))) == (
        "Mul(Pow(Add(Mul(Float('1.0', precision=136), Symbol('x')), "
        f"{FA136}), Integer(-1)), Add(Mul(Float('1.0', precision=136), "
        "Pow(Symbol('x'), Integer(2))), "
        "Float('-2.000000000000000000000000000029999999999991', precision=136)))"
    )
    # Machine precision is unchanged.
    assert srepr(cancel((Float(1.5) * x**2 - Float(1.5)) / (x - 1))) == (
        "Add(Mul(Float('1.5', precision=53), Symbol('x')), Float('1.5', precision=53))"
    )


def test_high_precision_through_polys_functions():
    assert srepr(together(A / x + B / y)) == (
        "Mul(Pow(Symbol('x'), Integer(-1)), Pow(Symbol('y'), Integer(-1)), "
        f"Add(Mul({FB}, Symbol('x')), Mul({FA}, Symbol('y'))))"
    )
    assert srepr(factor(A * x**2 - A)) == (
        f"Mul({FA}, Add(Symbol('x'), Integer(-1)), Add(Symbol('x'), Integer(1)))"
    )
    assert srepr(gcd(A * x**2 - A, x - 1)) == (
        "Add(Mul(Float('1.0', precision=103), Symbol('x')), Float('-1.0', precision=103))"
    )
    assert srepr(div(A * x**2 + B, x - 1)) == (
        f"(Add(Mul({FA136}, Symbol('x')), {FA136}), {F3})"
    )
    assert srepr(rem(B * x**2 + A, x - 1)) == F3
    assert srepr(simplify(A * x + B * x)) == f"Mul({F3}, Symbol('x'))"
    assert srepr(ratsimp(A / x + B / x)) == f"Mul({F3}, Pow(Symbol('x'), Integer(-1)))"
    assert srepr(lcm(A * x, x**2)) == "Mul(Float('1.0', precision=103), Pow(Symbol('x'), Integer(2)))"
    assert srepr(terms_gcd(A * x**2 + B * x)) == (
        f"Mul(Symbol('x'), Add(Mul({FA136}, Symbol('x')), {FB}))"
    )
    assert srepr(sqf(A * x**2 + 2 * A * x + A)) == (
        f"Mul({FA}, Pow(Add(Mul(Float('1.0', precision=103), Symbol('x')), "
        "Float('1.0', precision=103)), Integer(2)))"
    )


def test_poly_with_high_precision_coefficients():
    P = Poly(A * x**2 + B, x)
    assert str(P.domain) == "RR"
    assert [c._prec for c in P.coeffs()] == [136, 136]
    assert srepr(P.as_expr()) == (
        f"Add(Mul({FA136}, Pow(Symbol('x'), Integer(2))), {FB})"
    )
    R = RealField(136)
    Q = Poly(A * x**2 - A, x, domain=R)
    g = Q.gcd(Poly(x - 1, x, domain=R))
    assert srepr(g.as_expr()) == (
        "Add(Mul(Float('1.0', precision=136), Symbol('x')), Float('-1.0', precision=136))"
    )


def test_add_is_zero_with_single_imaginary_term():
    assert (r + I).is_zero is False
    assert (p + I).is_zero is False
    assert (1 + I * r).is_zero is False
    assert (x + I).is_zero is None


def test_pow_complex_base_zero_and_finite():
    # SymPy 1.14 values (an unassumed exponent may be infinite).
    assert ((1 + I) ** x).is_zero is None
    assert ((1 + I) ** r).is_zero is False
    assert ((1 + I) ** r).is_finite is True
    assert ((r + I) ** r).is_zero is False
    assert ((r + I) ** p).is_finite is True
    assert ((p + I) ** (I * r)).is_zero is False
    assert ((Rational(1, 2) + I / 2) ** ip).is_zero is False
    assert ((3 + 4 * I) ** x).is_finite is None
