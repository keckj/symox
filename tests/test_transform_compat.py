"""Transformations used by Trisplyne, compared with SymPy as a test oracle."""

import pytest

import symox as sx

sp = pytest.importorskip("sympy")


@pytest.mark.parametrize("expression,old", [
    ("x + y + z", "x + y"),
    ("x - y + z", "y - x"),
    ("x + y", "x + y + 1"),
    ("x + y + 2", "x + y - 1"),
    ("2*x + 2*y + z", "x + y"),
    ("2*x + 3*y", "x + y"),
    ("(x + y + z)**2", "x + y"),
    ("x**2*y", "x**2"),
    ("x**4*y", "x**2"),
    ("2*x + y", "2*x"),
    ("2*x*y", "2*x"),
    ("x**2*y*z", "x*y"),
    ("x**3*y**2*z", "x*y"),
    ("8*x**3*y**3", "2*x*y"),
    ("6*x**3*y**2", "2*x*y"),
    ("x**2*y**3", "2*x*y"),
    ("3*x**3*y**3/2", "x*y/2"),
    ("3*x**3*y**3/2", "x*y/4"),
    ("2*x**3*y**3", "x*y/2"),
    ("-x*y*z", "-x*y"),
    ("x*y*z", "-x*y"),
    ("z/(x*y)", "x*y"),
    ("z/(x**3*y**2)", "x*y"),
    ("x*z/y", "x/y"),
    ("x**2*y + x*y**2 + x*y", "x*y"),
    ("x**4", "x**2"),
    ("x**3", "x**2"),
    ("x**(-2)", "x**2"),
])
def test_algebraic_substitution_matches_sympy(expression, old):
    expected = sp.sympify(expression).subs(sp.sympify(old), sp.Symbol("w"))
    actual = sx.sympify(expression).subs(sx.sympify(old), sx.Symbol("w"))
    assert actual == sx.sympify(str(expected))


@pytest.mark.parametrize("expression,divisors,gens", [
    ("x**2 + y**2", ["x + y"], ()),
    ("x**2 + y**2", ["x + y"], ("y", "x")),
    ("x*y**2 + x*y + x", ["x*y - 1", "y**2 - 1"], ()),
    ("x*y**2 + x*y + x", ["y**2 - 1", "x*y - 1"], ()),
    ("2*x + 1", ["3*x"], ()),
    ("x", ["a*x"], ("x",)),
    ("x*y + y", ["y"], ("x",)),
    ("edge1**2 + edge1*edge2 + edge3", ["edge1 + edge2"], ()),
    ("x + y", [], ()),
    ("0", ["x + y"], ()),
])
def test_reduced_matches_sympy(expression, divisors, gens):
    expected_q, expected_r = sp.reduced(
        sp.sympify(expression), [sp.sympify(d) for d in divisors],
        *(sp.Symbol(g) for g in gens),
    )
    actual_q, actual_r = sx.reduced(
        sx.sympify(expression), [sx.sympify(d) for d in divisors],
        *(sx.Symbol(g) for g in gens),
    )
    assert isinstance(actual_q, list)
    assert actual_q == [sx.sympify(str(q)) for q in expected_q]
    assert actual_r == sx.sympify(str(expected_r))
    # Certify the division identity independently of quotient formatting.
    original = sx.sympify(expression)
    reconstructed = sum((q * sx.sympify(d) for q, d in zip(actual_q, divisors)), actual_r)
    assert sx.expand(reconstructed - original) == 0


def test_reduced_rejects_zero_divisor():
    with pytest.raises(ZeroDivisionError):
        sx.reduced(sx.Symbol("x"), [sx.Integer(0)])


def test_edge_substitution_inside_barycentric_numerator():
    v1x, v2x, v1y, v2y, edge = sx.symbols("v1x v2x v1y v2y E")
    expr = v2x - v1x + v2y - v1y
    assert expr.subs({v2x - v1x: edge}) == edge + v2y - v1y


@pytest.mark.parametrize("options", [
    {"domain": "QQ"}, {"modulus": 2}, {"polys": True}, {"auto": False},
])
def test_reduced_rejects_unsupported_options(options):
    x = sx.Symbol("x")
    with pytest.raises(NotImplementedError):
        sx.reduced(x, [x], **options)
