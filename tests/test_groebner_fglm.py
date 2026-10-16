"""FGLM conversion and polynomial division across native coefficient domains."""

import itertools
import random
import subprocess
import sys

import pytest
import sympy as sp
import symox as sm


ORDERS = ("lex", "grlex", "grevlex")


def expr(value):
    return value.as_expr() if isinstance(value, sm.Poly) else value


def oracle(value):
    return sp.sympify(str(expr(value)))


@pytest.mark.parametrize("source,target", list(itertools.product(ORDERS, ORDERS)))
@pytest.mark.parametrize("options", [{}, {"domain": "QQ"}, {"modulus": 2}, {"modulus": 5}])
@pytest.mark.parametrize("system", [
    ["x*y-1", "y**2-x"],
    ["x**2-3*y-x+1", "y**2-2*x+y-1"],
    ["x**2", "y**3"],  # nonradical quotient with nilpotent multiplication matrices
    ["x**2+y**2-1", "x*y-2"],
])
def test_fglm_matches_sympy_and_direct_basis(system, source, target, options):
    x, y = sm.symbols("x y")
    sx, sy = sp.symbols("x y")
    f = [sm.sympify(e) for e in system]
    g = sm.groebner(f, x, y, order=source, **options)
    expected = sp.groebner([sp.sympify(e) for e in system], sx, sy, order=source, **options).fglm(target)
    converted = g.fglm(target)
    assert converted.order == target
    assert converted.gens == (x, y)
    assert converted.domain == str(expected.domain)
    assert [oracle(e) for e in converted] == list(expected)
    assert converted.exprs == sm.groebner(f, x, y, order=target, **options).exprs
    assert converted.fglm(source).exprs == g.exprs
    for e in converted:
        assert g.contains(e)
    for e in g:
        assert converted.contains(e)


def test_fglm_three_variables_and_poly_return_mode():
    x, y, z = sm.symbols("x y z")
    f = [x+y+z, x*y+y*z+z*x, x*y*z-1]
    g = sm.groebner(f, x, y, z, order="grevlex", polys=True)
    lex = g.fglm("lex")
    assert isinstance(lex[0], sm.Poly)
    assert lex.exprs == sm.groebner(f, x, y, z).exprs
    assert lex.fglm("grevlex").exprs == g.exprs
    assert g.order == "grevlex"  # conversion preserves the source object
    assert lex.is_zero_dimensional


@pytest.mark.parametrize("system", [[], [0], [1], ["x*y-1"]])
def test_fglm_positive_dimensional_ideals_are_explicit(system):
    x, y = sm.symbols("x y")
    g = sm.groebner([sm.sympify(e) for e in system], x, y, order="grlex")
    assert g.fglm("grlex").exprs == g.exprs
    with pytest.raises(NotImplementedError, match="positive dimension"):
        g.fglm("lex")
    with pytest.raises(ValueError, match="order"):
        g.fglm("bogus")


@pytest.mark.parametrize("order", ORDERS)
def test_fglm_seeded_zero_dimensional_quotients(order):
    rng = random.Random(17567)
    x, y = sm.symbols("x y")
    for _ in range(8):
        # Monic triangular source systems guarantee a finite quotient algebra.
        f = [x*x + rng.randint(-3, 3)*y + rng.randint(-2, 2), y**3 + rng.randint(-3, 3)*y + rng.randint(-2, 2)]
        g = sm.groebner(f, x, y, order=order, domain="QQ")
        for target in ORDERS:
            assert g.fglm(target).exprs == sm.groebner(f, x, y, order=target, domain="QQ").exprs


@pytest.mark.parametrize("order", ORDERS)
@pytest.mark.parametrize("options", [
    {}, {"domain": "QQ"}, {"domain": "ZZ"}, {"auto": False},
    {"domain": "ZZ", "auto": True}, {"domain": "QQ", "auto": False},
    {"modulus": 5}, {"modulus": 2}, {"polys": True},
    {"polys": True, "domain": "QQ"}, {"polys": True, "modulus": 5},
    {"polys": True, "auto": False},
])
def test_reduced_options_match_sympy(order, options):
    x, y = sm.symbols("x y")
    sx, sy = sp.symbols("x y")
    f, divisors = x**3*y + y**4 + 1, [3*x-y*y, y*y-1]
    expected = sp.reduced(oracle(f), [oracle(e) for e in divisors], sx, sy, order=order, **options)
    q, r = sm.reduced(f, divisors, x, y, order=order, **options)
    assert [oracle(e) for e in q] == [e.as_expr() if isinstance(e, sp.Poly) else e for e in expected[0]]
    assert oracle(r) == (expected[1].as_expr() if isinstance(expected[1], sp.Poly) else expected[1])
    if options.get("polys"):
        assert all(isinstance(e, sm.Poly) for e in q+[r])
        assert [e.domain for e in q+[r]] == [str(e.domain) for e in expected[0]+[expected[1]]]
    delta = sp.expand(oracle(f)-sum(oracle(a)*oracle(b) for a, b in zip(q, divisors))-oracle(r))
    if "modulus" in options:
        assert sp.Poly(delta, sx, sy, modulus=options["modulus"]).is_zero
    else:
        assert delta == 0


def test_reduced_domain_selection_and_coercion():
    x = sm.Symbol("x")
    assert sm.reduced(x, [2*x], x) == ([sm.Rational(1, 2)], 0)
    assert sm.reduced(x, [2*x], x, domain="ZZ") == ([0], x)
    assert sm.reduced(x, [2*x], x, auto=False) == ([0], x)
    assert sm.reduced(x, [2*x], x, modulus=5) == ([-2], 0)
    for options in ({"domain": "ZZ"}, {"modulus": 5}):
        with pytest.raises(ValueError, match="integers"):
            sm.reduced(x/2, [x], x, **options)
    q, r = sm.reduced(sm.Poly(x, x, modulus=5), [sm.Poly(2*x, x, modulus=5)])
    assert q[0].as_expr() == -2 and q[0].domain == "GF(5)" and r.domain == "GF(5)"
    for options in ({}, {"polys": True}, {"modulus": 5}, {"auto": False}):
        q, r = sm.reduced(0, [x], x, **options)
        assert q == [] and expr(r) == 0
    with pytest.raises(ZeroDivisionError):
        sm.reduced(x, [5*x], x, modulus=5)
    with pytest.raises(ValueError, match="generators"):
        sm.reduced(1, [])


def test_fglm_without_sympy_runtime():
    subprocess.run([sys.executable, "-c", """
import sys
sys.modules['sympy'] = None
import symox as s
x, y = s.symbols('x y')
g = s.groebner([x*y-1,y*y-x],x,y,order='grevlex')
h = g.fglm('lex')
assert h.exprs == [x-y*y,y**3-1]
assert s.reduced(x,[2*x],x,domain='QQ') == ([s.Rational(1,2)],0)
"""], check=True, timeout=30)
