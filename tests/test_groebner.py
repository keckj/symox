"""Native Groebner bases: oracle comparisons and algebraic certificates."""

import itertools
import random
import subprocess
import sys

import pytest
import sympy as sp

import symox as sm


ORDERS = ("lex", "grlex", "grevlex")
DOMAINS = ({}, {"domain": "QQ"}, {"modulus": 2}, {"modulus": 7})


def oracle(expr, names="x y z"):
    return sp.sympify(str(expr), locals=dict(zip(names.split(), sp.symbols(names))))


def assert_same_basis(expressions, names="x y", **options):
    sx = sm.symbols(names)
    px = sp.symbols(names)
    if not isinstance(sx, tuple):
        sx, px = (sx,), (px,)
    native = sm.groebner([sm.sympify(f) for f in expressions], *sx, **options)
    expected = sp.groebner([sp.sympify(f) for f in expressions], *px, **options)
    assert len(native) == len(expected)
    assert all(sp.expand(oracle(a) - b) == 0 for a, b in zip(native, expected))
    assert native.gens == sx
    assert native.order == str(expected.order)
    assert native.domain == str(expected.domain)
    assert native.is_zero_dimensional == expected.is_zero_dimensional
    return native, expected


@pytest.mark.parametrize("order", ORDERS)
@pytest.mark.parametrize("options", DOMAINS)
@pytest.mark.parametrize("system", [
    ["x*y - 2*y", "2*y**2 - x**2"],
    ["2*x-y", "3*y-1"],
    ["x**2+y**2-1", "x*y-2"],
    ["x-y**2"],  # graded leading term has opposite sign to the lex one
])
def test_reduced_basis_matches_sympy(system, order, options):
    assert_same_basis(system, order=order, **options)


@pytest.mark.parametrize("order", ORDERS)
@pytest.mark.parametrize("options", [{"domain": "QQ"}, {"modulus": 5}])
def test_cyclic_three_basis_and_certificates(order, options):
    expressions = ["x+y+z", "x*y+y*z+z*x", "x*y*z-1"]
    native, expected = assert_same_basis(expressions, names="x y z", order=order, **options)
    x, y, z = sm.symbols("x y z")
    basis = native.exprs
    # Original generators must reduce to zero, not merely agree with an oracle.
    for f in expressions:
        assert native.contains(sm.sympify(f))
        assert native.reduce(sm.sympify(f))[1] == 0
    # Every S-pair reduces to zero: Buchberger's independent certificate.
    symbols = sp.symbols("x y z")
    polys = [sp.Poly(oracle(e), *symbols, **{k: v for k, v in options.items() if k != "order"}) for e in basis]
    for a, b in itertools.combinations(polys, 2):
        ea, ca = a.terms(order=order)[0]
        eb, cb = b.terms(order=order)[0]
        lcm = tuple(max(i, j) for i, j in zip(ea, eb))
        monomial = lambda e: sp.prod(v**k for v, k in zip(symbols, e))
        if "modulus" in options:
            p = options["modulus"]
            ia, ib = sp.invert(int(ca), p), sp.invert(int(cb), p)
        else:
            ia, ib = 1/ca, 1/cb
        s = sp.expand(ia * monomial(tuple(i-j for i, j in zip(lcm, ea))) * a.as_expr()
                      - ib * monomial(tuple(i-j for i, j in zip(lcm, eb))) * b.as_expr())
        assert native.reduce(sm.sympify(str(s)))[1] == 0
    f = x**5*y + y**4 + z + 2
    q, r = native.reduce(f)
    delta = sp.expand(oracle(f) - sum(oracle(a)*oracle(b) for a, b in zip(q, basis)) - oracle(r))
    if "modulus" in options:
        assert sp.Poly(delta, *symbols, modulus=options["modulus"]).is_zero
    else:
        assert delta == 0
    assert sp.expand(oracle(r) - expected.reduce(oracle(f))[1]) == 0


def test_reduction_integer_coefficients_and_poly_return_mode():
    x, y = sm.symbols("x y")
    basis = sm.groebner([2*x-y], x, y)
    assert basis.exprs == [2*x-y]
    q, r = basis.reduce(x)
    assert q == [sm.Rational(1, 2)]
    assert r == y/2
    assert basis.reduce(x, auto=False) == ([0], x)
    assert basis.reduce(0) == ([], 0)
    for auto in (True, False):
        with pytest.raises(ValueError, match="integers"):
            basis.reduce(x/3, auto=auto)
    with pytest.raises(ValueError, match="integers"):
        basis.contains(x/3)
    assert basis.contains(2*x-y)
    assert not basis.contains(x)
    g = sm.GroebnerBasis([2*x-y], x, y, polys=True)
    assert isinstance(g[0], sm.Poly)
    assert isinstance(g.polys[0], sm.Poly)
    q, r = g.reduce(x)
    assert q[0].domain == "QQ" and r.domain == "QQ"
    assert q[0].as_expr() == sm.Rational(1, 2)
    assert r.as_expr() == y/2
    assert isinstance(sm.groebner([sm.Poly(x+y, x, y)])[0], sm.Poly)
    assert sm.groebner([sm.Poly(x, x, domain="QQ")]).domain == "ZZ"
    assert sm.groebner([sm.Poly(2*x+1, x, modulus=5)]).domain == "ZZ"


def test_empty_unit_and_positive_dimensional_ideals():
    x, y = sm.symbols("x y")
    for expressions in ([], [0], [0, 0]):
        g = sm.groebner(expressions, x, y)
        assert list(g) == []
        assert g.reduce(x+y) == ([], x+y)
        assert g.contains(0) and not g.contains(x)
        assert not g.is_zero_dimensional
    for expressions in ([1], [7], [x, x*y-1]):
        g = sm.groebner(expressions, x, y)
        assert list(g) == [1]
        assert g.reduce(x+y) == ([x+y], 0)
        assert g.contains(x+y)
        assert not g.is_zero_dimensional
    assert not sm.is_zero_dimensional([x*y-1], x, y)
    assert sm.is_zero_dimensional([x*x-1, y*y-1], x, y)


def test_inference_assumptions_indexing_and_determinism():
    x = sm.Symbol("x", positive=True)
    y = sm.Symbol("y", real=True)
    f = [x*y-1, y*y-x]
    g = sm.groebner(f)
    assert g.gens == (x, y)
    assert g[0] == g.exprs[0] and g[-1] == g.exprs[-1]
    assert g[:] == g.exprs
    with pytest.raises(IndexError):
        _ = g[len(g)]
    assert "GroebnerBasis" in repr(g) and "order='lex'" in repr(g)
    for permutation in (f[::-1], f+[f[0]], [0]+f):
        assert sm.groebner(permutation).exprs == g.exprs
    reversed_gens = sm.groebner(f, y, x)
    assert reversed_gens.gens == (y, x)
    assert all(reversed_gens.contains(e) for e in f)


@pytest.mark.parametrize("order", ORDERS)
def test_seeded_quadratic_systems(order):
    rng = random.Random(19485)
    x, y = sp.symbols("x y")
    monomials = [x*x, x*y, y*y, x, y, 1]
    for _ in range(8):
        fs = [str(sp.Add(*(rng.randint(-2, 2)*m for m in monomials))) for _ in range(2)]
        assert_same_basis(fs, order=order, domain="QQ")


def test_unsupported_options_and_invalid_inputs_are_explicit():
    x, y = sm.symbols("x y")
    for fs in ([], [0], [1]):
        with pytest.raises(ValueError, match="generators"):
            sm.groebner(fs)
    with pytest.raises(ValueError, match="distinct"):
        sm.groebner([x], x, x)
    with pytest.raises(ValueError, match="nonnumeric"):
        sm.groebner([x], 1)
    with pytest.raises(ValueError, match="order"):
        sm.groebner([x], x, order="ilex")
    with pytest.raises(ValueError, match="buchberger"):
        sm.groebner([x], x, method="f4")
    with pytest.raises(NotImplementedError, match="exact rational"):
        sm.groebner([1.5*x], x)
    with pytest.raises(NotImplementedError, match="extensions"):
        sm.groebner([sm.sqrt(2)*x], x)
    with pytest.raises(NotImplementedError, match="ZZ, QQ"):
        sm.groebner([x], x, domain="EX")
    with pytest.raises(ValueError, match="integers"):
        sm.groebner([x/2], x, domain="ZZ")
    with pytest.raises(ValueError, match="prime"):
        sm.groebner([x], x, modulus=6)
    with pytest.raises(ValueError, match="denominator"):
        sm.groebner([x/5], x, modulus=5)
    g = sm.groebner([x], x, modulus=5)
    with pytest.raises(ValueError, match="denominator"):
        g.reduce(x/5)
    with pytest.raises(ValueError, match="coefficients"):
        g.reduce(x+y)
    assert list(sm.groebner([x], x, extension=True)) == [x]


def test_groebner_does_not_import_sympy_at_runtime():
    subprocess.run([sys.executable, "-c", """
import sys
sys.modules['sympy'] = None
import symox as sm
x, y = sm.symbols('x y')
g = sm.groebner([x*y-1, y*y-x], x, y)
assert g.contains(x*y-1)
assert g.is_zero_dimensional
"""], check=True, timeout=30)


@pytest.mark.parametrize("order", ORDERS)
def test_reduced_honors_monomial_order(order):
    x, y, z = sm.symbols("x y z")
    f = x**3*y + y**4 + z**3 + x
    divisors = [x-y*y, y-z*z]
    q, r = sm.reduced(f, divisors, x, y, z, order=order)
    oq, oremainder = sp.reduced(oracle(f), [oracle(g) for g in divisors], *sp.symbols("x y z"), order=order)
    assert [oracle(e) for e in q] == oq
    assert oracle(r) == oremainder
    assert sm.expand(sum((a*b for a, b in zip(q, divisors)), r)-f) == 0
    assert sm.reduced(0, divisors, x, y, z, order=order) == ([], 0)
    with pytest.raises(ZeroDivisionError):
        sm.reduced(f, [0], x, y, z, order=order)
