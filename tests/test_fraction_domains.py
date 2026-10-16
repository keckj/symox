"""Exact rational-function coefficient domains; no symbolic fallback engine."""

import subprocess
import sys
import pickle

import pytest
import sympy as sp
import symox as sm


ORDERS = ("lex", "grlex", "grevlex")


def expression(p):
    return p.as_expr() if isinstance(p, sm.Poly) else p


def oracle(p):
    return sp.sympify(str(expression(p)))


def equal(a, b):
    return sp.cancel(oracle(a)-b) == 0


@pytest.mark.parametrize("order", ORDERS)
@pytest.mark.parametrize("domain", ["QQ(a,b)", "ZZ(a,b)", "QQ[a,b]", "ZZ[a,b]", None])
@pytest.mark.parametrize("system", [
    ["a*x-1", "y-x"],
    ["(a+1)*x+y-1", "x+(b+1)*y"],
    ["x*y-a", "y**2-b*x"],
])
def test_fraction_field_and_coefficient_ring_bases_match_sympy(system, order, domain):
    x, y = sm.symbols("x y")
    sx, sy = sp.symbols("x y")
    options = {"order": order}
    if domain is not None:
        options["domain"] = domain
    g = sm.groebner([sm.sympify(e) for e in system], x, y, **options)
    expected = sp.groebner([sp.sympify(e) for e in system], sx, sy, **options)
    assert str(g.domain) == str(expected.domain)
    assert len(g) == len(expected)
    assert all(equal(a, b) for a, b in zip(g, expected))
    assert g.is_zero_dimensional == expected.is_zero_dimensional
    for f in system:
        assert g.contains(sm.sympify(f))
        assert g.reduce(sm.sympify(f))[1] == 0
    assert not g.contains(sm.Integer(1))
    f = x**3*y + x+y+1
    q, r = g.reduce(f)
    oq, oremainder = expected.reduce(oracle(f))
    assert len(q) == len(oq)
    assert all(equal(a, b) for a, b in zip(q, oq))
    assert equal(r, oremainder)
    assert sp.cancel(oracle(f)-sum(oracle(a)*oracle(b) for a, b in zip(q, g))-oracle(r)) == 0


@pytest.mark.parametrize("order", ORDERS)
@pytest.mark.parametrize("options", [
    {"domain": "QQ(a,b)"}, {"domain": "QQ(a,b)", "auto": False},
    {"domain": "ZZ(a,b)", "polys": True},
    {"domain": "QQ[a,b]", "auto": False},
    {"domain": "ZZ[a,b]", "auto": False},
    {"auto": False}, {}, {"polys": True},
])
def test_reduced_fraction_domains_match_sympy(order, options):
    x, y, a, b = sm.symbols("x y a b")
    sx, sy = sp.symbols("x y")
    f = (a+1)*x*x + b*y + 1
    divisors = [a*x-y, y*y-b]
    q, r = sm.reduced(f, divisors, x, y, order=order, **options)
    oq, oremainder = sp.reduced(oracle(f), [oracle(e) for e in divisors], sx, sy, order=order, **options)
    assert all(equal(a, b.as_expr() if isinstance(b, sp.Poly) else b) for a, b in zip(q, oq))
    assert equal(r, oremainder.as_expr() if isinstance(oremainder, sp.Poly) else oremainder)
    if options.get("polys"):
        assert all(isinstance(p, sm.Poly) for p in q+[r])
        assert [str(p.domain) for p in q+[r]] == [str(p.domain) for p in oq+[oremainder]]


def test_native_poly_fraction_arithmetic_and_metadata():
    x, a, b = sm.symbols("x a b")
    p = sm.Poly(x/a+1/(a+1), x, domain="QQ(a,b)")
    q = sm.Poly((a+1)*x-b, x, domain="QQ(a,b)")
    sx, sa, sb = sp.symbols("x a b")
    op = sp.Poly(sx/sa+1/(sa+1), sx, domain="QQ(a,b)")
    oq = sp.Poly((sa+1)*sx-sb, sx, domain="QQ(a,b)")
    assert isinstance(p, sm.Poly) and str(p.domain) == "QQ(a,b)"
    assert p.gens == (x,) and p.degree() == 1 and p.total_degree() == 1
    assert p.monoms() == [(1,), (0,)]
    assert equal(p.LC(), op.LC()) and equal(p.TC(), op.TC())
    assert all(equal(a, b) for a, b in zip(p.all_coeffs(), op.all_coeffs()))
    for result, expected in [(p+q, op+oq), (p-q, op-oq), (p*q, op*oq), (-p, -op),
                             (p**3, op**3), (p.monic(), op.monic()),
                             (p.diff(), op.diff()), (p.integrate(), op.integrate()),
                             (p.compose(q), op.compose(oq)), (p.shift(2), op.shift(2))]:
        assert isinstance(result, sm.Poly)
        assert str(result.domain) == "QQ(a,b)"
        assert equal(result, expected.as_expr())
    quotient, remainder = divmod(p*q, p)
    assert quotient == q and remainder.is_zero
    assert (p*q).exquo(p) == q
    for result, expected in [(p.gcd(p*q), op.gcd(op*oq)), (p.lcm(p*q), op.lcm(op*oq))]:
        assert equal(result, expected.as_expr())
    content, primitive = p.primitive()
    oc, opp = op.primitive()
    assert equal(content, oc) and equal(primitive, opp.as_expr())
    denominator, cleared = p.clear_denoms()
    od, ocp = op.clear_denoms()
    assert equal(denominator, od) and equal(cleared, ocp.as_expr())
    assert equal(p.content(), op.content())


def test_fraction_poly_factorization_certificates():
    x, a = sm.symbols("x a")
    p = sm.Poly((x/a+1/(a+1))**2*(x-a), x, domain="QQ(a)")
    for coefficient, factors in (p.factor_list(), p.sqf_list()):
        product = oracle(coefficient)
        for f, k in factors:
            assert isinstance(f, sm.Poly)
            product *= oracle(f)**k
        assert sp.cancel(product-oracle(p)) == 0
    assert not p.is_sqf
    assert not p.is_irreducible
    assert equal(p.sqf_part(), sp.Poly(oracle(p), sp.Symbol("x"), domain="QQ(a)").sqf_part().as_expr())


def test_fraction_domains_inference_assumptions_and_zero_unit_ideals():
    x = sm.Symbol("x", real=True)
    a = sm.Symbol("a", nonzero=True)
    g = sm.groebner([a*x-1], x, domain="QQ(a)", polys=True)
    assert g.gens == (x,) and isinstance(g[0], sm.Poly)
    assert str(g[0].domain) == "QQ(a)" and g[0].gens == (x,)
    assert g.contains(a*x-1)
    assert g.fglm("lex").exprs == g.exprs
    converted = g.fglm("grlex")
    assert converted.order == "grlex" and converted.domain == g.domain
    assert converted.exprs == g.exprs
    assert sm.groebner([a*x-1], domain="QQ(a)").gens == (x,)
    assert sm.Poly(x/a, x, domain="QQ(a)").gens == (x,)
    empty = sm.groebner([], x, domain="QQ(a)")
    assert list(empty) == [] and empty.reduce(x) == ([], x)
    unit = sm.groebner([a], x, domain="QQ(a)")
    assert list(unit) == [1] and unit.contains(x)
    assert not unit.is_zero_dimensional


def test_fraction_domains_reject_unknown_coefficients_and_zero_division():
    x, a, b = sm.symbols("x a b")
    for constructor in (lambda e: sm.Poly(e, x, domain="QQ(a)"),
                        lambda e: sm.groebner([e], x, domain="QQ(a)"),
                        lambda e: sm.reduced(e, [x], x, domain="QQ(a)")):
        with pytest.raises(ValueError, match="declared parameters"):
            constructor(b*x)
        with pytest.raises(ValueError, match="exact"):
            constructor(1.5*x)
        with pytest.raises(ValueError):
            constructor(sm.sqrt(2)*x)
    with pytest.raises(ValueError, match="ring"):
        sm.Poly(x/a, x, domain="QQ[a]")
    p = sm.Poly(x/a, x, domain="QQ(a)")
    with pytest.raises(ZeroDivisionError):
        p.div(0)
    with pytest.raises(TypeError, match="norm"):
        p.max_norm()
    with pytest.raises(TypeError):
        p.gcd(sm.Poly(sm.pi*x, x))
    with pytest.raises(ValueError, match="generators"):
        sm.Poly(x/a, x, x, domain="QQ(a)")
    with pytest.raises(ValueError, match="generators"):
        sm.Poly(x/a, 1, domain="QQ(a)")


def test_fraction_domains_have_no_sympy_runtime_dependency():
    subprocess.run([sys.executable, "-c", """
import sys
sys.modules['sympy'] = None
import symox as s
x,y,a=s.symbols('x y a')
g=s.groebner([a*x-1,y-x],x,y,domain='QQ(a)')
assert g.contains(a*x-1)
p=s.Poly(x/a+1/(a+1),x,domain='QQ(a)')
q,r=(p*p).div(p)
assert q==p and r.is_zero
"""], check=True, timeout=30)


@pytest.mark.parametrize("protocol", range(pickle.HIGHEST_PROTOCOL+1))
def test_pickle_preserves_coefficient_domains_and_unused_parameter_assumptions(protocol):
    x = sm.Symbol("x", real=True)
    a = sm.Symbol("a", positive=True)
    b = sm.Symbol("b", negative=True)
    # Bind b's assumptions through the input, then eliminate its occurrence.
    # A string domain alone cannot discover an unused Symbol's assumptions.
    p = sm.Poly(x/a+1+b, x, domain="QQ(a,b)")-b
    q = pickle.loads(pickle.dumps(p, protocol=protocol))
    assert q == p and q.domain == p.domain and q.gens == p.gens
    assert isinstance(q+b, sm.Poly)  # b is unused, but retains its assumptions
    g = sm.groebner([a*x-1, b*(a*x-1)], x, domain="QQ(a,b)", polys=True)
    h = pickle.loads(pickle.dumps(g, protocol=protocol))
    assert h.domain == g.domain and h.gens == g.gens and h.order == g.order
    assert isinstance(h[0], sm.Poly)
    quotients, remainder = h.reduce(b*x)
    assert quotients[0].as_expr() == b
    assert remainder.as_expr() == b/a
    for obj in (sm.Poly(x*x-1, x, domain="QQ"), sm.Poly(x*x+1, x, modulus=5),
                sm.groebner([x*x-1], x, modulus=5),
                sm.groebner([], x, domain="QQ(a,b)")):
        rebuilt = pickle.loads(pickle.dumps(obj, protocol=protocol))
        assert rebuilt.domain == obj.domain and rebuilt.gens == obj.gens


def test_poly_field_evaluation_and_multiple_fields_preserve_native_storage():
    x, y, a, b = sm.symbols("x y a b")
    p = sm.Poly(x/a+y/(a+1), x, y, domain="QQ(a,b)")
    evaluated = p.eval(x, 2)
    assert isinstance(evaluated, sm.Poly)
    assert evaluated.domain == p.domain and evaluated.gens == (y,)
    assert equal(evaluated, oracle(p).subs(sp.Symbol("x"), 2))
    first = sm.Poly(x/a+1, x, domain="QQ(a)")
    second = sm.Poly(x/b+1, x, domain="QQ(b)")
    result = first+second
    assert isinstance(result, sm.Poly) and str(result.domain) == "QQ(a,b)"
    assert equal(result, oracle(first)+oracle(second))
    assert first.termwise(lambda monomial, c: 2*c).domain == first.domain
    zz = sm.Poly(x/a+1, x, domain="ZZ(a)")
    assert str((zz+zz).domain) == "ZZ(a)"


@pytest.mark.parametrize("domain", ["QQ(a]", "QQ[a)", "QQ()", "QQ(a,a)", "QQ(a+x)"])
def test_malformed_domains_raise_python_errors(domain):
    x, a = sm.symbols("x a")
    for build in (lambda: sm.Poly(x/a, x, domain=domain),
                  lambda: sm.groebner([x/a], x, domain=domain),
                  lambda: sm.reduced(x/a, [x], x, domain=domain)):
        with pytest.raises((ValueError, NotImplementedError)):
            build()
