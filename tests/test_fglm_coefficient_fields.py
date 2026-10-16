"""Order conversion over exact symbolic coefficient fields and rings."""
import itertools
import pickle

import pytest
import sympy as sp
import symox as sm


ORDERS = ("lex", "grlex", "grevlex")


@pytest.mark.parametrize("source,target", list(itertools.permutations(ORDERS, 2)))
@pytest.mark.parametrize("domain", ["QQ(a,b)", "ZZ(a,b)", "QQ[a,b]", "ZZ[a,b]"])
@pytest.mark.parametrize("sources", [
    ["a*x*y-b", "y**2-x"],
    ["a*x**2+y-b", "x*y-1"],
])
def test_fglm_matches_direct_basis_and_sympy(source, target, domain, sources):
    a, b, x, y = sm.symbols("a b x y")
    values = [sm.sympify(e) for e in sources]
    original = sm.groebner(values, x, y, domain=domain, order=source)
    converted = original.fglm(target)
    expected = sp.groebner([e.to_sympy() for e in values], sp.Symbol("x"), sp.Symbol("y"), domain=domain, order=source).fglm(target)
    direct = sm.groebner(values, x, y, domain=domain, order=target)
    assert converted.order == target and converted.domain == original.domain
    assert converted.gens == original.gens
    assert len(converted) == len(expected)
    for actual, reference, recomputed in zip(converted, expected, direct):
        assert sp.cancel(actual.to_sympy()-reference) == 0
        assert sm.cancel(actual-recomputed) == 0
    for value in values:
        assert converted.contains(value)
    assert converted.fglm(source).exprs == original.exprs
    assert pickle.loads(pickle.dumps(converted)).exprs == converted.exprs


def test_fglm_preserves_assumptions_and_polynomial_output():
    a = sm.Symbol("a", positive=True)
    x, y = sm.symbols("x y", real=True)
    g = sm.groebner([a*x*y-1, y*y-x], x, y, domain="QQ(a)", order="grevlex", polys=True)
    converted = g.fglm("lex")
    assert all(isinstance(p, sm.Poly) for p in converted)
    assert converted.gens == (x, y)
    assert set().union(*(p.as_expr().free_symbols for p in converted)) <= {a, x, y}
    quotients, remainder = converted.reduce(a*x*y-1)
    assert all(isinstance(p, sm.Poly) for p in quotients)
    assert remainder.is_zero


def test_fglm_rejects_positive_dimension_and_preserves_same_order():
    a, x, y = sm.symbols("a x y")
    for domain in ("QQ(a)", "ZZ[a]"):
        g = sm.groebner([a*x*y], x, y, domain=domain, order="lex")
        assert g.fglm("lex").exprs == g.exprs
        with pytest.raises(NotImplementedError, match="positive dimension"):
            g.fglm("grevlex")
