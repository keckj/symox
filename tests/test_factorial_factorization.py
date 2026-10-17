"""Greatest factorial factors: exact shift chains, never ordinary powers."""
import pickle
import subprocess
import sys

import pytest
import symox as sm

from test_global_poly_division import selected


def shifted(p, amount):
    x, = p.gens
    # Exact polynomial composition is an independent public path from the
    # checked sparse shifts used by the factorial-factor algorithm.
    return p.compose(sm.Poly(x+amount, x, domain=p.domain))


def certify(p, factors):
    assert isinstance(factors, list)
    one = sm.Poly(1, *p.gens, domain=p.domain)
    product = one
    chains = []
    previous = 0
    for h, length in factors:
        assert isinstance(h, sm.Poly)
        assert h.gens == p.gens and h.domain == p.domain
        assert isinstance(length, int) and length > previous
        assert h.degree() > 0 and h.LC() == 1
        assert all(c in p.domain for c in h.coeffs())
        previous = length
        restored = pickle.loads(pickle.dumps(h))
        assert restored.domain == p.domain and restored.gens == p.gens
        assert restored.as_expr() == h.as_expr()
        chain = one
        for j in range(length):
            chain *= shifted(h, -j)
        product *= chain
        chains.append(chain)
    # Compare to f/LC(f) without asking monic() to promote an input ring.
    assert sm.Poly(sm.expand(product.as_expr()*p.LC()-p.as_expr()), *p.gens, domain=p.domain).is_zero
    assert sum(k*h.degree() for h, k in factors) == p.degree()
    for i, chain in enumerate(chains):
        for h, length in factors[i:]:
            assert chain.gcd(shifted(h, 1)).degree() == 0
            assert chain.gcd(shifted(h, -length)).degree() == 0


def oracle_factors(actual, expected):
    import sympy as sy
    assert [k for _, k in actual] == [k for _, k in expected]
    for (h, _), (oh, _) in zip(actual, expected):
        assert tuple(sm.to_sympy(g) for g in h.gens) == oh.gens
        assert sy.expand(sm.to_sympy(h.as_expr())-oh.as_expr()) == 0
        assert str(h.domain) == str(oh.domain)


@pytest.mark.parametrize('kind', ['ZZ', 'QQ', 'AN', 'QQRF', 'ZZRF', 'QQring', 'ZZring'])
@pytest.mark.parametrize('order', ['lex', 'grlex', 'grevlex'])
@pytest.mark.parametrize('case', range(4))
def test_typed_factorial_chains_oracle_certificates_and_metadata(kind, order, case):
    import sympy as sy
    x = sm.Symbol('x', real=True)
    k, ok, a, unused = selected(kind, order)
    c = a/(unused+1) if kind.endswith('RF') else a if kind.endswith('ring') else sm.sqrt(2) if kind == 'AN' else sm.Integer(2)
    y = x+c
    expression = [c*y*(y-1)*(y-2), y**3*(y-1)**2*(y-3),
                  (y*y+1)*((y-1)**2+1)*(y+4), y**2+1][case]
    p = sm.Poly(expression, x, domain=k)
    op = sy.Poly(sm.to_sympy(expression), sm.to_sympy(x), domain=ok)
    factors = p.gff_list()
    oracle_factors(factors, op.gff_list())
    certify(p, factors)
    oracle_factors(sm.gff_list(p), op.gff_list())
    oracle_factors(sm.gff_list(expression, x, domain=k, polys=True), op.gff_list())
    if kind.endswith('RF') or kind.endswith('ring'):
        for h, _ in factors:
            assert h.domain.symbols == (a, unused)
            assert h.domain.order == order


@pytest.mark.parametrize('case', range(5))
def test_documented_and_overlapping_chain_forms(case):
    import sympy as sy
    x = sm.Symbol('x')
    expressions = [x**5+2*x**4-x**3-2*x*x,
                   (x**3+7)*((x-1)**3+7)*(x*x+5*x)*((x-1)**2+5*(x-1))*((x-2)**2+5*(x-2)),
                   x*x*(x-sm.Rational(1,3))*(x-sm.Rational(1,2))*(x-1)*(x-sm.Rational(7,3))**3*(x-3)**2,
                   x**8*(x-1)**4*(x-2)**2, (x*x+1)**3*((x-1)**2+1)**2]
    p = sm.Poly(expressions[case], x, domain=sm.QQ)
    op = sy.Poly(sm.to_sympy(p.as_expr()), sm.to_sympy(x), domain=sy.QQ)
    fs = p.gff_list()
    oracle_factors(fs, op.gff_list())
    certify(p, fs)


@pytest.mark.parametrize('kind', ['ZZ', 'QQ', 'AN', 'QQRF', 'ZZRF', 'QQring', 'ZZring'])
@pytest.mark.parametrize('constant', [1, -7, 0])
def test_constants_and_zero(kind, constant):
    domain, _, _, _ = selected(kind)
    x = sm.Symbol('x')
    p = sm.Poly(constant, x, domain=domain)
    if constant:
        assert p.gff_list() == sm.gff_list(p) == []
    else:
        with pytest.raises(ValueError, match='zero polynomial'):
            p.gff_list()
        with pytest.raises(ValueError, match='zero polynomial'):
            sm.gff_list(p)


@pytest.mark.parametrize('unit', [-2, 2, sm.Rational(1,3), sm.Rational(-7,5)])
def test_scalar_content_is_excluded_from_monic_factorial_form(unit):
    x = sm.Symbol('x')
    p = sm.Poly(unit*x*(x-1), x, domain=sm.QQ)
    fs = p.gff_list()
    assert fs == [(sm.Poly(x, x, domain=sm.QQ), 2)]
    certify(p, fs)


@pytest.mark.parametrize('kind', ['ZZ', 'QQring', 'ZZring'])
def test_monic_normalization_must_be_exact_inside_original_ring(kind):
    import sympy as sy
    x = sm.Symbol('x')
    k, ok, a, _ = selected(kind)
    c = sm.Integer(2) if kind == 'ZZ' else a
    p = sm.Poly((c*x+1)*(c*x+1-c), x, domain=k)
    op = sy.Poly(sm.to_sympy(p.as_expr()), sm.to_sympy(x), domain=ok)
    with pytest.raises(sy.polys.polyerrors.ExactQuotientFailed):
        op.gff_list()
    with pytest.raises(sm.ExactQuotientFailed):
        p.gff_list()
    with pytest.raises(sm.ExactQuotientFailed):
        sm.gff_list(p)
    promoted = sm.gff_list(p, field=True)
    assert all(h.domain == k.get_field() for h, _ in promoted)
    certify(sm.Poly(p.as_expr(), x, domain=k.get_field()), promoted)


@pytest.mark.parametrize('generator', [sm.sqrt(2), sm.I, 2**sm.Rational(1,3)+sm.I,
                                      sm.sqrt(2)+sm.sqrt(3)+sm.sqrt(5)])
def test_number_field_embeddings_alias_and_scopes(generator):
    import sympy as sy
    x = sm.Dummy('x', positive=True)
    alias = sm.Symbol('alpha')
    k = sm.QQ.algebraic_field(generator, alias=alias)
    p = sm.Poly((x+generator)*(x+generator-1)*(x+generator-2), x, domain=k)
    fs = p.gff_list()
    assert [length for _, length in fs] == [3]
    assert fs[0][0].domain == k and fs[0][0].domain.ext.alias == alias
    certify(p, fs)
    op = sy.Poly(sm.to_sympy(p.as_expr()), sm.to_sympy(x), extension=sm.to_sympy(generator))
    assert sy.expand(sm.to_sympy(fs[0][0].as_expr())-op.gff_list()[0][0].as_expr()) == 0
    # A different embedded field and a plain rational operation must remain
    # independent of the preceding factory scope.
    assert sm.Poly(x*(x-1), x, domain=sm.QQ).gff_list()[0][1] == 2


@pytest.mark.parametrize('kind', ['expression', 'poly'])
@pytest.mark.parametrize('polys', [None, False, True])
def test_global_return_rendering_options(kind, polys):
    import sympy as sy
    x = sm.Symbol('x')
    expression = x*(x-1)*(x-2)
    f = sm.Poly(expression, x) if kind == 'poly' else expression
    of = sy.Poly(sm.to_sympy(expression), sm.to_sympy(x)) if kind == 'poly' else sm.to_sympy(expression)
    actual = sm.gff_list(f, x, polys=polys)
    expected = sy.gff_list(of, sm.to_sympy(x), polys=polys)
    assert len(actual) == len(expected) == 1 and actual[0][1] == expected[0][1] == 3
    assert isinstance(actual[0][0], sm.Poly) == isinstance(expected[0][0], sy.Poly)
    ah = actual[0][0].as_expr() if isinstance(actual[0][0], sm.Poly) else actual[0][0]
    eh = expected[0][0].as_expr() if isinstance(expected[0][0], sy.Poly) else expected[0][0]
    assert sm.to_sympy(ah) == eh


@pytest.mark.parametrize('kwargs', [{'auto':True}, {'auto':'bad'}, {'bogus':True}, {'polys':2}])
def test_options_are_applied_or_explicitly_rejected(kwargs):
    x = sm.Symbol('x')
    expected = sm.FlagError if 'auto' in kwargs else sm.OptionError
    with pytest.raises(expected):
        sm.gff_list(x, x, **kwargs)


@pytest.mark.parametrize('prime', [2, 3, 5, 1000000007])
def test_positive_characteristic_is_rejected_before_nondescending_gcd(prime):
    # SymPy 1.14 recurses forever on x**p-x over GF(p), since shifting by one
    # fixes that polynomial. Characteristic-zero GFF must reject this input.
    x = sm.Symbol('x')
    p = sm.Poly(x**prime-x, x, modulus=prime)
    with pytest.raises(ValueError, match='characteristic zero'):
        p.gff_list()


def test_multivariate_opaque_inexact_and_huge_sparse_inputs_are_explicit():
    x, y, a = sm.symbols('x y a')
    for p in [sm.Poly(x*y+1, x, y), sm.Poly(x+1, x, domain=sm.GF(5))]:
        with pytest.raises(ValueError):
            p.gff_list()
    for p in [sm.Poly(0.5*x+1, x), sm.Poly(sm.sin(a)*x+1, x)]:
        with pytest.raises(ValueError, match='exact coefficient domain'):
            p.gff_list()
    p = sm.Poly(x**1000000000+1, x)
    with pytest.raises(ValueError, match='resource limit'):
        p.gff_list()
    assert sm.Poly(x*(x-1), x).gff_list()[0][1] == 2


def test_native_factorial_factorization_without_sympy_runtime():
    program = r'''
import sys, importlib.abc, pickle
class NoSymPy(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == 'sympy' or fullname.startswith('sympy.'):
            raise AssertionError('runtime SymPy import')
sys.meta_path.insert(0, NoSymPy())
import symox as s
from symox.polys.polytools import gff_list
x=s.Dummy('x',real=True)
a=s.Symbol('a',positive=True)
b=s.Dummy('unused')
for k,c in [(s.ZZ,s.Integer(2)),(s.QQ,s.Rational(1,3)),(s.QQ.frac_field(a,b,order='grevlex'),a/(b+1)),(s.ZZ.poly_ring(a,b,order='grlex'),a),(s.QQ.algebraic_field(s.I),s.I)]:
    p=s.Poly(c*(x+c)*(x+c-1)*(x+c-2),x,domain=k)
    fs=gff_list(p)
    assert len(fs)==1 and fs[0][1]==3
    h=pickle.loads(pickle.dumps(fs[0][0]))
    assert h.domain==k and h.gens==p.gens
    chain=s.Poly(1,x,domain=k)
    for j in range(3):chain*=h.compose(s.Poly(x-j,x,domain=k))
    assert s.Poly(s.expand(chain.as_expr()*p.LC()-p.as_expr()),x,domain=k).is_zero
    assert all(t in k for t in h.coeffs())
assert not any(n=='sympy' or n.startswith('sympy.') for n in sys.modules)
'''
    subprocess.run([sys.executable, '-c', program], check=True, timeout=30)
