"""Expression grounds retain sparse support, exact stored values and fallbacks."""
import pickle
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx

N = 1 << 31
DOMAINS = ['EX', 'ComplexField(53)', 'ComplexField(140)',
           'RealField(140).poly_ring(a)', 'ComplexField(140).poly_ring(a)',
           'RealField(140).frac_field(a)', 'ComplexField(140).frac_field(a)']


def setup(name):
    x, y, z, a = sx.symbols('x y z a')
    K = eval(name, dict(vars(sx), a=a))
    value = sx.Float('1.23456789012345678901234567890123456789', precision=140)
    c = value*sx.I if name.startswith('Complex') else value
    if 'ring' in name:
        c *= a
    elif 'field(a)' in name:
        c = c*a/(a+1)
    if name == 'EX':
        c += sx.sin(a)
    return (x, y, z, a), K, c


def exact(value):
    expression = sx.to_sympy(value)
    return expression.xreplace({f: sp.Rational(f) for f in expression.atoms(sp.Float)})


def exact_dict(p):
    return {m: exact(c) for m, c in p.as_dict().items()}


@pytest.mark.parametrize('name', DOMAINS)
def test_large_expression_ground_construction_views_and_structural_maps(name):
    (x, y, z, a), K, c = setup(name)
    p = sx.Poly.from_dict({(N, N): c, (N, 0): 2}, x, y, domain=K)
    oracle_K = eval(name, dict(vars(sp), a=sp.Symbol('a')))
    oracle_c = oracle_K.to_sympy(oracle_K.from_sympy(sx.to_sympy(c)))
    oracle_c = oracle_c.xreplace({f: sp.Rational(f) for f in oracle_c.atoms(sp.Float)})
    expected = {(N, N): oracle_c, (N, 0): sp.Integer(2)}
    assert exact_dict(p) == expected
    assert p.degree_list() == (N, N) and p.total_degree() == 2*N
    assert p.gens == (x, y) and p.domain == K
    assert p.monoms() == [(N, N), (N, 0)]
    assert p.TC() == 0
    assert exact(p.as_expr()) == expected[(N, N)]*sp.Symbol('x')**N*sp.Symbol('y')**N + 2*sp.Symbol('x')**N
    reordered = p.reorder(y, x)
    assert exact_dict(reordered) == {(N, N): expected[(N, N)], (0, N): 2}
    assert reordered.reorder(x, y) == p
    h = p.homogenize(z)
    assert exact_dict(h) == {(N, N, 0): expected[(N, N)], (N, 0, N): 2}
    assert h.domain == K and h.total_degree() == 2*N
    assert pickle.loads(pickle.dumps(p)) == p
    rebuilt = sx.Poly(p.as_expr(), x, y, domain=K)
    assert rebuilt == p and exact_dict(rebuilt) == expected


@pytest.mark.parametrize('name', DOMAINS)
def test_sparse_additive_scalars_and_zero_protocols_have_exact_certificates(name):
    (x, y, z, a), K, c = setup(name)
    p = sx.Poly.from_dict({(N, N): c, (N, 0): 2}, x, y, domain=K)
    source = exact_dict(p)
    for result in (p+p, 2*p, p*2):
        assert result.domain == K
        assert exact_dict(result) == {m: 2*v for m, v in source.items()}
    assert exact_dict(-p) == {m: -v for m, v in source.items()}
    assert (p-p).is_zero and (p-p).domain == K
    assert (p*0).is_zero and (0*p).is_zero
    assert p**0 == sx.Poly(1, x, y, domain=K)
    assert p**1 == p
    if '140' in name:
        assert {f._prec for value in (p+p).coeffs() for f in value.atoms(sx.Float)} == {140}
    q = sx.Poly.from_dict({(N, 0): 2}, x, y, domain=K)
    assert exact_dict(p-q) == {(N, N): source[(N, N)]}


@pytest.mark.parametrize('name', DOMAINS)
def test_domain_conversion_preserves_large_support_and_selected_precision(name):
    (x, y, z, a), K, c = setup(name)
    plain = sx.Poly.from_dict({((1 << 32)-1,): 2}, x)
    p = plain.set_domain(K)
    assert p.monoms() == [((1 << 32)-1,)] and p.degree() == (1 << 32)-1
    assert p.domain == K and exact_dict(p) == {((1 << 32)-1,): 2}
    assert pickle.loads(pickle.dumps(p)) == p


@pytest.mark.parametrize('name', ['EX', 'ComplexField(53)', 'ComplexField(140)'])
@pytest.mark.parametrize('operation', ['mul', 'pow', 'div', 'diff', 'integrate', 'eval'])
def test_compact_inputs_keep_the_complete_dense_fallback(name, operation):
    (x, y, z, a), K, c = setup(name)
    p = sx.Poly(x**5 + 2*x + 1, x, domain=K)
    q = sx.Poly(x+1, x, domain=K)
    X = sp.Symbol('x')
    oracle_K = eval(name, vars(sp))
    F, G = sp.Poly(X**5+2*X+1, X, domain=oracle_K), sp.Poly(X+1, X, domain=oracle_K)
    actions = {'mul': lambda p, q: p*q, 'pow': lambda p, q: p**3,
               'div': lambda p, q: p.div(q), 'diff': lambda p, q: p.diff(),
               'integrate': lambda p, q: p.integrate(), 'eval': lambda p, q: p.eval(2)}
    got, wanted = actions[operation](p, q), actions[operation](F, G)
    def representation(value):
        if isinstance(value, tuple):
            return tuple(representation(v) for v in value)
        if hasattr(value, 'as_dict'):
            return {m: exact(c) if isinstance(c, sx.Expr) else c.xreplace({f: sp.Rational(f) for f in c.atoms(sp.Float)}) for m, c in value.as_dict().items()}
        if isinstance(value, sx.Expr):
            return exact(value)
        return value.xreplace({f: sp.Rational(f) for f in value.atoms(sp.Float)}) if isinstance(value, sp.Basic) else value
    assert representation(got) == representation(wanted)


def test_large_sparse_ex_groebner_roundtrip_and_explicit_overflow():
    source = '''import symox as s
x,y=s.symbols('x y');N=2**31
for method in ('buchberger','f5b'):
    g=s.groebner([x-y**N],x,y,domain=s.EX,method=method)
    assert g.contains(x-y**N)
    q,r=g.reduce(x-y**N)
    assert q==[1] and r==0
    try: g.contains(x*y**N)
    except NotImplementedError as error: assert 'exponent storage limit' in str(error)
    else: raise AssertionError('unrepresentable nonzero remainder returned membership')
'''
    result = subprocess.run([sys.executable, '-B', '-c', source], capture_output=True, text=True, timeout=20)
    assert result.returncode == 0, result.stderr
    assert result.stderr == ''
