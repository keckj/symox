"""Typed Brown--Traub PRS: exact normalization, elimination and certificates."""
import pickle
import subprocess
import sys

import pytest
import symox as sm

from test_global_poly_division import selected


def compare(actual, expected):
    import sympy as sy
    if isinstance(expected, (list, tuple)):
        assert isinstance(actual, type(expected)) and len(actual) == len(expected)
        for a, e in zip(actual, expected):
            compare(a, e)
    elif isinstance(expected, sy.Poly):
        assert isinstance(actual, sm.Poly)
        assert tuple(sm.to_sympy(g) for g in actual.gens) == expected.gens
        assert sy.expand(sm.to_sympy(actual.as_expr()) - expected.as_expr()) == 0
        assert str(actual.domain) == str(expected.domain)
    else:
        assert not isinstance(actual, sm.Poly)
        assert sy.expand(sm.to_sympy(actual) - expected) == 0


def expression(value):
    return value.as_expr() if isinstance(value, sm.Poly) else value


def in_domain_zero(value, generators, domain):
    assert sm.Poly(sm.expand(value), *generators, domain=domain).is_zero


def leading(poly):
    degree = int(poly.degree())
    return sm.Add(*(coefficient * sm.Mul(*(generator**power for generator, power in zip(poly.gens[1:], powers[1:])))
                    for powers, coefficient in poly.as_dict().items() if powers[0] == degree))


def determinant(matrix):
    if not matrix:
        return sm.Integer(1)
    return sm.Add(*((-1)**j * matrix[0][j] * determinant([row[:j]+row[j+1:] for row in matrix[1:]])
                    for j in range(len(matrix))))


def sylvester(poly, other):
    n, m = int(poly.degree()), int(other.degree())
    coefficient = lambda p, k: sm.Add(*(c * sm.Mul(*(v**e for v, e in zip(p.gens[1:], powers[1:])))
                                       for powers, c in p.as_dict().items() if powers[0] == k))
    matrix = [[sm.Integer(0) for _ in range(n+m)] for _ in range(n+m)]
    for row in range(m):
        for k in range(n+1):
            matrix[row][row+k] = coefficient(poly, n-k)
    for row in range(n):
        for k in range(m+1):
            matrix[m+row][row+k] = coefficient(other, m-k)
    return determinant(matrix)


@pytest.mark.parametrize('kind', ['ZZ', 'QQ', 'GF', 'AN', 'QQRF', 'ZZRF', 'QQring', 'ZZring'])
@pytest.mark.parametrize('order', ['lex', 'grlex', 'grevlex'])
@pytest.mark.parametrize('case', range(5))
def test_typed_prs_oracle_domains_and_elimination(kind, order, case):
    import sympy as sy
    x, y = sm.symbols('x y')
    domain, oracle_domain, a, b = selected(kind, order)
    c = a/(b+1) if kind.endswith('RF') else a if kind.endswith('ring') else sm.sqrt(2) if kind == 'AN' else sm.Integer(2)
    f, g = [
        (c*x*x+(y+1)*x+1, (c+1)*x+y),
        (x**3+c*x*x+y*x+1, x*x+c*x+y+1),
        (x**4+1, x*x+1),
        (c*(x*x+x+1), 2*c*(x+1)),
        (x*x+y*x+1, x*x+y*x+1),
    ][case]
    p, q = sm.Poly(f, x, y, domain=domain), sm.Poly(g, x, y, domain=domain)
    op = sy.Poly(sm.to_sympy(f), sm.to_sympy(x), sm.to_sympy(y), domain=oracle_domain)
    oq = sy.Poly(sm.to_sympy(g), sm.to_sympy(x), sm.to_sympy(y), domain=oracle_domain)
    sequence = p.subresultants(q)
    compare(sequence, op.subresultants(oq))
    compare(p.resultant(q), op.resultant(oq))
    compare(p.resultant(q, includePRS=True), op.resultant(oq, includePRS=True))
    compare(p.discriminant(), op.discriminant())
    compare(sm.subresultants(p, q), op.subresultants(oq))
    compare(sm.resultant(p, q, includePRS=True), op.resultant(oq, includePRS=True))
    compare(sm.discriminant(p), op.discriminant())
    compare(sm.subresultants(f, g, x, y, domain=domain, polys=True), sy.subresultants(sm.to_sympy(f), sm.to_sympy(g), sm.to_sympy(x), sm.to_sympy(y), domain=oracle_domain, polys=True))
    eliminated = [p.resultant(q), p.discriminant()]
    for value, generators in [(v, (x,y)) for v in sequence] + [(v, (y,)) for v in eliminated]:
        assert isinstance(value, sm.Poly) and value.domain == domain
        assert value.gens == generators
        assert all(coefficient in domain for coefficient in value.coeffs())
        restored = pickle.loads(pickle.dumps(value))
        assert restored.gens == value.gens and restored.domain == domain and restored.as_expr() == value.as_expr()
        if kind.endswith('RF') or kind.endswith('ring'):
            assert value.domain.symbols == (a, b) and value.domain.order == order
    for previous, current, following in zip(sequence, sequence[1:], sequence[2:]):
        assert int(following.degree()) < int(current.degree())
        # Each Brown normalization is exact in the original coefficient ring.
        assert previous.prem(current).rem(following, auto=False).is_zero


@pytest.mark.parametrize('kind', ['ZZ', 'QQ', 'GF', 'AN', 'QQRF', 'ZZring'])
def test_independent_sylvester_antisymmetry_and_discriminant_certificates(kind):
    x, y = sm.symbols('x y')
    domain, _, a, _ = selected(kind)
    c = a if kind in ('QQRF', 'ZZring') else sm.sqrt(2) if kind == 'AN' else sm.Integer(2)
    p = sm.Poly(x**3+c*x*x+y*x+1, x, y, domain=domain)
    q = sm.Poly(x*x+c*x+y+1, x, y, domain=domain)
    result = expression(p.resultant(q))
    in_domain_zero(result-sylvester(p, q), (y,), domain)
    in_domain_zero(expression(q.resultant(p))-result, (y,), domain)
    disc = expression(p.discriminant())
    in_domain_zero(disc*leading(p)+expression(p.resultant(p.diff())), (y,), domain)


@pytest.mark.parametrize('prime', [2, 3, 5, 7])
@pytest.mark.parametrize('f,g', [(0,0), (0,'x'), ('x',0), (2,3), ('x',2), ('x**2','2*x'), ('x','x**2'), ('x**5+1','x**2+1')])
def test_degenerate_cases_and_characteristic(prime, f, g):
    import sympy as sy
    x = sm.Symbol('x')
    convert = lambda value: sm.sympify(value)
    p, q = sm.Poly(convert(f), x, modulus=prime), sm.Poly(convert(g), x, modulus=prime)
    op = sy.Poly(sm.to_sympy(p.as_expr()), sm.to_sympy(x), modulus=prime)
    oq = sy.Poly(sm.to_sympy(q.as_expr()), sm.to_sympy(x), modulus=prime)
    compare(p.subresultants(q), op.subresultants(oq))
    expected = op.resultant(oq, includePRS=True)
    if not op.is_zero and not oq.is_zero and op.degree() < oq.degree() and op.degree()%2 and oq.degree()%2:
        expected = (-expected[0], expected[1])
    compare(p.resultant(q, includePRS=True), expected)
    compare(p.discriminant(), op.discriminant())


def test_resultant_restores_original_order_sign_even_when_sympy_does_not():
    import sympy as sy
    x = sm.Symbol('x')
    p, q = sm.Poly(x+1, x), sm.Poly(x**3+2, x)
    assert p.resultant(q) == 1 == sylvester(p, q)
    assert q.resultant(p) == -1
    assert p.subresultants(q) == q.subresultants(p)
    # SymPy 1.14 reorders the PRS and loses (-1)^(deg(f)*deg(g)).
    # Certify the mathematical result independently rather than reproducing it.
    assert sm.to_sympy(p.resultant(q)) == sy.expand(sm.to_sympy(q.as_expr()).subs(sm.to_sympy(x), -1))


@pytest.mark.parametrize('name', ['resultant', 'discriminant', 'subresultants'])
@pytest.mark.parametrize('polys', [None, False, True])
@pytest.mark.parametrize('kind', ['expr', 'poly', 'mixed'])
def test_global_rendering_and_flags(name, polys, kind):
    import sympy as sy
    x, y = sm.symbols('x y')
    f, g = x*x+x*y+1, 2*x+y
    of, og = sm.to_sympy(f), sm.to_sympy(g)
    if kind != 'expr':
        f, of = sm.Poly(f,x,y), sy.Poly(of,sm.to_sympy(x),sm.to_sympy(y))
    if kind == 'poly':
        g, og = sm.Poly(g,x,y), sy.Poly(og,sm.to_sympy(x),sm.to_sympy(y))
    args = (f,x,y) if name == 'discriminant' else (f,g,x,y)
    oracle_args = (of,sm.to_sympy(x),sm.to_sympy(y)) if name == 'discriminant' else (of,og,sm.to_sympy(x),sm.to_sympy(y))
    compare(getattr(sm,name)(*args,polys=polys),getattr(sy,name)(*oracle_args,polys=polys))
    for value in [False,True,None,2,'bad']:
        with pytest.raises(sm.FlagError):getattr(sm,name)(*args,auto=value)
    with pytest.raises(sm.OptionError):getattr(sm,name)(*args,bogus=True)


@pytest.mark.parametrize('name', ['resultant', 'discriminant', 'subresultants'])
@pytest.mark.parametrize('option', ['modulus','field','extension','domain','gens'])
def test_construction_options_are_applied(name, option):
    import sympy as sy
    x, y = sm.symbols('x y')
    f, g = x*x+y*x+1, 2*x+y
    opts = {'modulus':5} if option=='modulus' else {'field':True} if option=='field' else {'extension':sm.sqrt(2)} if option=='extension' else {'domain':sm.QQ} if option=='domain' else {'gens':(y,x)}
    oracle_opts = {'modulus':5} if option=='modulus' else {'field':True} if option=='field' else {'extension':sy.sqrt(2)} if option=='extension' else {'domain':sy.QQ} if option=='domain' else {'gens':(sm.to_sympy(y),sm.to_sympy(x))}
    generators = () if option=='gens' else (x,y)
    args = (f,*generators) if name=='discriminant' else (f,g,*generators)
    oracle_args = tuple(sm.to_sympy(v) for v in args)
    expected = getattr(sy,name)(*oracle_args,polys=True,**oracle_opts)
    compare(getattr(sm,name)(*args,polys=True,**opts),expected)


def test_aliases_assumptions_and_unused_parameters_remain_exact():
    x,y = sm.symbols('x y')
    a,b = sm.Dummy('a',positive=True), sm.Symbol('unused',real=True)
    for domain in [sm.QQ.algebraic_field(sm.sqrt(2),alias='alpha'),sm.QQ.frac_field(a,b,order='grevlex'),sm.ZZ.poly_ring(a,b,order='grlex')]:
        p,q = sm.Poly(x*x+y*x+1,x,y,domain=domain),sm.Poly(2*x+y,x,y,domain=domain)
        for result in p.subresultants(q)+[p.resultant(q),p.discriminant()]:
            assert result.domain == domain
            if domain.is_AlgebraicField:
                assert result.domain.ext.alias == sm.Symbol('alpha')
            else:
                assert result.domain.symbols == (a,b)
            assert pickle.loads(pickle.dumps(result)).domain == domain


def test_sparse_giant_degrees_and_growth_budget():
    x = sm.Symbol('x')
    p = sm.Poly(x**(2**32-1),x)
    assert p.resultant(sm.Poly(1,x)) == 1
    assert p.resultant(sm.Poly(x,x)) == 0
    with pytest.raises(ValueError,match='resource limit'):
        p.resultant(sm.Poly(2*x+1,x))


def test_native_prs_without_sympy_runtime():
    code = '''
import importlib.abc,sys,pickle
class Block(importlib.abc.MetaPathFinder):
 def find_spec(self,name,path=None,target=None):
  if name=='sympy' or name.startswith('sympy.'):raise ImportError('blocked')
sys.meta_path.insert(0,Block())
import symox as s
x,y=s.symbols('x y');a=s.Dummy('a',positive=True)
for k in [s.ZZ,s.QQ,s.GF(5),s.QQ.algebraic_field(s.sqrt(2)),s.QQ.frac_field(a),s.ZZ.poly_ring(a)]:
 p=s.Poly(x*x+y*x+1,x,y,domain=k);q=s.Poly(2*x+y,x,y,domain=k)
 res,prs=s.resultant(p,q,includePRS=True)
 assert res.gens==(y,) and res.domain==k
 assert prs==p.subresultants(q)==s.subresultants(p,q)
 assert all(v.domain==k for v in prs)
 assert s.Poly(res.as_expr()+y*y-4,y,domain=k).is_zero
 assert s.Poly(s.discriminant(p).as_expr()-y*y+4,y,domain=k).is_zero
 for v in [res,*prs]:assert pickle.loads(pickle.dumps(v)).domain==k
assert s.Poly(x+1,x).resultant(s.Poly(x**3+2,x))==1
assert s.Poly(x*x+1,x,modulus=5).discriminant()==1
assert not any(n=='sympy' or n.startswith('sympy.') for n in sys.modules)
'''
    result = subprocess.run([sys.executable,'-c',code],text=True,capture_output=True)
    assert result.returncode == 0,result.stderr


@pytest.mark.parametrize('generator', ['sqrt2', 'imaginary', 'cubic_quadratic', 'three_quadratics'])
def test_algebraic_embeddings_and_closed_scalar_results(generator):
    x = sm.Symbol('x')
    c = {'sqrt2':sm.sqrt(2), 'imaginary':sm.I,
         'cubic_quadratic':2**sm.Rational(1,3)+sm.I,
         'three_quadratics':sm.sqrt(2)+sm.sqrt(3)+sm.sqrt(5)}[generator]
    domain = sm.QQ.algebraic_field(c, alias='theta')
    p, q = sm.Poly(x*x+c*x+1, x, domain=domain), sm.Poly(2*x+1, x, domain=domain)
    assert domain.convert(p.resultant(q)) == domain.convert(5-2*c)
    assert domain.convert(p.discriminant()) == domain.convert(c*c-4)
    assert domain.convert(sm.resultant(p,q)) == domain.convert(5-2*c)
    assert domain.convert(sm.discriminant(p)) == domain.convert(c*c-4)
    for stage in p.subresultants(q):
        assert stage.domain == domain and stage.domain.ext.alias == sm.Symbol('theta')
        assert all(coefficient in domain for coefficient in stage.coeffs())


@pytest.mark.parametrize('name', ['resultant', 'discriminant', 'subresultants'])
@pytest.mark.parametrize('coefficient', ['radical', 'function', 'imaginary', 'float'])
def test_historical_opaque_coefficients_still_compute(name, coefficient):
    import sympy as sy
    x,a = sm.symbols('x a')
    c = {'radical':sm.sqrt(2), 'function':sm.sin(a), 'imaginary':sm.I, 'float':0.5*a}[coefficient]
    f,g = x*x+c*x+1, 2*x+1
    args = (f,x) if name == 'discriminant' else (f,g,x)
    compare(getattr(sm,name)(*args),getattr(sy,name)(*(sm.to_sympy(v) for v in args)))


def test_method_second_argument_keyword_matches_sympy():
    x = sm.Symbol('x')
    p,q = sm.Poly(x*x+1,x),sm.Poly(2*x+1,x)
    assert p.resultant(g=q,includePRS=True) == p.resultant(q,True)
    assert p.subresultants(g=q) == p.subresultants(q)
