"""Exact quadratic coefficient domains; SymPy is an oracle only in tests."""
import pickle
import pytest
import sympy as sp
import symox as sm


def pair_symbols():
    return sm.symbols('x y'), sp.symbols('x y')


def oracle(e):
    return sp.sympify(str(e), locals={'I': sp.I})


def same(a, b):
    assert sp.simplify(oracle(a) - b) == 0


@pytest.mark.parametrize('gen', ['sqrt2', 'I', 'sqrt3', 'scaled', 'negative', 'rational'])
@pytest.mark.parametrize('option', ['extension', 'list', 'domain', 'infer'])
def test_domain_construction_and_exact_arithmetic(gen, option):
    (x, y), (sx, sy) = pair_symbols()
    a, sa = {'sqrt2': (sm.sqrt(2), sp.sqrt(2)), 'I': (sm.I, sp.I),
             'sqrt3': (sm.sqrt(3), sp.sqrt(3)),
             'scaled': (2*sm.sqrt(2), 2*sp.sqrt(2)),
             'negative': (-sm.sqrt(2), -sp.sqrt(2)),
             'rational': (sm.sqrt(sm.Rational(2, 3)), sp.sqrt(sp.Rational(2, 3)))}[gen]
    opts = {'extension': a} if option == 'extension' else {'extension': [a]} if option == 'list' else {'domain': f'QQ<{a}>'} if option == 'domain' else {'extension': True}
    p = sm.Poly((1+a)*x*x + x/(1+a) + a, x, **opts)
    q = sp.Poly((1+sa)*sx*sx + sx/(1+sa) + sa, sx, extension=sa)
    assert p.domain.startswith('QQ<')
    same(p.as_expr(), q.as_expr())
    same(p.monic().as_expr(), q.monic().as_expr())
    same((p*p-p).as_expr(), (q*q-q).as_expr())
    same(p.diff().as_expr(), q.diff().as_expr())
    same(p.integrate().as_expr(), q.integrate().as_expr())
    assert p.content() == 1
    assert p.clear_denoms()[0] == 1
    d, r = divmod(p*p+1, p)
    same(d.as_expr(), q.as_expr())
    assert r.as_expr() == 1


@pytest.mark.parametrize('kind', ['sqrt2', 'I'])
@pytest.mark.parametrize('order', ['lex', 'grlex', 'grevlex'])
@pytest.mark.parametrize('option', ['extension', 'domain', 'infer', 'poly'])
def test_groebner_domains_orders_and_membership(kind, order, option):
    (x, y), (sx, sy) = pair_symbols()
    a, sa = (sm.sqrt(2), sp.sqrt(2)) if kind == 'sqrt2' else (sm.I, sp.I)
    f = [x*y-a, y*y-x]
    sf = [sx*sy-sa, sy*sy-sx]
    opts = {'extension': a} if option == 'extension' else {'domain': f'QQ<{a}>'} if option == 'domain' else {'extension': True} if option == 'infer' else {}
    if option == 'poly':
        f = [sm.Poly(e, x, y, extension=a) for e in f]
    g = sm.groebner(f, x, y, order=order, **opts)
    expected = sp.groebner(sf, sx, sy, order=order, extension=sa)
    assert g.domain == str(expected.domain)
    assert len(g) == len(expected)
    for actual, wanted in zip(g.exprs, expected.exprs):
        same(actual, wanted)
    assert g.is_zero_dimensional
    assert g.contains((x+y)*(x*y-a))
    assert not g.contains(x)
    inp = x*x*y + a*y + 1
    quotients, remainder = g.reduce(inp)
    expr = lambda e: e.as_expr() if isinstance(e, sm.Poly) else e
    identity = sum((expr(q)*b for q, b in zip(quotients, g.exprs)), sm.Integer(0)) + expr(remainder)
    same(identity, oracle(inp))
    sq, sr = expected.reduce(sx*sx*sy + sa*sy + 1)
    for q, eq in zip(quotients, sq):
        same(expr(q), eq)
    same(expr(remainder), sr)


@pytest.mark.parametrize('order', ['lex', 'grlex', 'grevlex'])
@pytest.mark.parametrize('polys', [False, True])
def test_reduced_exact_quotient_identity(order, polys):
    (x, y), (sx, sy) = pair_symbols()
    a, sa = sm.sqrt(2), sp.sqrt(2)
    f = x*x + a*x*y + 3
    divisors = [(1+a)*x+y, y*y-a]
    q, r = sm.reduced(f, divisors, x, y, order=order, extension=a, polys=polys)
    sq, sr = sp.reduced(sx*sx + sa*sx*sy+3, [(1+sa)*sx+sy, sy*sy-sa], sx, sy, order=order, extension=sa)
    expr = lambda e: e.as_expr() if isinstance(e, sm.Poly) else e
    for p, expected in zip(q, sq):
        same(expr(p), expected)
        if polys: assert p.domain == 'QQ<sqrt(2)>'
    same(expr(r), sr)
    same(sum((expr(p)*d for p, d in zip(q, divisors)), sm.Integer(0))+expr(r), oracle(f))


@pytest.mark.parametrize('kind', ['sqrt2', 'I'])
@pytest.mark.parametrize('expression', ['split', 'irreducible', 'multiple', 'coefficient'])
def test_exact_factorization_and_squarefree(kind, expression):
    x = sm.Symbol('x'); sx = sp.Symbol('x')
    a, sa = (sm.sqrt(2), sp.sqrt(2)) if kind == 'sqrt2' else (sm.I, sp.I)
    f, sf = {'split': (x*x-a*a, sx*sx-sa*sa),
             'irreducible': (x*x+3, sx*sx+3),
             'multiple': ((x-a)**3*(x+1)**2, (sx-sa)**3*(sx+1)**2),
             'coefficient': ((1+a)*x*x+2, (1+sa)*sx*sx+2)}[expression]
    p, spoly = sm.Poly(f, x, extension=a), sp.Poly(sf, sx, extension=sa)
    c, fs = p.factor_list()
    sc, sfs = spoly.factor_list()
    same(c, sc)
    assert sorted((sp.degree(oracle(q.as_expr()), sx), k) for q, k in fs) == sorted((q.degree(), k) for q, k in sfs)
    product = c
    for q, k in fs:
        assert q.domain == p.domain
        product *= q.as_expr()**k
    same(product, sf)
    c, fs = p.sqf_list()
    product = c
    for q, k in fs: product *= q.as_expr()**k
    same(product, sf)
    assert p.is_irreducible == spoly.is_irreducible
    assert p.is_sqf == spoly.is_sqf


def test_gcd_compose_shift_eval_and_compatible_subfield():
    (x,y), (sx,sy) = pair_symbols(); a = sm.sqrt(2)
    p = sm.Poly((x-a)*(x+1), x, extension=a)
    q = sm.Poly((x-a)*(x+2), x, extension=a)
    same(p.gcd(q).as_expr(), sx-sp.sqrt(2))
    same(p.lcm(q).as_expr(), (sx-sp.sqrt(2))*(sx+1)*(sx+2))
    same(p.compose(sm.Poly(x*x+a, x, extension=a)).as_expr(), ((sx*sx+sp.sqrt(2))-sp.sqrt(2))*((sx*sx+sp.sqrt(2))+1))
    same(p.shift(a).as_expr(), sx*(sx+sp.sqrt(2)+1))
    assert p.eval(a) == 0
    r = sm.Poly(x*y+a*x, x, y, extension=a).eval(x,a)
    assert isinstance(r,sm.Poly) and r.domain == p.domain
    same(r.as_expr(), sp.sqrt(2)*sy+2)
    same((p+sm.Poly(x+sm.sqrt(8),x,extension=sm.sqrt(8))).as_expr(), oracle(p.as_expr())+sx+sp.sqrt(8))


@pytest.mark.parametrize('kind',['poly','basis','empty','unit'])
def test_algebraic_pickle_and_unused_generators(kind):
    (x,y), _ = pair_symbols(); a = sm.sqrt(2)
    value = sm.Poly(x+a,x,extension=a) if kind=='poly' else sm.groebner([] if kind=='empty' else [1] if kind=='unit' else [x-a,y-x],x,y,extension=a,polys=True)
    restored = pickle.loads(pickle.dumps(value))
    assert restored.domain == value.domain
    assert restored.gens == value.gens
    if kind=='poly': assert restored.as_expr() == value.as_expr()
    else:
        assert restored.exprs == value.exprs
        assert restored.is_zero_dimensional == value.is_zero_dimensional


@pytest.mark.parametrize('order',['lex','grlex','grevlex'])
def test_algebraic_fglm(order):
    (x,y), (sx,sy) = pair_symbols();a=sm.sqrt(2)
    g=sm.groebner([x*y-a,y*y-x],x,y,extension=a,order='grevlex')
    converted=g.fglm(order)
    expected=sp.groebner([sx*sy-sp.sqrt(2),sy*sy-sx],sx,sy,extension=sp.sqrt(2),order=order)
    assert converted.domain==g.domain
    for p,q in zip(converted.exprs,expected.exprs):same(p,q)


@pytest.mark.parametrize('call',[
    lambda x,a: sm.Poly(x+a,x,extension=False),
    lambda x,a: sm.Poly(x+a,x,extension=a,domain='QQ'),
    lambda x,a: sm.Poly(x+a,x,extension=a,modulus=5),
    lambda x,a: sm.Poly(x+sm.sqrt(3),x,extension=a),
    lambda x,a: sm.Poly(x+a,x,domain='QQ<sqrt(2)'),
    lambda x,a: sm.groebner([x+a],x,extension=a,modulus=5),
    lambda x,a: sm.reduced(x+a,[x],x,extension=a,modulus=5),
])
def test_unsupported_fields_fail_explicitly(call):
    with pytest.raises((ValueError,TypeError,NotImplementedError)):
        call(sm.Symbol('x'),sm.sqrt(2))


def test_plain_extension_inference_does_not_change_rational_domain():
    x=sm.Symbol('x')
    assert sm.Poly(x*x-2,x,extension=True).domain=='ZZ'
    assert sm.groebner([x*x-2],x,extension=True).domain=='ZZ'


def test_top_level_algebraic_factor_and_gcd_options():
    x=sm.Symbol('x');sx=sp.Symbol('x');a=sm.sqrt(2)
    f=(x-a)**2*(x+1)
    same(sm.factor(x*x-2,x,extension=a),sx*sx-2)
    c,fs=sm.factor_list(f,x,extension=a)
    product=c
    for q,k in fs:product*=q**k
    same(product,(sx-sp.sqrt(2))**2*(sx+1))
    same(sm.gcd(f,(x-a)*(x+3),x,extension=a),sx-sp.sqrt(2))
    p=sm.Poly(f,x,extension=a);q=sm.Poly((x-a)*(x+3),x,extension=a)
    assert sm.gcd(p,q).domain==p.domain
    same(sm.gcd(p,q).as_expr(),sx-sp.sqrt(2))
    same(p.discriminant(),0)
    same(p.resultant(q),0)
    same(sm.Poly(x*x+a*x+3,x,extension=a).discriminant(),-10)


def test_runtime_algebraic_fields_have_no_sympy_dependency():
    import subprocess
    import sys
    code = '''
import builtins
original=builtins.__import__
def guarded(name,*args,**kwargs):
    if name=='sympy' or name.startswith('sympy.'):
        raise AssertionError('runtime SymPy import')
    return original(name,*args,**kwargs)
builtins.__import__=guarded
import symox as s
x,y=s.symbols('x y');a=s.sqrt(2)
p=s.Poly((x-a)**2,x,extension=a)
assert p.sqf_list()[1][0][1]==2
assert len(s.Poly(x*x-2,x,extension=a).factor_list()[1])==2
g=s.groebner([x*y-a,y*y-x],x,y,extension=a,order='grevlex')
assert g.contains((x+y)*(x*y-a))
assert g.fglm('lex').domain=='QQ<sqrt(2)>'
'''
    result=subprocess.run([sys.executable,'-c',code],text=True,capture_output=True)
    assert result.returncode==0,result.stderr


@pytest.mark.parametrize('operation', [lambda p,q:p+q,lambda p,q:p-q,lambda p,q:p*q,lambda p,q:p.div(q)])
def test_independent_number_fields_do_not_silently_degrade_to_expressions(operation):
    x=sm.Symbol('x')
    p=sm.Poly(x+sm.sqrt(2)+sm.sqrt(3),x,extension=[sm.sqrt(2),sm.sqrt(3)])
    b=sm.Integer(2)**sm.Rational(1,3)
    q=sm.Poly(x+b,x,extension=b)
    with pytest.raises((TypeError,ValueError,NotImplementedError)):
        operation(p,q)


@pytest.mark.parametrize('value',[0,1,2])
def test_algebraic_constant_properties(value):
    x=sm.Symbol('x');p=sm.Poly(value,x,extension=sm.sqrt(2))
    q=sp.Poly(value,sp.Symbol('x'),extension=sp.sqrt(2))
    assert p.is_irreducible==q.is_irreducible
    assert p.is_sqf==q.is_sqf
    same(p.sqf_part().as_expr(),q.sqf_part().as_expr())


@pytest.mark.parametrize('base',[sm.Integer(4),sm.Rational(4,9)])
def test_reducible_held_quadratic_generators_are_rejected(base):
    x=sm.Symbol('x');held=sm.Pow(base,sm.Rational(1,2),evaluate=False)
    with pytest.raises((ValueError,TypeError,NotImplementedError)):
        sm.Poly(x+held,x,extension=held)
