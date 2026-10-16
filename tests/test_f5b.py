"""Differential and exact certificates for the native signature F5B engine."""
import itertools
import pickle
import subprocess
import sys

import pytest
import sympy as sp
import symox as sm

ORDERS = ['lex', 'grlex', 'grevlex']


def oracle(expr):
    if isinstance(expr, sm.Poly):
        expr = expr.as_expr()
    return sp.sympify(str(expr), locals={'I': sp.I})


def expressions(strings, module):
    symbols = module.symbols('x y z a b')
    env = dict(zip(['x', 'y', 'z', 'a', 'b'], symbols))
    env.update(sqrt=module.sqrt, I=module.I)
    return [eval(text, {'__builtins__': {}}, env) for text in strings], symbols


@pytest.mark.parametrize('order', ORDERS)
@pytest.mark.parametrize('system', [
    ['x**2+2*x*y**2', 'x*y+2*y**3-1'],
    ['x**2', 'x*y', 'y**2'],
    ['x**2-y', 'x**3-z'],
    ['x*y-1', 'y*z-1', 'z*x-1'],
    ['x**2-y', 'y**2-z', 'z**2-2'],
    ['x*y', 'x*y+1'],
    ['0', 'x+y', 'x+y', '2*x+2*y'],
    ['0', '0'], [], ['1'],
    ['x**2/2+y/3', 'x*y-1/5'],
])
def test_f5b_rational_bases_match_sympy_and_buchberger(system, order):
    fs, (x,y,z,a,b) = expressions(system, sm)
    sf, (sx,sy,sz,sa,sb) = expressions(system, sp)
    # Keep exact rational coefficients when Python source contains division
    # between literal integers.
    fs = [sm.sympify(text) for text in system]
    sf = [sp.sympify(text) for text in system]
    g = sm.groebner(fs, x, y, z, order=order, method='f5b', domain='QQ')
    expected = sp.groebner(sf, sx, sy, sz, order=order, method='f5b', domain='QQ')
    assert g.exprs == sm.groebner(fs, x, y, z, order=order, method='buchberger', domain='QQ').exprs
    assert len(g) == len(expected)
    for actual, wanted in zip(g.exprs, expected.exprs):
        assert sp.Poly(oracle(actual)-wanted, sx,sy,sz).is_zero
    assert all(g.contains(f) for f in fs)


@pytest.mark.parametrize('order', ORDERS)
@pytest.mark.parametrize('modulus', [2,3,5,101,4294967291])
@pytest.mark.parametrize('kind', ['expression', 'poly'])
def test_f5b_prime_fields_including_characteristic_two(order, modulus, kind):
    (x,y),(sx,sy) = sm.symbols('x y'), sp.symbols('x y')
    fs = [2*x*x+y-1, x*y+3*x+2, y*y+x]
    sf = [2*sx*sx+sy-1,sx*sy+3*sx+2,sy*sy+sx]
    if kind=='poly':
        fs = [sm.Poly(f,x,y,modulus=modulus) for f in fs]
    opts = {'modulus':modulus}
    g = sm.groebner(fs,x,y,order=order,method='f5b',**opts)
    expected = sp.groebner(sf,sx,sy,order=order,method='f5b',modulus=modulus)
    assert g.domain == f'GF({modulus})'
    assert len(g) == len(expected)
    for actual,wanted in zip(g.exprs,expected.exprs):
        assert sp.Poly(oracle(actual)-wanted,sx,sy,modulus=modulus).is_zero
    assert all(g.contains(f) for f in fs)
    assert g.exprs == sm.groebner(fs,x,y,order=order,method='buchberger',**opts).exprs


@pytest.mark.parametrize('order',ORDERS)
@pytest.mark.parametrize('domain',['QQ(a,b)','ZZ(a,b)','QQ[a,b]','ZZ[a,b]'])
def test_f5b_parameter_coefficient_fields_and_rings(order,domain):
    (x,y,a,b),(sx,sy,sa,sb) = sm.symbols('x y a b'),sp.symbols('x y a b')
    fs = [a*x*y-1,x*x-y*(b+1)]
    sf = [sa*sx*sy-1,sx*sx-sy*(sb+1)]
    g = sm.groebner(fs,x,y,order=order,method='f5b',domain=domain)
    expected = sp.groebner(sf,sx,sy,order=order,method='f5b',domain=domain)
    assert g.domain==str(expected.domain)
    assert len(g)==len(expected)
    for actual,wanted in zip(g.exprs,expected.exprs):
        assert sp.cancel(oracle(actual)-wanted)==0
    assert g.exprs==sm.groebner(fs,x,y,order=order,domain=domain).exprs


@pytest.mark.parametrize('order',ORDERS)
@pytest.mark.parametrize('extension',['sqrt2','I'])
@pytest.mark.parametrize('poly_input',[False,True])
def test_f5b_algebraic_fields(order,extension,poly_input):
    (x,y),(sx,sy)=sm.symbols('x y'),sp.symbols('x y')
    a,sa=(sm.sqrt(2),sp.sqrt(2)) if extension=='sqrt2' else (sm.I,sp.I)
    fs=[a*x*y-1,x*x-y]
    if poly_input: fs=[sm.Poly(f,x,y,extension=a) for f in fs]
    kwargs={} if poly_input else {'extension':a}
    g=sm.groebner(fs,x,y,order=order,method='f5b',**kwargs)
    expected=sp.groebner([sa*sx*sy-1,sx*sx-sy],sx,sy,order=order,extension=sa,method='f5b')
    assert g.domain==str(expected.domain)
    assert len(g)==len(expected)
    for actual,wanted in zip(g.exprs,expected.exprs): assert sp.simplify(oracle(actual)-wanted)==0
    assert g.exprs==sm.groebner(fs,x,y,order=order,**kwargs).exprs
    assert g.contains(a*x*y-1)


@pytest.mark.parametrize('order',ORDERS)
@pytest.mark.parametrize('gens_reversed',[False,True])
def test_f5b_deterministic_permutations_and_zz_normalization(order,gens_reversed):
    (x,y,z),(sx,sy,sz)=sm.symbols('x y z'),sp.symbols('x y z')
    fs=[2*x*x+3*y,3*x*y+z,x*z-1]
    sf=[2*sx*sx+3*sy,3*sx*sy+sz,sx*sz-1]
    gens=(z,y,x) if gens_reversed else (x,y,z)
    sgens=(sz,sy,sx) if gens_reversed else (sx,sy,sz)
    expected=sp.groebner(sf,*sgens,order=order,method='f5b',domain='ZZ')
    for perm in itertools.permutations(fs):
        g=sm.groebner(perm,*gens,order=order,method='f5b',domain='ZZ')
        assert g.gens==gens and g.domain=='ZZ'
        assert len(g)==len(expected)
        for actual,wanted in zip(g.exprs,expected.exprs): assert sp.expand(oracle(actual)-wanted)==0


@pytest.mark.parametrize('order',ORDERS)
def test_f5b_cyclic_four_benchmark_and_s_pairs(order):
    a,b,c,d=sm.symbols('a b c d');sa,sb,sc,sd=sp.symbols('a b c d')
    fs=[a+b+c+d,a*b+b*c+c*d+d*a,a*b*c+a*b*d+a*c*d+b*c*d,a*b*c*d-1]
    sf=[sa+sb+sc+sd,sa*sb+sb*sc+sc*sd+sd*sa,sa*sb*sc+sa*sb*sd+sa*sc*sd+sb*sc*sd,sa*sb*sc*sd-1]
    g=sm.groebner(fs,a,b,c,d,order=order,method='f5b',domain='QQ')
    expected=sp.groebner(sf,sa,sb,sc,sd,order=order,method='f5b',domain='QQ')
    assert [oracle(e) for e in g.exprs]==expected.exprs
    # Construct every S-pair using the requested order and independently
    # check exact ideal reduction, including grevlex non-storage leaders.
    terms=[sp.Poly(oracle(e),sa,sb,sc,sd) for e in g.exprs]
    native_gens=(a,b,c,d)
    for i,p in enumerate(terms):
        for q in terms[i+1:]:
            ep,cp=p.terms(order=order)[0];eq,cq=q.terms(order=order)[0]
            lcm=tuple(max(u,v) for u,v in zip(ep,eq))
            mp=sm.Integer(1);mq=sm.Integer(1)
            for x,k,u,v in zip(native_gens,lcm,ep,eq):mp*=x**(k-u);mq*=x**(k-v)
            spoly=mp*g.exprs[i]/sm.sympify(str(cp))-mq*sm.sympify(str(q.as_expr()))/sm.sympify(str(cq))
            assert g.reduce(spoly)[1]==0


def test_f5b_api_validation_reduction_fglm_pickle_and_no_runtime_sympy():
    x,y=sm.symbols('x y')
    g=sm.GroebnerBasis([x*y-1,x*x-y],x,y,method='f5b',polys=True,order='grevlex')
    assert all(isinstance(p,sm.Poly) for p in g)
    assert sm.is_zero_dimensional([x*y-1,x*x-y],x,y,method='f5b')
    converted=g.fglm('lex');assert converted.exprs==sm.groebner([x*y-1,x*x-y],x,y).exprs
    restored=pickle.loads(pickle.dumps(g));assert restored.exprs==g.exprs
    assert restored.domain==g.domain and restored.order==g.order
    with pytest.raises(ValueError,match='valid Groebner algorithm'):sm.groebner([x],x,method='f4')
    with pytest.raises(ValueError,match='denominator'):sm.groebner([x/5],x,method='f5b',modulus=5)
    with pytest.raises(ValueError,match='prime'):sm.groebner([x],x,method='f5b',modulus=4)
    result=subprocess.run([sys.executable,'-c', '''
import sys
sys.modules['sympy']=None
import symox as s
x,y,a=s.symbols('x y a')
for options in ({'domain':'QQ'},{'modulus':2},{'modulus':101},{'domain':'QQ(a)'},{'extension':s.sqrt(2)}):
    c=a if 'a' in options.get('domain','') else s.sqrt(2) if 'extension' in options else s.Integer(1)
    g=s.groebner([c*x*y-1,x*x-y],x,y,method='f5b',**options)
    assert g.contains(c*x*y-1)
'''],capture_output=True,text=True)
    assert result.returncode==0,result.stderr


@pytest.mark.parametrize('method',['buchberger','f5b'])
@pytest.mark.parametrize('order',ORDERS)
def test_input_poly_keeps_fraction_domain_and_unused_parameters(method,order):
    x,y=sm.symbols('x y');a=sm.Symbol('a',nonzero=True);b=sm.Symbol('b',real=True)
    fs=[sm.Poly(a*x*y-1,x,y,domain='QQ(a,b)'),sm.Poly(x*x-y,x,y,domain='QQ(a,b)')]
    g=sm.groebner(fs,x,y,order=order,method=method)
    assert g.domain=='QQ(a,b)'
    assert all(isinstance(p,sm.Poly) and p.domain==g.domain for p in g)
    assert g.contains(a*x*y-1)
    assert g.exprs==sm.groebner([p.as_expr() for p in fs],x,y,domain='QQ(a,b)',order=order,method=method).exprs


@pytest.mark.parametrize('method',['buchberger','f5b'])
def test_actual_parameter_assumption_conflicts_remain_explicit(method):
    x=sm.Symbol('x');a=sm.Symbol('a',nonzero=True);other_a=sm.Symbol('a')
    first=sm.Poly(a*x-1,x,domain='QQ(a)')
    second=sm.Poly(other_a*x+1,x,domain='QQ(a)')
    with pytest.raises(ValueError,match='ambiguous parameter assumptions'):
        sm.groebner([first,second],x,method=method)


@pytest.mark.parametrize('method',['buchberger','f5b'])
def test_finite_poly_input_without_modulus_uses_sympy_domain_inference(method):
    x=sm.Symbol('x');sx=sp.Symbol('x')
    g=sm.groebner([sm.Poly(2*x+1,x,modulus=5)],x,method=method)
    expected=sp.groebner([sp.Poly(2*sx+1,sx,modulus=5)],sx,method=method)
    assert g.domain==str(expected.domain)=='ZZ'
    assert [oracle(e) for e in g.exprs]==expected.exprs
