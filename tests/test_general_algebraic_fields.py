"""Certified cubic/biquadratic fields, proper subfields and primitive maps."""
import pickle
import subprocess
import sys

import pytest
import sympy as sp
import symox as sm

GENERATORS = ['(2/3)**(1/3)', '(5/7)**(1/3)', '2**(1/3)', '1+2**(1/3)', '1/(1+2**(1/3))',
              'sqrt(2)+I', 'sqrt(2)+sqrt(3)', '(sqrt(2)+I)**2',
              '1+sqrt(2)', '1/(sqrt(2)+I)']


def expr(text,module=sm):
    return module.sympify(text)


def oracle(value):
    if isinstance(value,sm.Poly):value=value.as_expr()
    return sp.sympify(str(value),locals={'I':sp.I})


def same(actual,expected,*generators):
    delta=oracle(actual)-expected
    if generators:assert sp.Poly(delta,*generators,extension=True).is_zero
    else:assert sp.simplify(delta)==0


@pytest.mark.parametrize('generator',GENERATORS)
@pytest.mark.parametrize('option',['extension','domain'])
def test_exact_field_arithmetic_for_transformed_generators(generator,option):
    x,sx=sm.Symbol('x'),sp.Symbol('x');a,sa=expr(generator),expr(generator,sp)
    options={'extension':a} if option=='extension' else {'domain':f'QQ<{a}>'}
    p=sm.Poly((1+a)*x*x+(a*a-a)*x+1/(a+2),x,**options)
    q=sp.Poly((1+sa)*sx*sx+(sa*sa-sa)*sx+1/(sa+2),sx,extension=sa)
    assert p.domain.startswith('QQ<')
    same(p,q.as_expr(),sx)
    same(p.monic(),q.monic().as_expr(),sx)
    same((p*p-p), (q*q-q).as_expr(),sx)
    quotient,remainder=(p*p+1).div(p)
    same(quotient,q.as_expr(),sx);assert remainder.as_expr()==1
    same(p.diff(),q.diff().as_expr(),sx)
    same(p.integrate(),q.integrate().as_expr(),sx)
    same(p.shift(a),q.shift(sa).as_expr(),sx)
    same(p.eval(a),q.eval(sa))
    assert pickle.loads(pickle.dumps(p)).as_expr()==p.as_expr()


@pytest.mark.parametrize('generator',GENERATORS+['2/3','0'])
@pytest.mark.parametrize('polys',[False,True])
def test_minimal_polynomial_is_irreducible_and_matches_oracle(generator,polys):
    x,sx=sm.Symbol('t'),sp.Symbol('t');a,sa=expr(generator),expr(generator,sp)
    actual=sm.minimal_polynomial(a,x,polys=polys)
    wanted=sp.minimal_polynomial(sa,sx)
    assert sp.expand(oracle(actual)-wanted)==0
    assert sp.Poly(oracle(actual),sx).is_irreducible
    assert sp.simplify(oracle(actual).subs(sx,sa))==0
    assert sm.minpoly(a,x,polys=polys)==actual
    if polys:assert actual.domain=='QQ' and actual.gens==(x,)


@pytest.mark.parametrize('generators',[
    ['sqrt(2)','I'],['sqrt(2)','sqrt(3)'],['sqrt(2)','sqrt(8)'],
    ['2**(1/3)','4**(1/3)'],['1/2','sqrt(2)'],
    ['sqrt(2)','sqrt(3)','sqrt(6)'],
    ['sqrt(2)+sqrt(3)','sqrt(2)-sqrt(3)'],
])
@pytest.mark.parametrize('polys',[False,True])
def test_primitive_element_linear_combination_and_exact_embeddings(generators,polys):
    x,sx=sm.Symbol('t'),sp.Symbol('t');values=[expr(g) for g in generators]
    minimal,weights,reps=sm.primitive_element(values,x,ex=True,polys=polys)
    assert len(weights)==len(values)==len(reps)
    primitive=sum((c*a for c,a in zip(weights,values)),sm.Integer(0))
    wanted=sp.minimal_polynomial(oracle(primitive),sx)
    assert sp.expand(oracle(minimal)-wanted)==0
    for original,coefficients in zip(values,reps):
        reconstructed=sm.Integer(0)
        for coefficient in coefficients:reconstructed=reconstructed*primitive+coefficient
        same(reconstructed,oracle(original))
    plain=sm.primitive_element(values,x,ex=False,polys=polys)
    assert plain==(minimal,weights)
    if polys:assert minimal.domain=='QQ'


@pytest.mark.parametrize('generator',['2**(1/3)','sqrt(2)+I','sqrt(2)+sqrt(3)','1+sqrt(2)','(sqrt(2)+I)**2'])
def test_norms_and_true_factorization_in_supported_fields(generator):
    x,sx=sm.Symbol('x'),sp.Symbol('x');a,sa=expr(generator),expr(generator,sp)
    p=sm.Poly((x-a)/(1+a),x,extension=a)
    expected=sp.Poly((sx-sa)/(1+sa),sx,extension=sa)
    norm=p.norm();assert norm.domain=='QQ'
    assert sp.expand(oracle(norm)-expected.norm().as_expr())==0
    minimal=sm.minimal_polynomial(a,x)
    target=sm.Poly(minimal,x,extension=a)
    c,factors=target.factor_list();sc,sf=sp.Poly(sp.minpoly(sa,sx),sx,extension=sa).factor_list()
    same(c,sc)
    assert sorted((f.degree(),k) for f,k in factors)==sorted((f.degree(),k) for f,k in sf)
    product=sm.Poly(c,x,extension=a)
    for f,k in factors:
        assert f.domain==target.domain
        assert target.rem(f).is_zero
        product*=f**k
    assert product==target
    repeated=sm.Poly((x-a)**3*(x+1)**2,x,extension=a)
    c,parts=repeated.sqf_list();product=sm.Poly(c,x,extension=a)
    for f,k in parts:product*=f**k
    assert product==repeated
    same(repeated.gcd(sm.Poly((x-a)*(x+3),x,extension=a)),sx-sa,sx)


@pytest.mark.parametrize('generator',['2**(1/3)','sqrt(2)+I','sqrt(2)+sqrt(3)','(sqrt(2)+I)**2'])
@pytest.mark.parametrize('order',['lex','grlex','grevlex'])
@pytest.mark.parametrize('method',['buchberger','f5b'])
def test_groebner_reduced_and_fglm_use_actual_coefficient_fields(generator,order,method):
    (x,y),(sx,sy)=sm.symbols('x y'),sp.symbols('x y');a,sa=expr(generator),expr(generator,sp)
    fs=[a*x*y-1,x*x-y]
    g=sm.groebner(fs,x,y,extension=a,order=order,method=method,polys=True)
    expected=sp.groebner([sa*sx*sy-1,sx*sx-sy],sx,sy,extension=sa,order=order,method=method)
    assert len(g)==len(expected) and g.is_zero_dimensional
    for actual,wanted in zip(g.exprs,expected.exprs):same(actual,wanted,sx,sy)
    assert all(p.domain==g.domain for p in g)
    assert g.contains((x+y)*(a*x*y-1))
    q,r=g.reduce(x*x*y+a*y+1)
    same(sum((p.as_expr()*f for p,f in zip(q,g.exprs)),sm.Integer(0))+r.as_expr(),sx*sx*sy+sa*sy+1,sx,sy)
    lex=g.fglm('lex');lex_expected=expected.fglm('lex')
    for actual,wanted in zip(lex.exprs,lex_expected.exprs):same(actual,wanted,sx,sy)
    quotients,remainder=sm.reduced(x*x+a*y,[a*x+y,y*y-a],x,y,order=order,extension=a,polys=True)
    sq,sr=sp.reduced(sx*sx+sa*sy,[sa*sx+sy,sy*sy-sa],sx,sy,order=order,extension=sa)
    for actual,wanted in zip(quotients,sq):same(actual,wanted,sx,sy)
    same(remainder,sr,sx,sy)
    restored=pickle.loads(pickle.dumps(g));assert restored.domain==g.domain and restored.exprs==g.exprs


def test_joint_coercion_inferred_fields_and_dependent_extensions():
    x,y=sm.symbols('x y');sx,sy=sp.symbols('x y');u,v=sm.sqrt(2),sm.I
    p,q=sm.Poly(x+u,x,extension=u),sm.Poly(x+v,x,extension=v)
    for actual,wanted in [(p+q,2*sx+sp.sqrt(2)+sp.I),(p*q,(sx+sp.sqrt(2))*(sx+sp.I)),(p-q,sp.sqrt(2)-sp.I)]:
        assert isinstance(actual,sm.Poly) and actual.domain.startswith('QQ<')
        same(actual,wanted,sx)
    joint=sm.Poly(x+u+v,x,extension=[u,v])
    same(joint,sx+sp.sqrt(2)+sp.I,sx)
    inferred=sm.Poly(x+u+v,x,extension=True);same(inferred,sx+sp.sqrt(2)+sp.I,sx)
    g=sm.groebner([sm.Poly(x-u,x,y,extension=u),sm.Poly(y-v,x,y,extension=v)],x,y,method='f5b')
    assert g.contains(x-u) and g.contains(y-v)
    assert sm.Poly(x+sm.sqrt(8),x,extension=[u,sm.sqrt(8)]).domain=='QQ<sqrt(2)>'


def test_subfield_membership_is_not_guessed_from_an_annihilating_polynomial():
    x=sm.Symbol('x');theta=(sm.sqrt(2)+sm.I)**2
    assert sp.degree(oracle(sm.minpoly(theta,x)),sp.Symbol('x'))==2
    with pytest.raises(ValueError):sm.Poly(x+sm.I,x,extension=theta)
    with pytest.raises(ValueError):sm.Poly(x+sm.sqrt(2),x,extension=theta)
    p=sm.Poly(x+sm.I*sm.sqrt(2),x,extension=theta)
    same(p,sp.Symbol('x')+sp.I*sp.sqrt(2),sp.Symbol('x'))


@pytest.mark.parametrize('generator',[
    '2**(1/5)','(-2)**(1/3)','sqrt(1+sqrt(2))',
    '2**(1/3)+sqrt(2)','sqrt(2)+sqrt(3)+sqrt(5)','pi',
])
def test_outside_supported_towers_fails_explicitly(generator):
    x=sm.Symbol('x');a=expr(generator)
    with pytest.raises((ValueError,NotImplementedError)):sm.Poly(x+a,x,extension=a)
    with pytest.raises((ValueError,NotImplementedError)):sm.minpoly(a,x)


def test_numberfield_import_paths_options_and_no_sympy_runtime():
    from symox.polys.numberfields.minpoly import minimal_polynomial
    from symox.polys.numberfields.subfield import primitive_element
    from symox.polys.polytools import Poly
    x=sm.Symbol('x');a=2**sm.Rational(1,3)
    assert minimal_polynomial(a,x)==x**3-2
    assert primitive_element([sm.sqrt(2),sm.I],x)[0]==x**4-2*x*x+9
    assert Poly is sm.Poly
    with pytest.raises(NotImplementedError):sm.minpoly(a,x,compose=False)
    with pytest.raises(NotImplementedError):sm.minpoly(a,x,domain='QQ(a)')
    with pytest.raises(ValueError):sm.minpoly(a,1)
    result=subprocess.run([sys.executable,'-c','''
import sys
sys.modules['sympy']=None
import symox as s
from symox.polys.numberfields import minimal_polynomial,primitive_element
x,y=s.symbols('x y');a=2**s.Rational(1,3)
assert minimal_polynomial(a,x)==x**3-2
assert len(s.Poly(x*x*x-2,x,extension=a).factor_list()[1])==2
b=s.sqrt(2)+s.I
assert s.Poly(x-b,x,extension=b).norm().as_expr()==x**4-2*x*x+9
p=s.groebner([a*x*y-1,x*x-y],x,y,extension=a,method='f5b',order='grevlex')
assert p.fglm('lex').contains(a*x*y-1)
assert len(primitive_element([s.sqrt(2),s.I],x,ex=True))==3
'''],capture_output=True,text=True)
    assert result.returncode==0,result.stderr


@pytest.mark.parametrize('a,b,c,divisor',[(2,3,6,1),(6,15,10,3)])
@pytest.mark.parametrize('order',['lex','grlex','grevlex'])
def test_biquadratic_field_closure_under_canonical_products(a,b,c,divisor,order):
    x=sm.Symbol('x');u,v,target=sm.sqrt(a),sm.sqrt(b),sm.sqrt(c)
    options={'extension':[u,v]}
    p=sm.Poly(x+target,x,**options)
    assert sm.Poly(p.as_expr(),x,domain=p.domain)==p
    q=sm.Poly(x+u*v/divisor,x,**options)
    assert p==q
    g=sm.groebner([x-target],x,order=order,method='f5b',**options)
    assert g.contains(x*x-c)
    assert g.reduce(x*x-c)[1]==0
