"""Exact algebraic extension families, including closure and ideal certificates."""
import pickle

import pytest
import sympy as sp
import symox as sm

FIELDS = ['2**(1/3)+sqrt(2)', '2**(1/3)+I', '2**(1/3)+I*sqrt(3)',
          '(2/3)**(1/3)+sqrt(5)', 'sqrt(2)+sqrt(3)+sqrt(5)',
          'sqrt(2)+sqrt(3)+I', '2**(1/6)', '1/(1+2**(1/3)+sqrt(2))']


def same(a,b,*gens):
    delta=sp.sympify(str(a))-sp.sympify(str(b))
    if gens:
        assert sp.Poly(delta,*[sp.Symbol(str(g)) for g in gens],extension=True).is_zero,(a,b)
    else:
        assert sp.expand(delta)==0 or sp.simplify(delta)==0,(a,b)


@pytest.mark.parametrize('field',FIELDS)
@pytest.mark.parametrize('mode',['explicit','inferred','domain'])
def test_domains_minpoly_arithmetic_and_exact_norm(field,mode):
    x=sm.Symbol('x');sx=sp.Symbol('x');theta=sm.sympify(field);eta=sp.sympify(field)
    kwargs={'extension':theta} if mode=='explicit' else {'extension':True} if mode=='inferred' else {'domain':f'QQ<{theta}>'}
    p=sm.Poly(x-theta,x,**kwargs)
    oracle=sp.Poly(sx-eta,sx,extension=eta)
    assert p.domain.startswith('QQ<')
    assert sm.minpoly(theta,x,polys=True).degree() == sp.degree(sp.minimal_polynomial(eta,sx))
    # SymPy's integral defining polynomial can be nonmonic; its norm()
    # retains that resultant scaling. Compare the mathematical field norm,
    # whose multiplication determinant is monic for x-theta.
    expected_norm=sp.Poly(oracle.norm().as_expr(),sx).monic().as_expr()
    same(p.norm().as_expr(),expected_norm,x)
    same(p.norm().as_expr(),sp.Poly(sp.minimal_polynomial(eta,sx),sx).monic().as_expr(),x)
    r=p*p+1
    quotient,remainder=r.div(p)
    same((quotient*p+remainder).as_expr(),r.as_expr(),x)
    assert remainder.as_expr()==1
    restored=pickle.loads(pickle.dumps(p))
    assert restored.domain==p.domain and restored.as_expr()==p.as_expr()
    if mode=='explicit':
        f=sm.Poly((x-theta)*(x-1),x,**kwargs)
        c,factors=f.factor_list()
        product=sm.Poly(c,x,**kwargs)
        for q,k in factors:
            assert f.rem(q).is_zero
            product*=q**k
        same(product.as_expr(),f.as_expr(),x)


@pytest.mark.parametrize('field',FIELDS[:3]+FIELDS[4:7])
@pytest.mark.parametrize('order',['lex','grlex','grevlex'])
@pytest.mark.parametrize('method',['buchberger','f5b'])
def test_groebner_reduction_and_fglm_over_new_fields(field,order,method):
    x,y=sm.symbols('x y');sx,sy=sp.symbols('x y');a=sm.sympify(field);sa=sp.sympify(field)
    inputs=[a*x*y-1,x*x-y]
    basis=sm.groebner(inputs,x,y,extension=a,order=order,method=method)
    oracle=sp.groebner([sa*sx*sy-1,sx*sx-sy],sx,sy,extension=sa,order=order,method=method)
    assert len(basis)==len(oracle)
    for actual,expected in zip(basis.exprs,oracle.exprs):same(actual,expected,x,y)
    q,r=basis.reduce(x**4+y+1)
    identity=r+sum((coefficient*polynomial for coefficient,polynomial in zip(q,basis.exprs)),sm.Integer(0))
    same(identity,x**4+y+1,x,y)
    assert all(basis.contains(e) for e in inputs)
    converted=basis.fglm('lex')
    expected=oracle.fglm('lex')
    for actual,target in zip(converted.exprs,expected.exprs):same(actual,target,x,y)
    q,r=sm.reduced(x**4+y+1,basis.exprs,x,y,extension=a,order=order)
    same(r,basis.reduce(x**4+y+1)[1],x,y)


@pytest.mark.parametrize('items,member,expected_degree',[
    (['2**(1/3)','sqrt(2)','4**(1/3)','sqrt(8)'],'2**(7/6)',6),
    (['2**(1/3)','I*sqrt(3)','I*sqrt(12)'],'2**(1/3)*I*sqrt(3)',6),
    (['3**(1/3)','I*sqrt(3)','I*sqrt(12)'],'3**(5/6)*I',6),
    (['sqrt(2)','sqrt(3)','sqrt(5)'],'sqrt(30)',8),
    (['sqrt(6)','sqrt(15)','sqrt(7)','sqrt(70)'],'sqrt(10)',8),
])
def test_canonical_products_dependencies_and_roundtrip(items,member,expected_degree):
    x=sm.Symbol('x');extensions=[sm.sympify(item) for item in items]
    a=sm.sympify(member)
    p=sm.Poly(x+a,x,extension=extensions)
    assert sm.minpoly(sm.sympify(p.domain[3:-1]),x,polys=True).degree()==expected_degree
    assert p==sm.Poly(p.as_expr(),x,domain=p.domain)
    assert sm.Poly(x*x-a*a,x,extension=extensions).rem(sm.Poly(x-a,x,extension=extensions)).is_zero


@pytest.mark.parametrize('field',[
    '2**(1/3)+sqrt(2)+sqrt(3)','sqrt(2)+sqrt(3)+sqrt(5)+sqrt(7)',
    '2**(1/3)+3**(1/3)','sqrt(1+2**(1/3))','(-2)**(1/6)',
])
def test_outside_families_explicitly_unsupported(field):
    x=sm.Symbol('x');theta=sm.sympify(field)
    with pytest.raises((ValueError,TypeError,NotImplementedError)):
        sm.Poly(x+theta,x,extension=theta)


@pytest.mark.parametrize('field',FIELDS)
def test_norm_is_exact_multiplication_determinant_including_denominator_scaling(field):
    x=sm.Symbol('x');sx=sp.Symbol('x');t=sp.Symbol('t')
    theta=sm.sympify(field);eta=sp.sympify(field)
    minimal=sp.Poly(sp.minimal_polynomial(eta,t),t).monic()
    degree=minimal.degree()
    companion=sp.zeros(degree)
    for j in range(degree-1):companion[j+1,j]=1
    for i in range(degree):companion[i,degree-1]=-minimal.nth(i)
    determinant=((sp.Rational(2,3)*sx)*sp.eye(degree)-companion/5).det(method='domain-ge')
    p=sm.Poly(sm.Rational(2,3)*x-theta/5,x,extension=theta)
    same(p.norm().as_expr(),determinant,x)


@pytest.mark.parametrize('coefficient',['2**(1/3)*I*sqrt(3)','3**(5/6)*I','sqrt(6)'])
def test_inference_recombines_numeric_coefficients_without_splitting_the_field(coefficient):
    x=sm.Symbol('x',positive=True)
    a=sm.sympify(coefficient)
    p=sm.Poly(a*x+x*x,x,extension=True)
    assert p.gens==(x,)
    expected_degree=2 if coefficient=='sqrt(6)' else 6
    assert sm.minpoly(sm.sympify(p.domain[3:-1]),sm.Symbol('t'),polys=True).degree()==expected_degree
    squared=p*p
    restored=sm.Poly(squared.as_expr(),x,extension=True)
    assert restored.gens==(x,)
    same(restored.as_expr(),squared.as_expr(),x)
    assert sm.Poly(restored.as_expr()-squared.as_expr(),x,domain=p.domain).is_zero
    assert p.rem(sm.Poly(x,x,extension=True)).is_zero


def test_multivariate_inference_preserves_generator_identity_and_assumptions():
    x=sm.Symbol('x',positive=True)
    y=sm.Symbol('y',real=True)
    p=sm.Poly(sm.sqrt(2)*x+sm.sqrt(3)*y,x,y,extension=True)
    assert p.gens==(x,y)
    assert p.total_degree()==1
    for exponent,expected in [((1,0),sm.sqrt(2)),((0,1),sm.sqrt(3))]:
        actual=p.as_dict()[exponent]
        same(actual,expected)
        assert sm.Poly(actual-expected,sm.Symbol('t'),domain=p.domain).is_zero
