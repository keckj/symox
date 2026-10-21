"""Native constructors must preserve domains and polynomial certificates."""
import pickle
import pytest
import symox as s
import sympy as sp
from symox.polys.rings import ring,xring,vring,sring
from symox.polys.fields import field,xfield,vfield,sfield


@pytest.mark.parametrize('domain',[s.ZZ,s.QQ])
@pytest.mark.parametrize('order',['lex','grlex','grevlex'])
def test_sparse_ring_arithmetic_and_division(domain,order):
    R,x,y=ring('x,y',domain,order)
    f=(x+y)*(x-y)+2*y
    g=x-y
    q,r=f.div(g)
    assert q*g+r==f
    assert f.exquo(1)==f
    assert ((x+y)*(x-y)).exquo(x-y)==x+y
    assert ((x+y)*(x-y))%(x-y)==0
    a=R.from_dict({(2,0):1,(0,2):-1,(0,1):2})
    assert a==f and a.parent()==R
    assert all(c.parent()==domain for c in a.coeffs())
    assert f.degree(x)==2 and f.degree(y)==2 and f.total_degree()==2
    assert f.diff(x)==2*x
    assert f.evaluate(x,2).as_expr()==R._native.convert(4-y.as_expr()**2+2*y.as_expr()).as_expr()
    assert R.to_domain()==domain.poly_ring(*R.symbols,order=order)
    assert pickle.loads(pickle.dumps(f))==f
    sx,sy=sp.symbols('x y')
    oracle=sp.Poly(sx*sx-sy*sy+2*sy,sx,sy,domain=str(domain))
    assert sp.expand(f.as_expr().to_sympy()-oracle.as_expr())==0


def test_gcd_multi_divisor_and_wrong_domains():
    R,x,y=ring('x,y',s.QQ)
    assert ((x+y)*(x-1)).gcd((x+y)*(y+1))==x+y
    q,r=(x*x+y*y).div([x,y])
    assert q[0]*x+q[1]*y+r==x*x+y*y
    with pytest.raises(s.ExactQuotientFailed):x.exquo(y)
    with pytest.raises(ZeroDivisionError):x.div(0)
    with pytest.raises(ValueError):R.from_dict({(-1,0):1})
    with pytest.raises(ValueError):ring('x,x',s.QQ)
    with pytest.raises(s.CoercionFailed):ring('x',s.ZZ)[0](s.Rational(1,2))


@pytest.mark.parametrize('domain',[s.ZZ,s.QQ])
def test_fraction_cancellation_and_quotient_rule(domain):
    F,x,y=field('x,y',domain)
    f=(x*x-y*y)/(x-y)
    assert f==x+y and f.denom==1
    g=(x+1)/(y+1)
    assert g.numer.as_expr()==x.as_expr()+1
    assert g.denom.as_expr()==y.as_expr()+1
    assert g.diff(x)==1/(y+1)
    assert g.diff(y)==-(x+1)/(y+1)**2
    assert g*g**-1==1
    assert F((g.numer,g.denom))==g
    assert g.evaluate(x,2).as_expr()==F(3/(y.as_expr()+1)).as_expr()
    assert pickle.loads(pickle.dumps(g))==g
    with pytest.raises(ZeroDivisionError):x/0
    with pytest.raises(ValueError):g.to_poly()


def test_constructor_shapes_inference_and_namespace_injection():
    R,gens=xring('a,b',s.ZZ)
    assert gens==R.gens
    F,gens=xfield('a,b',s.QQ)
    assert gens==F.gens
    x,y=s.symbols('x y')
    R,fs=sring([x/2+y,2*x])
    assert R.domain==s.QQ and R.symbols==(x,y) and len(fs)==2
    F,f=sfield((x*x-1)/(x-1))
    assert F.domain==s.ZZ and f.as_expr()==x+1
    namespace={'s':s,'vring':vring,'vfield':vfield}
    exec("R=vring('fresh_r,fresh_s',s.QQ)\nF=vfield('fresh_t',s.ZZ)",namespace)
    assert namespace['fresh_r']==namespace['R'].gens[0]
    assert namespace['fresh_t']==namespace['F'].gens[0]

def test_composite_domain_backend_limits_are_explicit():
    domain=s.GF(7)
    with pytest.raises(NotImplementedError,match="not implemented"):
        ring("x",domain)
    with pytest.raises(NotImplementedError,match="not implemented"):
        field("x",domain)
    # The reference QQ<sqrt(2)>[x] is a supported nested ring; its fraction
    # field has no element backend yet and says so.
    R,x=ring("x",s.QQ.algebraic_field(s.sqrt(2)))
    assert str(R)=="Polynomial ring in x over QQ<sqrt(2)> with lex order"
    assert str(x**2+x)=="x**2 + x"
    with pytest.raises(NotImplementedError,match="not implemented"):
        field("x",s.QQ.algebraic_field(s.sqrt(2)))


def test_constant_rings_and_empty_construction():
    R,f=sring(s.Rational(2,3))
    assert R.ngens==0 and R.domain==s.QQ and f.to_dict()=={():s.QQ(2,3)}
    assert f.div(R(2))==(R(s.Rational(1,3)),R.zero)
    R,values=sring([])
    assert values==[] and R.ngens==0 and R.domain==s.ZZ
    F,f=sfield(s.Rational(2,3))
    assert F.ngens==0 and f.numer==2 and f.denom==3
    with pytest.raises(NotImplementedError):ring('x',s.QQ,lambda m:m)
