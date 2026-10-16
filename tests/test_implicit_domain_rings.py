"""Implicit symbolic coefficients select an actual exact ring backend."""
import pickle
import subprocess
import sys
import pytest
import symox as sm


@pytest.mark.parametrize('rational',[False,True])
@pytest.mark.parametrize('dummy',[False,True])
def test_implicit_parameter_ring_coefficients_and_pickle(rational,dummy):
    sy=pytest.importorskip('sympy')
    x,y=sm.symbols('x y')
    factory=sm.Dummy if dummy else sm.Symbol
    a=factory('a',positive=True);b=factory('b',nonzero=True)
    c=sm.Rational(2,3) if rational else sm.Integer(2)
    expr=c*a*x*x+(a+b)*y+b
    p=sm.Poly(expr,x,y)
    oracle=sy.Poly(sm.to_sympy(expr),sm.to_sympy(x),sm.to_sympy(y))
    assert p.gens==(x,y) and p.domain.is_PolynomialRing
    assert set(p.domain.symbols)=={a,b}
    assert p.domain.dom==(sm.QQ if rational else sm.ZZ)
    assert str(p.domain)==str(oracle.domain)
    assert sm.to_sympy(p.as_expr())==oracle.as_expr()
    assert [sm.to_sympy(c) for c in p.coeffs()]==oracle.coeffs()
    assert p.coeff_monomial(x*x)==c*a
    assert all(p.domain.convert(c).parent()==p.domain for c in p.as_dict().values())
    q=pickle.loads(pickle.dumps(p))
    assert q.domain==p.domain and q.gens==p.gens and q==p
    g=sm.groebner([p],x,y)
    assert g.domain==p.domain and all(q.domain==p.domain for q in g)
    assert g.contains(p)
    primitive=p.integrate()
    assert primitive.domain==p.domain.get_field()
    assert primitive.diff(x)==p.set_domain(primitive.domain)


@pytest.mark.parametrize('rational',[False,True])
def test_list_coefficients_and_monomial_access_use_ring_coefficients(rational):
    x=sm.Symbol('x');a=sm.Symbol('a',positive=True);b=sm.Symbol('b',negative=True)
    c=sm.Rational(2,3) if rational else sm.Integer(2)
    p=sm.Poly([c*a,a+b,b],x)
    assert p.gens==(x,)
    assert p.all_coeffs()==[c*a,a+b,b]
    assert p.coeff_monomial(x*x)==c*a
    assert p.domain==((sm.QQ if rational else sm.ZZ).poly_ring(a,b))
    q,r=p.div(sm.Poly(x+1,x))
    assert q.domain==p.domain and r.domain==p.domain
    assert sm.Poly(q.as_expr()*(x+1)+r.as_expr()-p.as_expr(),x,domain=p.domain).is_zero
    assert all(p.domain.convert(c).parent()==p.domain for c in q.coeffs()+r.coeffs())


@pytest.mark.parametrize('rational',[False,True])
@pytest.mark.parametrize('pair',range(6))
def test_migration_gcd_keeps_coefficient_ring_content(rational,pair):
    sy=pytest.importorskip('sympy')
    x=sm.Symbol('x');a=sm.Symbol('a',positive=True)
    examples=[(a*x,a*x*x),(2*a*x,4*a*x*x),((a+1)*x,(a+1)*x*x),(2*(a+1)*x,4*(a+1)*x*x),(2*a*x+2*a,4*a*x*x+8*a*x+4*a),(a*x+1,2*a*x+2)]
    lhs,rhs=examples[pair]
    if rational:lhs=lhs/3;rhs=rhs/5
    p=sm.Poly(lhs,x);q=sm.Poly(rhs,x)
    oracle_p=sy.Poly(sm.to_sympy(lhs),sm.to_sympy(x));oracle_q=sy.Poly(sm.to_sympy(rhs),sm.to_sympy(x))
    g=p.gcd(q);expected=oracle_p.gcd(oracle_q)
    assert g.domain==p.domain
    assert sy.expand(sm.to_sympy(g.as_expr())-expected.as_expr())==0
    assert p.rem(g).is_zero and q.rem(g).is_zero
    assert all(g.domain.convert(c).parent()==g.domain for c in g.coeffs())


def test_legacy_implicit_ring_pickle_is_upgraded_without_string_identity():
    x=sm.Symbol('x');a=sm.Symbol('a',positive=True)
    p=sm._core._rebuild_poly(a*x+1,[x],None,'ZZ[a]',None)
    assert p.domain==sm.ZZ.poly_ring(a)
    assert p.integrate().domain==sm.ZZ.frac_field(a)
    assert p.coeff_monomial(x)==a


@pytest.mark.parametrize('domain',['ZZ','QQ',sm.ZZ,sm.QQ])
def test_explicit_ground_domain_cannot_hide_symbolic_coefficients(domain):
    x,a=sm.symbols('x a')
    with pytest.raises(sm.CoercionFailed):sm.Poly(a*x+1,x,domain=domain)


@pytest.mark.parametrize('coefficient',[sm.sin(sm.Symbol('a')),sm.sqrt(2),sm.I])
def test_legacy_expression_coefficients_keep_calculations_and_pickle(coefficient):
    x=sm.Symbol('x')
    p=sm.Poly(coefficient*x+1,x)
    assert p.as_expr()==coefficient*x+1
    assert sm.expand(p.integrate().diff().as_expr()-p.as_expr())==0
    assert pickle.loads(pickle.dumps(p))==p
    with pytest.raises(NotImplementedError):p.domain


def test_implicit_domain_ring_without_runtime_sympy():
    subprocess.run([sys.executable,'-c',r'''
import sys,pickle
class NoSymPy:
 def find_spec(self,name,*args):
  if name=='sympy' or name.startswith('sympy.'):raise RuntimeError('runtime SymPy')
sys.meta_path.insert(0,NoSymPy())
import symox as s
x=s.Symbol('x');a=s.Symbol('a',positive=True)
p=s.Poly(a*x,x);q=s.Poly(a*x*x,x)
assert p.domain==s.ZZ.poly_ring(a)
assert p.gcd(q).as_expr()==a*x
assert p.integrate().domain==s.ZZ.frac_field(a)
assert pickle.loads(pickle.dumps(p)).domain==p.domain
assert p.domain.convert(p.coeff_monomial(x)).parent()==p.domain
assert 'sympy' not in sys.modules
'''],check=True)


def test_inexact_symbolic_coefficients_keep_legacy_behavior():
    x,a=sm.symbols('x a')
    p=sm.Poly(sm.Float(0.5)*a*x+1,x)
    assert sm.simplify(p.integrate().diff().as_expr()-p.as_expr())==0
    assert pickle.loads(pickle.dumps(p))==p
    with pytest.raises(NotImplementedError):p.domain
