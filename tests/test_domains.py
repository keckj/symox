"""Native domains preserve exact coefficients, embeddings and symbol identity."""
import pickle
import subprocess
import sys

import pytest
import symox as sm


def domains():
    a=sm.Symbol('a', positive=True)
    b=sm.Symbol('b', real=True)
    return [sm.ZZ, sm.QQ, sm.GF(5),sm.GF(7),sm.QQ.algebraic_field(sm.sqrt(2)),sm.QQ.algebraic_field(2**sm.Rational(1,3)+sm.I),sm.QQ.frac_field(a,b),sm.QQ.poly_ring(a,b),sm.ZZ.poly_ring(a)]


@pytest.mark.parametrize('index',range(9))
def test_descriptor_roundtrip_and_closed_arithmetic(index):
    k=domains()[index]
    assert isinstance(k,sm.Domain) and k.is_Ring and k.is_Exact
    assert k != str(k) and str(k)==repr(k)
    assert k==pickle.loads(pickle.dumps(k))
    assert hash(k)==hash(pickle.loads(pickle.dumps(k)))
    a,b=k(2),k(3)
    assert isinstance(a,k.dtype) and k.of_type(a)
    assert a.parent()==k
    assert (a+b).parent()==k and (a*b).parent()==k
    assert k.to_sympy(a+b)==k.to_sympy(k(5))
    assert a-a==k.zero and a**0==k.one
    assert pickle.loads(pickle.dumps(a))==a
    if k.is_Field:
        assert a*(a**-1)==k.one
        assert k.exquo(a,a)==k.one
    assert k.characteristic() == (5 if index==2 else 7 if index==3 else 0)


def test_integer_rational_conversion_and_numeric_hashes():
    assert sm.ZZ.get_field()==sm.QQ and sm.QQ.get_ring()==sm.ZZ
    assert sm.QQ(1,2)==sm.Rational(1,2)
    assert sm.ZZ(2)==sm.QQ(2)==2
    assert hash(sm.QQ(1,2))==hash(sm.Rational(1,2))
    assert {sm.ZZ(2),sm.QQ(2),2}=={2}
    assert sm.QQ.convert(sm.ZZ(2),sm.ZZ)==sm.QQ(2)
    assert int(sm.QQ.numer(sm.QQ(2,3)))==2
    assert int(sm.QQ.denom(sm.QQ(2,3)))==3
    assert sm.ZZ.div(-7,3)==(sm.ZZ(-3),sm.ZZ(2))
    assert sm.QQ.div(2,3)==(sm.QQ(2,3),sm.QQ.zero)
    for k,value in [(sm.ZZ,sm.Rational(1,2)),(sm.QQ,sm.sqrt(2)),(sm.GF(5),sm.Rational(1,2))]:
        with pytest.raises(sm.CoercionFailed):k.convert(value)
    with pytest.raises(sm.CoercionFailed):sm.ZZ(2)**-1
    with pytest.raises(ZeroDivisionError):sm.QQ(1)/sm.QQ.zero


@pytest.mark.parametrize('prime',[2,3,5,101,65537])
def test_finite_fields_exact_inverse_and_conversion(prime):
    k=sm.GF(prime)
    assert k.is_FF and k.is_FiniteField and k.is_Field
    assert int(k(-1))==prime-1
    assert k(prime)==k.zero
    assert k(2 if prime>2 else 1)**(prime-1)==k.one
    assert sm.ZZ.convert(k(-1))==sm.ZZ(1 if prime==2 else -1)
    with pytest.raises(sm.DomainError):k.get_ring()
    with pytest.raises(sm.CoercionFailed):k(1)+sm.GF(103)(1)


@pytest.mark.parametrize('modulus',[0,-4])
def test_invalid_characteristics(modulus):
    with pytest.raises(ValueError,match='modulus must be a positive integer'):sm.GF(modulus)


@pytest.mark.parametrize('modulus',[1,4,6,4294967296])
def test_composite_and_large_characteristics_are_rings_modulo_n(modulus):
    # The reference GF(n) is Z/nZ for every positive n; it is a field only
    # for prime n.
    k=sm.GF(modulus)
    assert k.is_Field is False and k.mod==modulus and str(k)=='GF(%d)'%modulus


def test_function_fields_assumptions_order_membership_and_fraction_parts():
    a,b=sm.symbols('a b',real=True)
    k=sm.QQ.frac_field(a,b,order='grevlex')
    ring=k.get_ring()
    assert k.symbols==(a,b) and ring.symbols==(a,b)
    assert k.dom==sm.QQ and ring.dom==sm.QQ
    assert k.order=='grevlex' and ring.is_PolynomialRing
    assert k!=sm.QQ.frac_field(sm.Symbol('a'),b,order='grevlex')
    assert k!=sm.QQ.frac_field(a,b,order='lex')
    u=k.convert((a+1)/(b-1))
    assert u*u**-1==k.one
    assert k.numer(u).parent()==ring and k.denom(u).parent()==ring
    assert k.convert(k.numer(u))/k.convert(k.denom(u))==u
    assert a in ring and 1/a not in ring
    with pytest.raises(sm.CoercionFailed):ring.convert(1/a)
    assert sm.QQ[a,b]==sm.QQ.poly_ring(a,b)


@pytest.mark.parametrize('generator',[sm.sqrt(2),sm.I,2**sm.Rational(1,3),sm.sqrt(2)+sm.sqrt(3),2**sm.Rational(1,3)+sm.I,sm.sqrt(2)+sm.sqrt(3)+sm.sqrt(5)])
def test_algebraic_field_embedding_and_membership(generator):
    k=sm.QQ.algebraic_field(generator,alias='alpha')
    assert str(k)=='QQ<alpha>' and k.is_AlgebraicField and k.dom==sm.QQ
    assert k.ext.alias==sm.Symbol('alpha')
    a=k.from_sympy(generator)
    assert a==k.unit and k.to_sympy(a)==generator
    assert a*a**-1==k.one
    assert all(c.parent()==sm.QQ for c in a.to_list())
    assert k != sm.QQ.algebraic_field(generator,alias='beta')
    with pytest.raises(sm.DomainError):k.get_ring()


@pytest.mark.parametrize('order',['lex','grlex','grevlex'])
@pytest.mark.parametrize('kind',['QQ','GF','AN','RF','ring'])
@pytest.mark.parametrize('method',['buchberger','f5b'])
def test_poly_groebner_reduced_typed_domains(kind,order,method):
    x,y=sm.symbols('x y')
    a=sm.Symbol('a',positive=True)
    unused=sm.Symbol('unused',negative=True)
    k={'QQ':sm.QQ,'GF':sm.GF(5),'AN':sm.QQ.algebraic_field(sm.sqrt(2),alias='alpha'),'RF':sm.QQ.frac_field(a,unused,order='grlex'),'ring':sm.QQ.poly_ring(a,unused,order='grlex')}[kind]
    c=sm.sqrt(2) if kind=='AN' else a if kind in ['RF','ring'] else sm.Integer(2)
    f=sm.Poly(x*x+c*y,x,y,domain=k)
    assert f.domain==k and f.get_domain()==k
    assert pickle.loads(pickle.dumps(f)).domain==k
    g=sm.groebner([f,sm.Poly(y*y-1,x,y,domain=k)],x,y,domain=k,order=order,method=method,polys=True)
    assert g.domain==k and all(p.domain==k for p in g)
    h=pickle.loads(pickle.dumps(g))
    assert h.domain==k and h.exprs==g.exprs
    q,r=sm.reduced(f,[f],x,y,domain=k,order=order,polys=True)
    assert r.is_zero and all(p.domain==k for p in q+[r])
    if kind!='ring':
        assert g.fglm('lex').domain==k


def test_poly_domain_changes_are_checked_and_preserve_selected_field():
    x=sm.Symbol('x')
    p=sm.Poly(x/2+1,x,domain=sm.QQ)
    assert p.to_field().domain==sm.QQ
    with pytest.raises(sm.CoercionFailed):p.to_ring()
    assert sm.Poly(x+1,x).to_field().domain==sm.QQ
    k=sm.QQ.algebraic_field(sm.sqrt(2),alias='a')
    l=sm.QQ.algebraic_field(sm.sqrt(2),alias='b')
    p=sm.Poly(x+sm.sqrt(2),x,domain=k)
    q=sm.Poly(x-sm.sqrt(2),x,domain=l)
    assert (p*q).as_expr()==x*x-2
    assert (p+q).domain==k
    assert p.set_domain(l).domain==l
    coeffs=sm.Poly([k.convert(sm.sqrt(2)),k.one],x,domain=k)
    assert coeffs.as_expr()==sm.sqrt(2)*x+1


@pytest.mark.parametrize('kind',['QQ','ZZ','GF','AN','RF','ring'])
def test_differential_coefficients_and_domain_operations(kind):
    sy=pytest.importorskip('sympy')
    x=sm.Symbol('x');ox=sy.Symbol('x')
    k={'QQ':sm.QQ,'ZZ':sm.ZZ,'GF':sm.GF(5),'AN':sm.QQ.algebraic_field(sm.sqrt(2)),'RF':sm.QQ.frac_field(x),'ring':sm.QQ.poly_ring(x)}[kind]
    ok={'QQ':sy.QQ,'ZZ':sy.ZZ,'GF':sy.GF(5),'AN':sy.QQ.algebraic_field(sy.sqrt(2)),'RF':sy.QQ.frac_field(ox),'ring':sy.QQ.poly_ring(ox)}[kind]
    for name in ['is_QQ','is_ZZ','is_Field','is_Ring','is_Exact','is_AlgebraicField','is_FractionField','is_PolynomialRing']:
        assert getattr(k,name)==getattr(ok,name)
    for a,b in [(2,3),(0,3),(-6,4)]:
        for op in ['gcd','lcm']:
            actual=k.to_sympy(getattr(k,op)(k(a),k(b)))
            expected=ok.to_sympy(getattr(ok,op)(ok(a),ok(b)))
            assert sy.simplify(sy.sympify(str(actual))-expected)==0


def test_import_paths_and_no_sympy_runtime():
    code=r'''
import sys,pickle
class NoSymPy:
    def find_spec(self,name,*args):
        if name=='sympy' or name.startswith('sympy.'):raise RuntimeError('runtime SymPy')
sys.meta_path.insert(0,NoSymPy())
import symox as s
from symox.polys.domains import QQ,ZZ,GF,Domain
from symox.polys.polyerrors import DomainError,CoercionFailed
x=s.Symbol('x');a=s.Symbol('a',positive=True);b=s.Symbol('b',negative=True)
for k,c in [(QQ,s.Integer(2)),(GF(7),s.Integer(2)),(QQ.algebraic_field(s.sqrt(2),alias='alpha'),s.sqrt(2)),(QQ.frac_field(a,b),1/a)]:
    v=k.convert(c);assert v*v**-1==k.one
    p=s.Poly(x*x+c*x+1,x,domain=k)
    q=pickle.loads(pickle.dumps(p));assert q==p and q.domain==k
    g=s.groebner([p],x,domain=k,method='f5b')
    assert pickle.loads(pickle.dumps(g)).domain==k
assert 'sympy' not in sys.modules
'''
    subprocess.run([sys.executable,'-c',code],check=True)


@pytest.mark.parametrize('kind',['RF','ring','AN'])
@pytest.mark.parametrize('method',['buchberger','f5b'])
def test_inferred_domains_keep_unused_parameters_aliases_and_orders(kind,method):
    x=sm.Symbol('x');a=sm.Symbol('a',positive=True);b=sm.Symbol('b',negative=True)
    k={'RF':sm.QQ.frac_field(a,b,order='grevlex'),'ring':sm.QQ.poly_ring(a,b,order='grlex'),'AN':sm.QQ.algebraic_field(sm.sqrt(2),alias='alpha')}[kind]
    p=sm.Poly(x+1,x,domain=k)
    for q in [p+p,p*p,p.gcd(p),p.diff()]:
        assert q.domain==k
    g=sm.groebner([p,p*p],x,method=method)
    assert g.domain==k and all(q.domain==k for q in g)
    q,r=sm.reduced(p,[p],x,polys=True)
    assert all(v.domain==k for v in q+[r])
    assert pickle.loads(pickle.dumps(g)).domain==k


@pytest.mark.parametrize('kind',['ZZ','QQ','GF'])
def test_numeric_domain_elements_interoperate_with_expressions(kind):
    x=sm.Symbol('x');k={'ZZ':sm.ZZ,'QQ':sm.QQ,'GF':sm.GF(5)}[kind]
    c=k(2)
    assert c+x==x+2 and x+c==x+2
    assert c-x==2-x and x-c==x-2
    assert c*x==2*x and x*c==2*x
    assert c/x==2/x and x/c==x/2
    assert x**c==x*x
    if kind!='GF':
        assert sm.Integer(2)==c and c==sm.Integer(2)
        assert {sm.Integer(2),c,2}=={2}
    else:
        # The reference ModularInteger compares integers modulo p.
        assert c==2 and 2==c and c==7 and not c!=2 and c==sm.Integer(2)


@pytest.mark.parametrize('ground',['QQ','ZZ'])
@pytest.mark.parametrize('ring',[False,True])
def test_differential_nonconstant_gcd_lcm(ground,ring):
    sy=pytest.importorskip('sympy')
    x=sm.Symbol('x');ox=sy.Symbol('x')
    native_ground=getattr(sm,ground);oracle_ground=getattr(sy,ground)
    k=(native_ground.poly_ring if ring else native_ground.frac_field)(x)
    ok=(oracle_ground.poly_ring if ring else oracle_ground.frac_field)(ox)
    pairs=[('2*x','4*x'),('2*x+2','4*x+4'),('0','-3*x-3'),('2*x**2+4*x','3*x**2+6*x')]
    if not ring:pairs += [('2*x/(2*x+1)','3*x/(3*x+1)'),('2*x/3','4*x/5')]
    for a,b in pairs:
        for op in ['gcd','lcm']:
            actual=k.to_sympy(getattr(k,op)(k(a),k(b)))
            expected=ok.to_sympy(getattr(ok,op)(ok.from_sympy(sy.sympify(a)),ok.from_sympy(sy.sympify(b))))
            assert sy.simplify(sy.sympify(str(actual))-expected)==0,(k,a,b,op,actual,expected)


def test_bezout_modular_inverse_and_polynomial_ring_euclidean_certificates():
    k=sm.ZZ
    for a,b in [(30,21),(-30,21),(0,7),(7,0),(0,0)]:
        u,v,g=k.gcdex(k(a),k(b))
        assert u*k(a)+v*k(b)==g==k.gcd(k(a),k(b))
    assert int(k.invert(k(3),k(7)))==5
    with pytest.raises(sm.NotInvertible):k.invert(k(2),k(4))
    x=sm.Symbol('x');r=sm.QQ.poly_ring(x)
    a=r(x+1);b=r(x*x+1)
    u,v,g=r.gcdex(a,b)
    assert u*a+v*b==g==r.one
    inverse=r.invert(a,b)
    assert sm.Poly(r.to_sympy(a*inverse-r.one),x,domain=sm.QQ).rem(sm.Poly(x*x+1,x,domain=sm.QQ)).is_zero


def test_inference_from_plain_poly_matches_sympy_reinterpretation():
    x=sm.Symbol('x')
    for k in [sm.QQ,sm.GF(5)]:
        p=sm.Poly(2*x+1,x,domain=k)
        g=sm.groebner([p],x)
        assert g.domain==sm.ZZ and all(q.domain==sm.ZZ for q in g)
        assert g.contains(p)


@pytest.mark.parametrize('index',range(9))
def test_numpy_object_arrays_keep_exact_domain_elements(index):
    np=pytest.importorskip('numpy')
    k=domains()[index]
    a=np.array([k.one,k(2)],dtype=object)
    b=np.array([k(2),k.one],dtype=object)
    result=np.dot(a,b)
    assert result.parent()==k and result==k(4)
    assert all(v.parent()==k for v in a+b)


def test_algebraic_element_aliases_preserve_embedding_arithmetic_and_hash():
    k=sm.QQ.algebraic_field(sm.sqrt(2),alias='alpha')
    l=sm.QQ.algebraic_field(sm.sqrt(2),alias='beta')
    a=k.convert(sm.sqrt(2));b=l.convert(sm.sqrt(2))
    assert k!=l and a==b and b==a
    assert hash(a)==hash(b) and len({a,b})==1
    assert (a+b).parent()==k and (b+a).parent()==l
    assert k.to_sympy(a+b)==2*sm.sqrt(2)
    assert (a*b).parent()==k and a/b==k.one
    j=sm.QQ.algebraic_field(sm.sqrt(8))
    c=j.convert(sm.sqrt(2))
    assert c!=a
    with pytest.raises(sm.CoercionFailed):a+c
    assert k.convert(c)==a
    assert pickle.loads(pickle.dumps(a))==b


def test_legacy_unbound_parameter_names_do_not_erase_assumptions():
    x=sm.Symbol('x');a=sm.Symbol('a',nonzero=True);ordinary=sm.Symbol('a')
    first=sm.Poly(a*x-1,x,domain='QQ(a,b)')
    second=sm.Poly(x*x+1,x,domain='QQ(a,b)')
    assert a in first.domain.symbols and ordinary in second.domain.symbols
    g=sm.groebner([first,second],x)
    assert len(g.domain.symbols)==3
    assert a in g.domain.symbols and ordinary in g.domain.symbols
    assert len([s for s in g.domain.symbols if s.name=='a'])==2
    assert sm.groebner([first,second],x).domain==g.domain
    assert pickle.loads(pickle.dumps(g)).domain==g.domain


@pytest.mark.parametrize('kind',['GF','RF'])
def test_unconnected_expression_operations_reject_typed_domains(kind):
    x=sm.Symbol('x');a=sm.Symbol('a')
    k=sm.GF(5) if kind=='GF' else sm.QQ.frac_field(a)
    if kind=='GF':
        # The reference treats domain=GF(p) as modulus=p.
        assert sm.factor(x*x+1,x,domain=k)==(x-2)*(x+2)
        assert sm.invert(x,x*x+1,x,domain=k)==-x
        return
    with pytest.raises(NotImplementedError):sm.factor(x*x+1,x,domain=k)
    with pytest.raises(NotImplementedError):sm.invert(x,x*x+1,x,domain=k)


@pytest.mark.parametrize('order',['lex','grlex','grevlex'])
def test_declared_parameter_order_survives_polynomial_and_basis_operations(order):
    a=sm.Symbol('a',positive=True);b=sm.Symbol('b',negative=True)
    x,y=sm.symbols('x y')
    k=sm.QQ.frac_field(b,a,order=order)
    p=sm.Poly(x*x+a*y+1,x,y,domain=k)
    for value in [p+p,p*p,p.diff(x),p.integrate(x),p.gcd(p*p),p.clear_denoms()[1]]:
        assert value.domain==k
        assert value.domain.symbols==(b,a)
    g=sm.groebner([p,sm.Poly(y*y-1,x,y,domain=k)],x,y,domain=k,order='grevlex')
    for value in [g,g.fglm('lex'),pickle.loads(pickle.dumps(g))]:
        assert value.domain==k and value.domain.symbols==(b,a)


def test_dummy_parameter_identity_survives_typed_domain_and_basis_pickles():
    x=sm.Symbol('x');a=sm.Dummy('a',positive=True);b=sm.Dummy('a',positive=True)
    assert a!=b
    k=sm.QQ.frac_field(a,b,order='grlex')
    p=sm.Poly(a*x+1,x,domain=k)
    assert p.domain.symbols==(a,b)
    q=pickle.loads(pickle.dumps(p))
    assert q.domain==k and q.domain.symbols==(a,b) and q==p
    g=sm.groebner([p],x,domain=k,method='f5b')
    restored=pickle.loads(pickle.dumps(g))
    assert restored.domain==k and restored.domain.symbols==(a,b)
    assert restored.exprs==g.exprs
