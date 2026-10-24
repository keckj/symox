"""Exact ground roots and included units retain typed domains and certificates."""
import pickle
import subprocess
import sys

import pytest
import symox as sm

from test_global_poly_division import selected


def reconstruct(p,factors):
    result=sm.Poly(1,*p.gens,domain=p.domain)
    for factor,multiplicity in factors:
        assert isinstance(factor,sm.Poly) and factor.gens==p.gens and factor.domain==p.domain
        assert isinstance(multiplicity,int) and multiplicity>=1
        assert all(coefficient in p.domain for coefficient in factor.coeffs())
        result *= factor**multiplicity
        restored=pickle.loads(pickle.dumps(factor))
        assert restored.gens==factor.gens and restored.domain==p.domain and restored.as_expr()==factor.as_expr()
    assert sm.Poly(sm.expand(result.as_expr()-p.as_expr()),*p.gens,domain=p.domain).is_zero


def compare_roots(actual,expected):
    assert isinstance(actual,dict)
    assert {sm.to_sympy(root):multiplicity for root,multiplicity in actual.items()}==expected
    assert all(isinstance(root,sm.Expr) and isinstance(multiplicity,int) and multiplicity>0 for root,multiplicity in actual.items())


@pytest.mark.parametrize('kind',['ZZ','QQ','GF','AN','QQRF','ZZRF','QQring','ZZring'])
@pytest.mark.parametrize('order',['lex','grlex','grevlex'])
@pytest.mark.parametrize('case',range(4))
def test_typed_factor_views_ground_oracle_metadata_and_exact_reconstruction(kind,order,case):
    import sympy as sy
    x=sm.Symbol('x')
    domain,odomain,a,b=selected(kind,order)
    c=a/(b+1) if kind.endswith('RF') else a if kind.endswith('ring') else sm.sqrt(2) if kind=='AN' else sm.Integer(2)
    expression=[(c*x-1)**2*(x-2)*(x*x+1),c*(x+1)**3,sm.Integer(0),7*c][case]
    p=sm.Poly(expression,x,domain=domain)
    op=sy.Poly(sm.to_sympy(expression),sm.to_sympy(x),domain=odomain)
    roots=p.ground_roots()
    compare_roots(roots,op.ground_roots())
    compare_roots(sm.ground_roots(p),op.ground_roots())
    compare_roots(sm.ground_roots(expression,x,domain=domain),op.ground_roots())
    associated=domain.get_field() if not domain.is_Field else domain
    for root,multiplicity in roots.items():
        assert root in associated
        residue=sm.Poly(p.as_expr(),x,domain=associated).eval(root)
        assert associated.convert(residue)==associated.zero
    for factors in [p.factor_list_include(),p.sqf_list_include(),p.sqf_list_include(all=True)]:
        reconstruct(p,factors)
        for factor,_ in factors:
            if kind.endswith('RF') or kind.endswith('ring'):
                assert factor.domain.symbols==(a,b) and factor.domain.order==order
    if p.is_zero or p.total_degree()==0:
        assert p.factor_list_include()==[(p,1)]
        assert p.sqf_list_include()==[(p,1)]
        assert roots=={}


@pytest.mark.parametrize('coefficient',[2,-2,sm.Rational(1,3),sm.Rational(-7,5)])
@pytest.mark.parametrize('power',[2,3,4])
def test_repeated_factor_units_are_never_raised_to_its_multiplicity(coefficient,power):
    import sympy as sy
    x=sm.Symbol('x')
    p=sm.Poly(coefficient*(x+1)**power,x,domain=sm.QQ)
    included=p.factor_list_include()
    reconstruct(p,included)
    assert any(int(factor.degree())==0 and multiplicity==1 for factor,multiplicity in included)
    if power==2 and coefficient==2:
        # SymPy 1.14 absorbs 2 into x+1 while retaining exponent 2. Its
        # returned product is 4*(x+1)^2, so the native form keeps (2,1).
        assert included==[(sm.Poly(2,x,domain=sm.QQ),1),(sm.Poly(x+1,x,domain=sm.QQ),2)]
        reference=sy.Poly(2*(sm.to_sympy(x)+1)**2,sm.to_sympy(x),domain=sy.QQ).factor_list_include()
        assert sy.expand(sy.prod(f.as_expr()**k for f,k in reference)-sm.to_sympy(p.as_expr())) != 0


@pytest.mark.parametrize('domain',[sm.ZZ,sm.QQ,sm.GF(5)])
def test_unit_absorption_when_multiplicity_one_and_squarefree_placeholders(domain):
    x=sm.Symbol('x')
    p=sm.Poly(2*(x-1)*(x+1),x,domain=domain)
    included=p.factor_list_include()
    reconstruct(p,included)
    assert len(included)==2 and all(m==1 for _,m in included)
    p=sm.Poly(2*(x-1)**3,x,domain=domain)
    included=p.sqf_list_include(all=True)
    reconstruct(p,included)
    assert [m for _,m in included]==[1,2,3]
    assert included[1][0]==sm.Poly(1,x,domain=domain)


@pytest.mark.parametrize('generator',['sqrt2','I','cubic_quadratic','three_quadratics'])
def test_algebraic_embeddings_and_aliases_remain_exact(generator):
    import sympy as sy
    x=sm.Dummy('x',real=True)
    c={'sqrt2':sm.sqrt(2),'I':sm.I,'cubic_quadratic':2**sm.Rational(1,3)+sm.I,'three_quadratics':sm.sqrt(2)+sm.sqrt(3)+sm.sqrt(5)}[generator]
    domain=sm.QQ.algebraic_field(c,alias='alpha')
    p=sm.Poly(3*(x-c)**2*(x+1),x,domain=domain)
    roots=p.ground_roots()
    assert roots=={c:2,sm.Integer(-1):1}
    for factors in [p.factor_list_include(),p.sqf_list_include(all=True)]:
        reconstruct(p,factors)
        assert all(f.domain.ext.alias==sm.Symbol('alpha') and f.gens==(x,) for f,_ in factors)
    assert all(v.gens==(x,) for v,_ in pickle.loads(pickle.dumps(p.factor_list_include())))


def test_associated_fraction_fields_allow_roots_outside_the_original_ring():
    x,a=sm.symbols('x a')
    for domain,expression,expected in [(sm.ZZ,(2*x-1)**3,{sm.Rational(1,2):3}),(sm.ZZ.poly_ring(a),(a*x-1)**2,{1/a:2}),(sm.QQ.poly_ring(a),(a*x-1)**2,{1/a:2})]:
        p=sm.Poly(expression,x,domain=domain)
        assert p.ground_roots()==expected
        assert all(root in domain.get_field() for root in expected)
        reconstruct(p,p.factor_list_include())
        reconstruct(p,p.sqf_list_include())


def test_finite_inseparability_and_large_prime_ground_roots():
    x=sm.Symbol('x')
    for prime,expression,expected in [(2,(x+1)**8,{sm.Integer(1):8}),(5,x**5-x,{sm.Integer(v):1 for v in [-2,-1,0,1,2]}),(1000000007,(x-3)**2*(x+7),{sm.Integer(3):2,sm.Integer(-7):1})]:
        p=sm.Poly(expression,x,modulus=prime)
        assert p.ground_roots()==expected
        reconstruct(p,p.factor_list_include())
        reconstruct(p,p.sqf_list_include(all=True))


@pytest.mark.parametrize('name',['factor_list_include','sqf_list_include'])
def test_multivariate_rational_and_function_coefficient_factor_metadata(name):
    x,y,a=sm.symbols('x y a')
    for domain in [sm.ZZ,sm.QQ,sm.ZZ.poly_ring(a),sm.QQ.frac_field(a)]:
        coefficient=a if domain not in [sm.ZZ,sm.QQ] else sm.Integer(2)
        p=sm.Poly(coefficient*(x+y)**2*(x-y),x,y,domain=domain)
        reconstruct(p,getattr(p,name)())
        # the reference MultivariatePolynomialError of Poly.ground_roots
        with pytest.raises(sm.MultivariatePolynomialError,match='Cannot compute ground roots of'):p.ground_roots()


@pytest.mark.parametrize('option',['modulus','extension','field','gens','domain'])
def test_ground_roots_global_construction_options_are_applied(option):
    import sympy as sy
    x,y=sm.symbols('x y')
    expression=x*x-2
    native={'modulus':7} if option=='modulus' else {'extension':sm.sqrt(2)} if option=='extension' else {'field':True} if option=='field' else {'gens':(x,)} if option=='gens' else {'domain':sm.QQ}
    oracle={'modulus':7} if option=='modulus' else {'extension':sy.sqrt(2)} if option=='extension' else {'field':True} if option=='field' else {'gens':(sm.to_sympy(x),)} if option=='gens' else {'domain':sy.QQ}
    generators=() if option=='gens' else (x,)
    compare_roots(sm.ground_roots(expression,*generators,**native),sy.ground_roots(sm.to_sympy(expression),*(sm.to_sympy(g) for g in generators),**oracle))


def test_unknown_flags_and_unsupported_domains_fail_explicitly():
    x,y=sm.symbols('x y')
    for flag in ['auto','polys']:
        for value in [True,False,None,2,'bad']:
            with pytest.raises(sm.FlagError):sm.ground_roots(x*x-1,x,**{flag:value})
    with pytest.raises(sm.OptionError):sm.ground_roots(x*x-1,x,bogus=True)
    import sympy as sy
    for p in [sm.Poly(0.5*x+1,x),sm.Poly(sm.sin(y)*x+1,x)]:
        if p.domain==sm.RR:
            # The reference reads the linear factors of factor_list over RR.
            q=sy.Poly(sm.to_sympy(p.as_expr()),sm.to_sympy(x))
            assert str(p.ground_roots())==str(q.ground_roots())
        else:
            with pytest.raises(ValueError,match='exact coefficient domain'):p.ground_roots()
        # The legacy RR and ZZ[sin(y)] domains take the factorization of the
        # reference.
        q=sy.Poly(sm.to_sympy(p.as_expr()),sm.to_sympy(x))
        for method in ['factor_list_include','sqf_list_include']:
            assert str(getattr(p,method)())==str(getattr(q,method)())
    p=sm.Poly(x*y+1,x,y,modulus=5)
    with pytest.raises(ValueError,match='multivariate finite-field'):p.factor_list_include()


def test_large_shapes_and_growth_fail_with_value_error_then_restore_scope():
    x=sm.Symbol('x')
    p=sm.Poly(x**1_000_000_000+1,x)
    for method in ['ground_roots','factor_list_include','sqf_list_include']:
        with pytest.raises(ValueError,match='resource limit') as failure:getattr(p,method)()
        assert type(failure.value) is ValueError
    p=sm.Poly(sm.Integer(1<<20000)*x+1,x)
    assert p.ground_roots()=={-sm.Rational(1,1<<20000):1}
    p=sm.Poly(sm.Integer(1<<(1<<21))*x+1,x)
    with pytest.raises(ValueError,match='coefficient size'):p.ground_roots()
    assert sm.Poly(x*x-1,x).ground_roots()=={sm.Integer(-1):1,sm.Integer(1):1}
    reconstruct(sm.Poly(2*(x+1)**2,x),sm.Poly(2*(x+1)**2,x).factor_list_include())


def test_checked_factor_views_without_sympy_runtime():
    code='''
import importlib.abc,sys,pickle
class Block(importlib.abc.MetaPathFinder):
 def find_spec(self,name,path=None,target=None):
  if name=='sympy' or name.startswith('sympy.'):raise ImportError('blocked')
sys.meta_path.insert(0,Block())
import symox as s
x=s.Dummy('x',real=True);a=s.Symbol('a',positive=True)
for k in [s.ZZ,s.QQ,s.GF(5),s.QQ.algebraic_field(s.sqrt(2)),s.QQ.frac_field(a),s.ZZ.poly_ring(a)]:
 p=s.Poly(2*(x+1)**3,x,domain=k)
 assert s.ground_roots(p)=={s.Integer(-1):3}
 for factors in [p.factor_list_include(),p.sqf_list_include(all=True)]:
  result=s.Poly(1,x,domain=k)
  for q,m in factors:
   assert q.domain==k and q.gens==(x,)
   result*=q**m
   assert pickle.loads(pickle.dumps(q)).domain==k
  assert s.Poly(result.as_expr()-p.as_expr(),x,domain=k).is_zero
assert s.Poly(2*x-1,x).ground_roots()=={s.Rational(1,2):1}
assert s.Poly(a*x-1,x,domain=s.ZZ.poly_ring(a)).ground_roots()=={1/a:1}
try:s.Poly(x**1000000000+1,x).ground_roots()
except ValueError:pass
else:raise AssertionError('no budget error')
assert s.ground_roots(x*x-1,x)=={s.Integer(-1):1,s.Integer(1):1}
assert not any(n=='sympy' or n.startswith('sympy.') for n in sys.modules)
'''
    result=subprocess.run([sys.executable,'-c',code],text=True,capture_output=True)
    assert result.returncode==0,result.stderr
