"""Sparse exponent transforms preserve exact coefficient descriptors."""
import os
import pickle
import subprocess
import sys

import pytest
import symox as sm

from test_global_poly_division import selected

KINDS=['ZZ','QQ','GF','AN','QQRF','ZZRF','QQring','ZZring']


def same_poly(actual,expected,domain=None):
    import sympy as sy
    assert isinstance(actual,sm.Poly)
    assert tuple(sm.to_sympy(g) for g in actual.gens)==expected.gens
    assert str(actual.domain)==str(expected.domain)
    if domain is not None:assert actual.domain==domain
    assert sy.expand(sm.to_sympy(actual.as_expr())-expected.as_expr())==0
    restored=pickle.loads(pickle.dumps(actual))
    assert restored.gens==actual.gens and restored.domain==actual.domain
    assert restored.as_dict()==actual.as_dict()


def inflate_deflated(value,factors,original_gens):
    result=sm.Integer(0)
    for powers,coefficient in value.terms():
        term=coefficient
        for g,p,d in zip(original_gens,powers,factors):term*=g**(p*d)
        result+=term
    return sm.expand(result)


@pytest.mark.parametrize('kind',KINDS)
@pytest.mark.parametrize('order',['lex','grlex','grevlex'])
@pytest.mark.parametrize('case',range(6))
def test_structural_methods_oracle_metadata_and_reconstruction(kind,order,case):
    import sympy as sy
    k,ok,a,unused=selected(kind,order)
    if kind=='AN':
        k=sm.QQ.algebraic_field(sm.sqrt(2),alias='alpha')
        ok=sy.QQ.algebraic_field(sy.sqrt(2),alias='alpha')
    x=sm.Dummy('x',real=True);y,z=sm.symbols('y z')
    c=a/(unused+1) if kind.endswith('RF') else a if kind.endswith('ring') else sm.sqrt(2) if kind=='AN' else sm.Integer(2)
    expression=[sm.Integer(0),c,c*x**6*y**4+3*x*x*y*y,c*x*x+y*y,c*x*x+x+y,c*x**4+3*x*x][case]
    p=sm.Poly(expression,x,y,z,domain=k)
    op=sy.Poly(sm.to_sympy(expression),*(sm.to_sympy(g) for g in (x,y,z)),domain=ok)
    for method in ['deflate','terms_gcd']:
        factors,q=getattr(p,method)();ofactors,oq=getattr(op,method)()
        assert isinstance(factors,tuple) and factors==ofactors
        same_poly(q,oq,k)
        if method=='deflate':reconstructed=inflate_deflated(q,factors,p.gens)
        else:
            monomial=sm.Integer(1)
            for g,d in zip(p.gens,factors):monomial*=g**d
            reconstructed=sm.expand(monomial*q.as_expr())
        assert sm.Poly(reconstructed-p.as_expr(),*p.gens,domain=k).is_zero
    same_poly(p.exclude(),op.exclude(),k)
    for gens in [(z,x,y),(y,z,x),()]:
        actual=p.reorder(*gens);expected=op.reorder(*(sm.to_sympy(g) for g in gens))
        same_poly(actual,expected,k)
        assert actual.reorder(*p.gens).as_dict()==p.as_dict()
    actual=p.homogeneous_order();expected=op.homogeneous_order()
    assert actual is None if expected is None else sm.to_sympy(actual)==expected
    for gens in [(),(x,),(y,z),(x,y),(z,y,x),(x,x)]:
        assert p.has_only_gens(*gens)==op.has_only_gens(*(sm.to_sympy(g) for g in gens))


@pytest.mark.parametrize('kind',KINDS)
@pytest.mark.parametrize('case',range(4))
def test_dense_views_have_identical_types_and_coefficient_membership(kind,case):
    import sympy as sy
    k,ok,a,b=selected(kind,'grevlex')
    x=sm.Dummy('x',positive=True)
    c=a/(b+1) if kind.endswith('RF') else a if kind.endswith('ring') else sm.sqrt(2) if kind=='AN' else sm.Integer(2)
    p=sm.Poly([sm.Integer(0),c,c*x**6+3*x**2+c*x,c*x][case],x,domain=k)
    op=sy.Poly(sm.to_sympy(p.as_expr()),sm.to_sympy(x),domain=ok)
    assert p.all_monoms()==op.all_monoms()
    assert [sm.to_sympy(c) for c in p.all_coeffs()]==op.all_coeffs()
    assert [(m,sm.to_sympy(c)) for m,c in p.all_terms()]==op.all_terms()
    assert all(isinstance(m,tuple) and len(m)==1 for m in p.all_monoms())
    assert all(c in k for c in p.all_coeffs())
    assert p.all_terms()==list(zip(p.all_monoms(),p.all_coeffs()))


@pytest.mark.parametrize('options',[{}, {'wrt':'a,x2'}, {'wrt':['a','x2']}, {'sort':'a>x>y'}, {'sort':['y','a','x']}, {'order':'grlex'}])
def test_reorder_default_wrt_sort_and_exact_large_numeric_suffixes(options):
    import sympy as sy
    gs=sm.symbols('a x2 x10 x999999999999999999999 x1000000000000000000000 y')
    p=sm.Poly(sum((i+1)*g for i,g in enumerate(gs)),*reversed(gs),domain=sm.QQ)
    op=sy.Poly(sm.to_sympy(p.as_expr()),*(sm.to_sympy(g) for g in reversed(gs)),domain=sy.QQ)
    same_poly(p.reorder(**options),op.reorder(**options),sm.QQ)


def test_reorder_uses_generator_identity_and_preserves_coefficient_parameters():
    x1=sm.Dummy('x',real=True);x2=sm.Dummy('x',real=True)
    a=sm.Symbol('a',positive=True);unused=sm.Dummy('unused')
    k=sm.QQ.frac_field(a,unused,order='grevlex')
    p=sm.Poly(a*x1**2+3*x2,x1,x2,domain=k)
    q=p.reorder(x2,x1)
    assert q.gens==(x2,x1) and q.domain==k
    assert q.as_dict()=={(0,2):a,(1,0):sm.Integer(3)}
    assert p.reorder(wrt=x2).gens==(x2,x1)
    assert p.reorder(wrt=x1).gens==(x1,x2)
    assert q.reorder(x1,x2).as_dict()==p.as_dict()
    lookalike=sm.Dummy('x',real=True)
    with pytest.raises(sm.GeneratorsError):p.has_only_gens(lookalike)
    with pytest.raises(sm.PolynomialError):p.reorder(x1,lookalike)


def test_reorder_never_introduces_duplicate_generators_as_sympy_114_does():
    import sympy as sy
    x,y=sm.symbols('x y');p=sm.Poly(x+y,x,y)
    corrupted=sy.Poly(sm.to_sympy(x+y),sm.to_sympy(x),sm.to_sympy(y)).reorder(sm.to_sympy(x),sm.to_sympy(x),sm.to_sympy(y))
    assert sy.expand(corrupted.as_expr()-sm.to_sympy(p.as_expr()))!=0
    for gens in [(x,x,y),(x,y,y),(x,x),(x,)]:
        with pytest.raises(sm.PolynomialError):p.reorder(*gens)


def test_exception_types_options_and_dense_bounds():
    from symox.polys.polyerrors import GeneratorsError,PolynomialError
    x,y,z=sm.symbols('x y z');p=sm.Poly(x+y,x,y)
    import sympy
    assert GeneratorsError is sm.GeneratorsError
    assert issubclass(GeneratorsError,PolynomialError)==issubclass(sympy.GeneratorsError,sympy.PolynomialError)
    for g in [z,0,[x]]:
        with pytest.raises(GeneratorsError):p.has_only_gens(g)
    for method in ['all_monoms','all_coeffs','all_terms']:
        with pytest.raises(PolynomialError,match='multivariate'):getattr(p,method)()
        huge=sm.Poly(x**10**9+1,x)
        with pytest.raises(ValueError,match='dense output length'):getattr(huge,method)()
    for opts in [{'unknown':True},{'wrt':'x,'},{'wrt':3},{'order':'invalid'},{'sort':3}]:
        with pytest.raises(sm.OptionError):p.reorder(**opts)
    with pytest.raises(NotImplementedError,match='domain'):p.reorder(domain=sm.QQ)
    # These stay sparse: the large exponent must not be used as a list size.
    d,q=huge.deflate();assert d==(10**9,) and q.as_expr()==x+1
    assert huge.homogeneous_order() is None
    assert huge.terms_gcd()[0]==(0,)
    restored=pickle.loads(pickle.dumps(GeneratorsError('sentinel')))
    assert type(restored) is GeneratorsError and restored.args==('sentinel',)


def test_unicode_decimal_generator_suffixes_match_the_oracle():
    import sympy as sy
    gs=sm.symbols('x٢ x10 x० x０ x0 x९ x٩ x١٢ x y')
    p=sm.Poly(sum((i+1)*g for i,g in enumerate(gs)),*gs)
    op=sy.Poly(sm.to_sympy(p.as_expr()),*(sm.to_sympy(g) for g in gs))
    same_poly(p.reorder(),op.reorder(),sm.ZZ)


def test_default_generator_label_work_is_bounded_but_explicit_permutation_stays_sparse():
    x,y=sm.symbols('x y');f=sm.Function('f');atom=x
    for _ in range(35):atom=f(atom,atom)
    p=sm.Poly(atom+1,atom,y)
    with pytest.raises(ValueError,match='unfolded generator label'):p.reorder()
    q=p.reorder(y,atom)
    assert q.gens==(y,atom) and q.as_dict()=={(0,1):sm.Integer(1),(0,0):sm.Integer(1)}
    big=sm.Symbol('a'*100_001)
    p=sm.Poly(big+y,big,y)
    with pytest.raises(ValueError,match='generator label'):p.reorder()
    assert p.reorder(y,big).gens==(y,big)


def test_legacy_expression_and_inexact_coefficients_are_carried_without_rebuilding():
    x,y,a=sm.symbols('x y a')
    for e in [sm.sin(a)*x**6+x*x,sm.Rational(1,4)*x**6+sm.Float(0.5)*x*x]:
        p=sm.Poly(e,x,y)
        d,q=p.deflate();assert d==(2,1)
        assert inflate_deflated(q,d,p.gens)==sm.expand(p.as_expr())
        assert p.reorder(y,x).as_expr()==p.as_expr()
        assert p.exclude().gens==(x,)
        assert p.terms_gcd()[0]==(2,0)


def test_no_sympy_runtime_and_metadata_pickle():
    script=r'''
import sys,importlib.abc,pickle
class NoSymPy(importlib.abc.MetaPathFinder):
    def find_spec(self,name,path=None,target=None):
        if name=='sympy' or name.startswith('sympy.'):
            raise AssertionError('runtime SymPy import')
sys.meta_path.insert(0,NoSymPy())
import symox as s
x=s.Dummy('x',real=True);y=s.Dummy('y');z=s.Symbol('z')
a=s.Symbol('a',positive=True);b=s.Dummy('unused')
for k in [s.ZZ,s.QQ,s.GF(5),s.QQ.frac_field(a,b,order='grevlex'),s.ZZ.poly_ring(a,b,order='grlex'),s.QQ.algebraic_field(s.I,alias='alpha')]:
    p=s.Poly(3*x**6*y**4+2*x*x*y*y,x,y,z,domain=k)
    for q in [p.deflate()[1],p.terms_gcd()[1],p.exclude(),p.reorder(z,y,x)]:
        r=pickle.loads(pickle.dumps(q));assert r.domain==k and r.gens==q.gens and r.as_dict()==q.as_dict()
    assert p.has_only_gens(x,y)
    assert not p.has_only_gens(x)
    assert p.homogeneous_order() is None
    q=s.Poly(x**3+1,x,domain=k)
    assert q.all_monoms()==[(3,),(2,),(1,),(0,)]
    assert q.all_terms()==list(zip(q.all_monoms(),q.all_coeffs()))
'''
    subprocess.run([sys.executable,'-c',script],check=True,timeout=30,env=os.environ.copy())
