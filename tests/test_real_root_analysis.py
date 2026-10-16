"""Exact real-root certificates, without numerical root classification."""
import pickle
import subprocess
import sys

import pytest
import symox as sm

CASES = ['0', '7', 'x', '2*x+1', '2*x**2+1', '(x-1)**3*(x+2)**2',
         '(x-1)**3*(x+2)**2*(x*x-2)', 'x**5-x-1',
         '(3*x-1)*(3*x+1)*(1000000*x-1)*(1000000*x-2)']


def pair(source, domain):
    import sympy as sy
    x = sm.Symbol('x',real=True)
    expr = sm.sympify(source).xreplace({sm.Symbol('x'):x})
    p = sm.Poly(expr,x,domain=domain)
    expected = sy.Poly(sm.to_sympy(expr),sm.to_sympy(x),domain=sy.QQ if domain==sm.QQ else sy.ZZ)
    return p,expected


def compare_sturm(actual, expected):
    import sympy as sy
    assert len(actual)==len(expected)
    for native,oracle in zip(actual,expected):
        expr = native.as_expr() if isinstance(native,sm.Poly) else native
        assert sy.expand(sm.to_sympy(expr)-oracle.as_expr())==0
        if isinstance(native,sm.Poly):
            assert native.domain==sm.QQ
            assert tuple(sm.to_sympy(g) for g in native.gens)==oracle.gens


@pytest.mark.parametrize('source',CASES)
@pytest.mark.parametrize('domain',[sm.ZZ,sm.QQ])
def test_sturm_exact_oracle_and_field_promotion(source,domain):
    import sympy as sy
    p,oracle = pair(source,domain)
    expected = oracle.sturm()
    compare_sturm(p.sturm(),expected)
    compare_sturm(sm.sturm(p),expected)
    compare_sturm(sm.sturm(p.as_expr(),*p.gens,domain=domain,polys=True),expected)
    compare_sturm(sm.sturm(p.as_expr(),*p.gens,domain=domain,polys=False),expected)
    if domain==sm.QQ:
        compare_sturm(p.sturm(auto=False),oracle.sturm(auto=False))
    else:
        with pytest.raises(sm.DomainError):p.sturm(auto=False)
        with pytest.raises(sm.DomainError):sm.sturm(p,auto=False)
    for native in p.sturm():
        restored=pickle.loads(pickle.dumps(native))
        assert restored.domain==sm.QQ and restored.gens==p.gens and restored.as_expr()==native.as_expr()
    for a,b,c in zip(p.sturm(),p.sturm()[1:],p.sturm()[2:]):
        assert a.rem(b).as_expr()==-c.as_expr()


@pytest.mark.parametrize('source',CASES)
@pytest.mark.parametrize('domain',[sm.ZZ,sm.QQ])
@pytest.mark.parametrize('bounds',[(None,None),(-2,2),(-2,1),(-1,0),(0,1),(1,1),(sm.Rational(1,3),sm.Rational(3,2)),(-sm.oo,sm.oo)])
def test_closed_counts_match_exact_oracle_and_ignore_multiplicity(source,domain,bounds):
    p,oracle=pair(source,domain)
    a,b=bounds
    expected=int(oracle.count_roots(None if a is None else sm.to_sympy(a),None if b is None else sm.to_sympy(b)))
    assert p.count_roots(inf=a,sup=b)==expected
    assert sm.count_roots(p,a,b)==expected
    if not p.is_zero and int(p.degree())>0:
        assert sm.count_roots(p.as_expr(),a,b)==expected


def certificate_intervals(p,oracle,actual,epsilon=None,sqf=False,lower=None,upper=None):
    import sympy as sy
    intervals = [(value,1) for value in actual] if sqf else actual
    expected_total=oracle.count_roots(None if lower is None else sm.to_sympy(lower),None if upper is None else sm.to_sympy(upper))
    assert len(intervals)==expected_total
    previous=None
    for (a,b),multiplicity in intervals:
        assert isinstance(a,sm.Expr) and isinstance(b,sm.Expr) and a.is_Rational and b.is_Rational
        assert a<=b
        if previous is not None:assert previous<a
        previous=b
        assert p.count_roots(a,b)==1==oracle.count_roots(sm.to_sympy(a),sm.to_sympy(b))
        if lower is not None:assert a>=lower
        if upper is not None:assert b<=upper
        if epsilon is not None:assert sm.to_sympy(b-a)<sy.Rational(sm.to_sympy(epsilon))
        if not sqf:
            matched=[m for factor,m in oracle.sqf_list()[1] if factor.count_roots(sm.to_sympy(a),sm.to_sympy(b))]
            assert matched==[multiplicity]


@pytest.mark.parametrize('source',CASES)
@pytest.mark.parametrize('domain',[sm.ZZ,sm.QQ])
@pytest.mark.parametrize('fast',[False,True])
@pytest.mark.parametrize('sqf',[False,True])
def test_isolation_is_complete_sorted_disjoint_and_has_exact_multiplicities(source,domain,fast,sqf):
    p,oracle=pair(source,domain)
    options={'fast':fast,'sqf':sqf,'eps':sm.Rational(1,1000)}
    actual=p.intervals(**options)
    certificate_intervals(p,oracle,actual,epsilon=sm.Rational(1,1000),sqf=sqf)
    assert sm.intervals(p,**options)==actual
    if not p.is_zero and int(p.degree())>0:
        assert sm.intervals(p.as_expr(),**options)==actual


@pytest.mark.parametrize('bounds',[(-2,1),(1,1),(sm.Rational(1,3),2),(None,0),(0,None),(10,20)])
@pytest.mark.parametrize('fast',[False,True])
def test_isolation_respects_closed_bounds(bounds,fast):
    p,oracle=pair('(x-1)**3*(x+2)**2*(x*x-2)',sm.ZZ)
    lower,upper=bounds
    actual=p.intervals(inf=lower,sup=upper,fast=fast)
    certificate_intervals(p,oracle,actual,lower=lower,upper=upper)


def test_rational_singletons_and_irrational_high_degree_root():
    p,_=pair('(3*x-1)*(3*x+1)*(1000000*x-1)*(1000000*x-2)',sm.QQ)
    assert p.intervals()==[((sm.Rational(-1,3),sm.Rational(-1,3)),1),((sm.Rational(1,1000000),sm.Rational(1,1000000)),1),((sm.Rational(1,500000),sm.Rational(1,500000)),1),((sm.Rational(1,3),sm.Rational(1,3)),1)]
    p,oracle=pair('x**5-x-1',sm.ZZ)
    roots=p.intervals(eps=sm.Rational(1,2**150))
    certificate_intervals(p,oracle,roots,epsilon=sm.Rational(1,2**150))
    assert roots[0][0][0]>1 and roots[0][0][1]<2
    assert sm.intervals(0)==[]==sm.intervals(7)


@pytest.mark.parametrize('source,bounds',[('x*x-2',(1,2)),('x*x-2',(-2,-1)),('x**5-x-1',(1,2)),('(x-1)**3',(0,2)),('3*x-1',(sm.Rational(1,3),sm.Rational(1,3)))])
@pytest.mark.parametrize('options',[{}, {'steps':0}, {'steps':3}, {'eps':sm.Rational(1,10**10)}, {'eps':sm.Rational(1,2**150),'fast':True}, {'eps':sm.Rational(1,2**150),'steps':1}, {'steps':3,'fast':True}])
def test_refinement_keeps_one_root_and_obeys_precision_or_iteration_limit(source,bounds,options):
    p,oracle=pair(source,sm.QQ)
    lower,upper=bounds
    a,b=p.refine_root(lower,upper,**options)
    assert lower<=a<=b<=upper
    assert p.count_roots(a,b)==1==oracle.count_roots(sm.to_sympy(a),sm.to_sympy(b))
    if options.get('steps')==0:assert (a,b)==bounds
    if 'eps' in options and 'steps' not in options:assert b-a<options['eps']
    assert sm.refine_root(p,lower,upper,**options)==(a,b)
    assert sm.refine_root(p.as_expr(),lower,upper,**options)==(a,b)


def test_refinement_with_float_precision_and_sub_float_resolution():
    x=sm.Symbol('x')
    p=sm.Poly(x*x-2,x)
    a,b=p.refine_root(1,2,eps=1e-8)
    assert p.count_roots(a,b)==1 and float(b-a)<1e-8
    a,b=p.refine_root(1,2,eps=sm.Rational(1,2**250))
    assert p.count_roots(a,b)==1 and b-a<sm.Rational(1,2**250)
    assert float(a)==float(b)
    assert p.count_roots(sm.Rational(1,10),sm.Rational(2,10))==0


def test_closed_count_additivity_with_shared_rational_root():
    p,_=pair('(x-1)**3*(x+2)**2*(x*x-2)',sm.ZZ)
    for left,middle,right in [(-3,1,3),(-3,-2,3),(-2,0,2),(1,1,2)]:
        assert p.count_roots(left,right)==p.count_roots(left,middle)+p.count_roots(middle,right)-p.count_roots(middle,middle)


@pytest.mark.parametrize('source', ['0','7','x*x-2'])
def test_invalid_intervals_fail_instead_of_guessing_or_copying_oracle_bugs(source):
    p,_=pair(source,sm.ZZ)
    with pytest.raises(ValueError,match='increasing order'):p.count_roots(1,0)
    with pytest.raises(ValueError,match='increasing order'):p.intervals(inf=1,sup=0)
    with pytest.raises(ValueError,match='increasing order'):p.refine_root(1,0)
    for epsilon in [0,-1]:
        with pytest.raises(ValueError,match='positive rational'):p.intervals(eps=epsilon)
        with pytest.raises(ValueError,match='positive rational'):p.refine_root(0,2,eps=epsilon)
    with pytest.raises(sm.RefinementFailed) as failure:p.refine_root(1,1)
    assert type(pickle.loads(pickle.dumps(failure.value))) is sm.RefinementFailed


def test_refinement_requires_a_unique_root_and_squarefree_when_requested():
    p,_=pair('x*x-2',sm.QQ)
    with pytest.raises(sm.RefinementFailed):p.refine_root(-2,2)
    with pytest.raises(sm.RefinementFailed):p.refine_root(2,3)
    p,_=pair('(x-1)**3',sm.QQ)
    with pytest.raises(sm.RefinementFailed):p.refine_root(0,2,check_sqf=True)
    assert p.refine_root(1,1)==(1,1)


@pytest.mark.parametrize('method', ['sturm','count_roots','intervals','refine_root'])
def test_unsupported_domains_and_complex_options_are_explicit(method):
    x,a,y=sm.symbols('x a y')
    args=(1,2) if method=='refine_root' else ()
    for p in [sm.Poly(x*x+1,x,modulus=5),sm.Poly(x*x+1,x,domain=sm.QQ.algebraic_field(sm.sqrt(2))),sm.Poly(x*x+a,x,domain=sm.QQ.frac_field(a)),sm.Poly(x*y+1,x,y),sm.Poly(x*x+sm.sin(a),x),sm.Poly(0.5*x+1,x)]:
        with pytest.raises(sm.DomainError):getattr(p,method)(*args)
    p=sm.Poly(x*x+1,x)
    with pytest.raises(NotImplementedError):p.intervals(all=True)
    with pytest.raises(NotImplementedError):sm.intervals([x*x+1,x*x-2])
    with pytest.raises(NotImplementedError):sm.intervals(x*x+1,strict=True)
    with pytest.raises(NotImplementedError):p.count_roots(0,1+sm.I)
    with pytest.raises(sm.CoercionFailed):p.count_roots(sm.sqrt(2),2)


def test_sturm_options_and_metadata_are_applied():
    x=sm.Dummy('x',positive=True)
    p=sm.Poly((x-1)**2*(x+2),x,domain=sm.ZZ)
    assert all(v.gens==(x,) and v.domain is sm.QQ for v in p.sturm())
    assert all(v.domain is sm.QQ for v in sm.sturm(p,domain=sm.QQ,polys=True,auto=False))
    with pytest.raises(sm.OptionError):sm.sturm(p,bogus=True)
    with pytest.raises(sm.OptionError):sm.sturm(p,auto='bad')
    with pytest.raises(sm.DomainError):sm.sturm(p,modulus=5)
    with pytest.raises(TypeError):p.count_roots(eps=sm.Rational(1,10))


def test_degree_and_coefficient_limits_reject_before_large_allocations():
    x=sm.Symbol('x')
    p=sm.Poly(x**1_000_000_000+1,x)
    for method,args in [('sturm',()),('count_roots',()),('intervals',()),('refine_root',(0,2))]:
        with pytest.raises(ValueError,match='resource limit'):getattr(p,method)(*args)
    p=sm.Poly(sm.Integer(1<<262144)*x+1,x)
    with pytest.raises(ValueError,match='resource limit'):p.count_roots()


def test_root_analysis_without_sympy_runtime():
    code='''
import importlib.abc,sys,pickle
class Block(importlib.abc.MetaPathFinder):
 def find_spec(self,name,path=None,target=None):
  if name=='sympy' or name.startswith('sympy.'):raise ImportError('blocked')
sys.meta_path.insert(0,Block())
import symox as s
x=s.Dummy('x',real=True);p=s.Poly((x-1)**3*(x+2)**2*(x*x-2),x)
assert p.count_roots()==4 and s.count_roots(p,-2,1)==3
assert all(q.gens==(x,) and q.domain==s.QQ for q in p.sturm())
roots=s.intervals(p,eps=s.Rational(1,1000));assert len(roots)==4
assert sorted(m for _,m in roots)==[1,1,2,3]
for (a,b),m in roots:assert p.count_roots(a,b)==1 and b-a<s.Rational(1,1000)
q=s.Poly(x**5-x-1,x);a,b=q.refine_root(1,2,eps=s.Rational(1,2**150))
assert s.count_roots(q,a,b)==1 and b-a<s.Rational(1,2**150)
try:q.refine_root(1,1)
except s.RefinementFailed as e:assert type(pickle.loads(pickle.dumps(e))) is s.RefinementFailed
else:raise AssertionError('invalid singleton accepted')
assert not any(n=='sympy' or n.startswith('sympy.') for n in sys.modules)
'''
    result=subprocess.run([sys.executable,'-c',code],text=True,capture_output=True)
    assert result.returncode==0,result.stderr
