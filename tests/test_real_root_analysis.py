"""Exact real-root counting, isolation and refinement, compared with SymPy."""
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


def outcome(module, call):
    """The printed result of ``call(module)``, or its exception class and text."""
    try:
        value = call(module)
    except Exception as error:  # noqa: BLE001 - both outcomes are compared
        return ('raises', type(error).__name__, str(error))
    return ('returns', str(value))


def same_as_reference(call):
    import sympy as sy
    expected = outcome(sy, call)
    assert outcome(sm, call) == expected
    return expected


def poly(module, source, domain='ZZ', symbol='x'):
    x = module.Symbol(symbol, real=True)
    expr = module.sympify(source, locals={'x': x})
    return module.Poly(expr, x, domain=domain)


@pytest.mark.parametrize('source',CASES)
@pytest.mark.parametrize('domain',['ZZ','QQ'])
@pytest.mark.parametrize('fast',[False,True])
@pytest.mark.parametrize('sqf',[False,True])
def test_isolation_matches_reference(source,domain,fast,sqf):
    eps = lambda m: m.Rational(1,1000)
    same_as_reference(lambda m: poly(m,source,domain).intervals(fast=fast,sqf=sqf,eps=eps(m)))
    same_as_reference(lambda m: m.intervals(poly(m,source,domain),fast=fast,sqf=sqf,eps=eps(m)))
    same_as_reference(lambda m: m.intervals(poly(m,source,domain).as_expr(),fast=fast,sqf=sqf,eps=eps(m)))


@pytest.mark.parametrize('bounds',[(-2,1),(1,1),((1,3),2),(None,0),(0,None),(10,20),(1,0),(2,-2)])
@pytest.mark.parametrize('fast',[False,True])
def test_isolation_bounds_match_reference(bounds,fast):
    def bound(m, b):
        return m.Rational(*b) if isinstance(b,tuple) else b
    lower,upper=bounds
    same_as_reference(lambda m: poly(m,'(x-1)**3*(x+2)**2*(x*x-2)').intervals(
        inf=bound(m,lower),sup=bound(m,upper),fast=fast))


def test_rational_singletons_and_irrational_high_degree_root():
    same_as_reference(lambda m: poly(m,'(3*x-1)*(3*x+1)*(1000000*x-1)*(1000000*x-2)','QQ').intervals())
    same_as_reference(lambda m: poly(m,'x**5-x-1').intervals(eps=m.Rational(1,2**150)))
    same_as_reference(lambda m: (m.intervals(0), m.intervals(7)))


@pytest.mark.parametrize('source,bounds',[('x*x-2',(1,2)),('x*x-2',(-2,-1)),('x**5-x-1',(1,2)),('(x-1)**3',(0,2)),('3*x-1',((1,3),(1,3)))])
@pytest.mark.parametrize('options',[{}, {'steps':0}, {'steps':3}, {'eps':(1,10**10)}, {'eps':(1,2**150),'fast':True}, {'eps':(1,2**150),'steps':1}, {'steps':3,'fast':True}])
def test_refinement_matches_reference(source,bounds,options):
    def value(m, b):
        return m.Rational(*b) if isinstance(b,tuple) else b
    def call(m, f):
        kwargs={k: value(m,v) for k,v in options.items()}
        lower,upper=(value(m,b) for b in bounds)
        return f(m)(lower,upper,**kwargs)
    same_as_reference(lambda m: call(m, lambda m: poly(m,source,'QQ').refine_root))
    same_as_reference(lambda m: call(m, lambda m: lambda a,b,**k: m.refine_root(poly(m,source,'QQ'),a,b,**k)))
    same_as_reference(lambda m: call(m, lambda m: lambda a,b,**k: m.refine_root(poly(m,source,'QQ').as_expr(),a,b,**k)))


def test_refinement_with_float_precision_and_sub_float_resolution():
    x=sm.Symbol('x')
    p=sm.Poly(x*x-2,x)
    a,b=p.refine_root(1,2,eps=1e-8)
    assert p.count_roots(a,b)==1 and float(b-a)<1e-8
    a,b=p.refine_root(1,2,eps=sm.Rational(1,2**250))
    assert p.count_roots(a,b)==1 and b-a<sm.Rational(1,2**250)
    assert float(a)==float(b)
    assert p.count_roots(sm.Rational(1,10),sm.Rational(2,10))==0
    same_as_reference(lambda m: m.Poly(m.Symbol('x')**2-2).refine_root(1,2,eps=1e-8))


def test_closed_count_additivity_with_shared_rational_root():
    p,_=pair('(x-1)**3*(x+2)**2*(x*x-2)',sm.ZZ)
    for left,middle,right in [(-3,1,3),(-3,-2,3),(-2,0,2),(1,1,2)]:
        assert p.count_roots(left,right)==p.count_roots(left,middle)+p.count_roots(middle,right)-p.count_roots(middle,middle)


@pytest.mark.parametrize('source', ['0','7','x*x-2'])
@pytest.mark.parametrize('call', [
    lambda p, m: p.count_roots(1,0),
    lambda p, m: p.intervals(inf=1,sup=0),
    lambda p, m: p.refine_root(1,0),
    lambda p, m: p.refine_root(2,1),
    lambda p, m: p.intervals(eps=0),
    lambda p, m: p.intervals(eps=-1),
    lambda p, m: p.refine_root(0,2,eps=0),
    lambda p, m: p.refine_root(0,2,eps=-1),
    lambda p, m: p.refine_root(1,1),
])
def test_degenerate_intervals_match_reference(source,call):
    same_as_reference(lambda m: call(poly(m,source),m))


def test_refinement_failure_pickles():
    p,_=pair('x*x-2',sm.ZZ)
    with pytest.raises(sm.RefinementFailed) as failure:p.refine_root(2,3)
    assert type(pickle.loads(pickle.dumps(failure.value))) is sm.RefinementFailed


@pytest.mark.parametrize('call', [
    lambda m: poly(m,'x*x-2','QQ').refine_root(-2,2),
    lambda m: poly(m,'x*x-2','QQ').refine_root(2,3),
    lambda m: poly(m,'(x-1)**3','QQ').refine_root(0,2,check_sqf=True),
    lambda m: poly(m,'(x-1)**3','QQ').refine_root(0,2),
    lambda m: poly(m,'(x-1)**3','QQ').refine_root(1,1),
])
def test_refinement_uniqueness_matches_reference(call):
    same_as_reference(call)


def _polys(m):
    x,a,y=m.symbols('x a y')
    return [m.Poly(x*x+1,x,modulus=5),
            m.Poly(x*x+1,x,domain=m.QQ.algebraic_field(m.sqrt(2))),
            m.Poly(x*x-2,x,domain=m.QQ.algebraic_field(m.sqrt(2))),
            m.Poly(x*x+a,x,domain=m.QQ.frac_field(a)),
            m.Poly(x*y+1,x,y),
            m.Poly(x*x+m.sin(a),x),
            m.Poly(0.5*x+1,x),
            m.Poly(0.5*x*x-1,x)]


@pytest.mark.parametrize('method,index', [
    (method, index) for method in ['sturm','count_roots','intervals','refine_root'] for index in range(8)])
def test_restricted_domains_match_reference(index,method):
    args=(1,2) if method=='refine_root' else ()
    same_as_reference(lambda m: getattr(_polys(m)[index],method)(*args))


@pytest.mark.parametrize('call', [
    lambda m, x: m.Poly(x*x+1,x).intervals(all=True),
    lambda m, x: m.Poly(x**3-x+1,x).intervals(all=True,eps=m.Rational(1,100)),
    lambda m, x: m.intervals([x*x+1,x*x-2]),
    lambda m, x: m.intervals(x*x+1,strict=True),
    lambda m, x: m.Poly(x*x+1,x).count_roots(0,1+m.I),
    lambda m, x: m.Poly(x*x+1,x).count_roots(-m.I,m.I),
    lambda m, x: m.Poly(x*x-2,x).count_roots(m.sqrt(2),2),
    lambda m, x: m.Poly(x*x-2,x).count_roots(m.Rational(1,10),0.5),
    lambda m, x: m.Poly(x*x-2,x).count_roots(m.oo,2),
])
def test_complex_and_list_options_match_reference(call):
    same_as_reference(lambda m: call(m, m.Symbol('x')))


@pytest.mark.parametrize('auto', [True, False])
@pytest.mark.parametrize('build', [
    lambda m, x: m.Poly((x-1)**5*(x+2),x,modulus=5),
    lambda m, x: m.Poly(x**10-1,x,modulus=5),
    lambda m, x: m.Poly(0,x,modulus=5),
    lambda m, x: m.Poly((x*x-2)**2*(x-1),x,domain=m.QQ.algebraic_field(m.sqrt(2))),
    lambda m, x: m.Poly(x**3+m.Symbol('a')*x,x,domain=m.QQ[m.Symbol('a')]),
    lambda m, x: m.Poly(x*m.Symbol('y')+1,x,m.Symbol('y'),domain='QQ'),
    lambda m, x: m.Poly(0.1*x**3-0.7*x+0.3,x),
    lambda m, x: m.Poly((0.1*x-0.3)**2*(x+0.7),x),
    lambda m, x: m.Poly(m.Float(1,30)*x**3/3-x+1,x),
    lambda m, x: m.Poly(0,x),
    lambda m, x: m.Poly(7,x),
    lambda m, x: m.Poly(7,x,domain='QQ'),
])
def test_sturm_matches_reference(build,auto):
    same_as_reference(lambda m: build(m, m.Symbol('x')).sturm(auto=auto))
    same_as_reference(lambda m: m.sturm(build(m, m.Symbol('x')), auto=auto))


def test_sturm_options_and_metadata_are_applied():
    x=sm.Dummy('x',positive=True)
    p=sm.Poly((x-1)**2*(x+2),x,domain=sm.ZZ)
    assert all(v.gens==(x,) and v.domain is sm.QQ for v in p.sturm())
    assert all(v.domain is sm.QQ for v in sm.sturm(p,domain=sm.QQ,polys=True,auto=False))
    for call in [lambda m, p: m.sturm(p,bogus=True), lambda m, p: m.sturm(p,auto='bad'),
                 lambda m, p: m.sturm(p,auto=1), lambda m, p: m.sturm(p,modulus=5),
                 lambda m, p: p.count_roots(eps=m.Rational(1,10))]:
        same_as_reference(lambda m: call(m, m.Poly((m.Symbol('x')-1)**2*(m.Symbol('x')+2),domain='ZZ')))


def test_degree_and_coefficient_limits_reject_before_large_allocations():
    x=sm.Symbol('x')
    p=sm.Poly(x**1_000_000_000+1,x)
    # The reference would allocate a dense list of a billion coefficients.
    for method,args in [('sturm',()),('count_roots',()),('intervals',()),('refine_root',(0,2))]:
        with pytest.raises(ValueError,match='resource limit'):getattr(p,method)(*args)
    # Large coefficients are computed exactly, as in the reference.
    p=sm.Poly(sm.Integer(1<<262144)*x+1,x)
    assert p.count_roots()==1


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
assert q.refine_root(1,1)==(1,1)
try:q.refine_root(2,3)
except s.RefinementFailed as e:assert type(pickle.loads(pickle.dumps(e))) is s.RefinementFailed
else:raise AssertionError('interval without a root accepted')
assert not any(n=='sympy' or n.startswith('sympy.') for n in sys.modules)
'''
    result=subprocess.run([sys.executable,'-c',code],text=True,capture_output=True)
    assert result.returncode==0,result.stderr
