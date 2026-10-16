"""Polynomial antiderivatives must stay in, or explicitly extend, their domain."""
import pickle
import subprocess
import sys

import pytest
import symox as sm

KINDS=['ZZ','QQ','ZZring','QQring','ZZRF','QQRF','GF','AN']


def context(module,kind):
    x,y=module.symbols('x y')
    a=module.Symbol('a',positive=True)
    if kind=='ZZ':k=module.ZZ;c=module.Integer(1)
    elif kind=='QQ':k=module.QQ;c=module.Integer(1)
    elif kind=='GF':k=module.GF(5);c=module.Integer(1)
    elif kind=='AN':k=module.QQ.algebraic_field(module.sqrt(2),alias='alpha');c=module.sqrt(2)
    else:
        ground=module.ZZ if kind.startswith('ZZ') else module.QQ
        k=(ground.poly_ring if kind.endswith('ring') else ground.frac_field)(a,order='grevlex');c=a
    return x,y,k,c


@pytest.mark.parametrize('kind',KINDS)
@pytest.mark.parametrize('auto',[True,False])
@pytest.mark.parametrize('expression',['zero','exact','fraction','multiple'])
def test_exact_domains_match_sympy_and_every_coefficient_belongs(kind,auto,expression):
    sy=pytest.importorskip('sympy')
    x,y,k,c=context(sm,kind);ox,oy,ok,oc=context(sy,kind)
    def expr(name,x,y,c):
        return {'zero':0,'exact':2*c*x,'fraction':c*x,'multiple':12*c*x*y+6*c*x*x*y*y}[name]
    p=sm.Poly(expr(expression,x,y,c),x,y,domain=k)
    op=sy.Poly(expr(expression,ox,oy,oc),ox,oy,domain=ok)
    try:expected=op.integrate(auto=auto)
    except Exception as failure:
        with pytest.raises(Exception) as native:p.integrate(auto=auto)
        assert type(native.value).__name__==type(failure).__name__
        return
    # SymPy 1.14 truncates nonexact multivariate integer integration. Keep
    # the independent derivative certificate and require an exact exception.
    if sy.expand(expected.diff(ox).as_expr()-op.as_expr())!=0:
        assert not auto and kind in ['ZZ','ZZring']
        with pytest.raises(sm.ExactQuotientFailed):p.integrate(auto=auto)
        return
    actual=p.integrate(auto=auto)
    assert str(actual.domain)==str(expected.domain)
    assert sy.simplify(sm.to_sympy(actual.as_expr())-expected.as_expr())==0
    for coefficient in actual.as_dict().values():
        value=actual.domain.convert(coefficient)
        assert value.parent()==actual.domain
    assert sm.Poly(actual.diff(x).as_expr()-p.as_expr(),x,y,domain=actual.domain).is_zero
    assert pickle.loads(pickle.dumps(actual)).domain==actual.domain
    if kind.endswith(('RF','ring')):assert actual.domain.order=='grevlex'


@pytest.mark.parametrize('kind',KINDS)
@pytest.mark.parametrize('auto',[True,False])
@pytest.mark.parametrize('specs',['repeated','two_generators','zero','negative','indices'])
def test_tuple_orders_and_multiple_generators(kind,auto,specs):
    sy=pytest.importorskip('sympy')
    x,y,k,c=context(sm,kind);ox,oy,ok,oc=context(sy,kind)
    p=sm.Poly(12*c*x*y,x,y,domain=k);op=sy.Poly(12*oc*ox*oy,ox,oy,domain=ok)
    spec=lambda x,y:{'repeated':[(x,2)],'two_generators':[(x,2),(y,1)],'zero':[(x,0)],'negative':[(x,-1)],'indices':[(0,1),(-1,1)]}[specs]
    try:expected=op.integrate(*spec(ox,oy),auto=auto)
    except Exception as failure:
        with pytest.raises(Exception) as native:p.integrate(*spec(x,y),auto=auto)
        assert type(native.value).__name__==type(failure).__name__
        return
    actual=p.integrate(*spec(x,y),auto=auto)
    assert str(actual.domain)==str(expected.domain)
    assert sy.simplify(sm.to_sympy(actual.as_expr())-expected.as_expr())==0
    assert all(actual.domain.convert(c).parent()==actual.domain for c in actual.as_dict().values())


@pytest.mark.parametrize('auto',[True,False])
def test_positive_characteristic_missing_antiderivative_is_not_a_panic(auto):
    x=sm.Symbol('x');k=sm.GF(5)
    p=sm.Poly(x**4+2*x,x,domain=k)
    with pytest.raises(sm.NotInvertible,match='zero'):p.integrate(auto=auto)
    q=sm.Poly(x,x,domain=k).integrate(auto=auto)
    assert q.as_expr()==-2*x*x and q.domain==k
    assert q.diff()==sm.Poly(x,x,domain=k)


def test_failed_ring_integration_is_transactional_and_preserves_metadata():
    x,y=sm.symbols('x y');a=sm.Symbol('a',positive=True);b=sm.Symbol('b',negative=True)
    k=sm.ZZ.poly_ring(a,b,order='grevlex')
    p=sm.Poly(a*x,x,y,domain=k)
    with pytest.raises(sm._core.ExactQuotientFailed):p.integrate(auto=False)
    assert p.domain==k and p.as_expr()==a*x
    q=p.integrate(auto=True)
    assert q.domain==k.get_field() and q.domain.symbols==(a,b)
    assert q.as_expr()==a*x*x/2
    assert q.domain.convert(q.coeff_monomial(x*x)).parent()==q.domain
    # A later operation can fail after an earlier one succeeded, with no mutation.
    with pytest.raises(sm._core.ExactQuotientFailed):sm.Poly(2*a*x*y,x,y,domain=k).integrate(x,y,auto=False)
    for spec in [(x,1,2),(3,1),(-3,1)]:
        with pytest.raises(ValueError):p.integrate(spec)


def test_native_domain_integration_without_sympy():
    code=r'''
import sys,pickle
class NoSymPy:
    def find_spec(self,name,*args):
        if name=='sympy' or name.startswith('sympy.'):raise RuntimeError('runtime SymPy')
sys.meta_path.insert(0,NoSymPy())
import symox as s
x=s.Symbol('x');a=s.Symbol('a',positive=True);b=s.Symbol('b',negative=True)
k=s.ZZ.poly_ring(a,b,order='grevlex');p=s.Poly(a*x,x,domain=k)
try:p.integrate(auto=False)
except s._core.ExactQuotientFailed:pass
else:raise AssertionError('nonexact ring integration accepted')
q=p.integrate();assert q.domain==k.get_field() and q.domain.symbols==(a,b)
assert q.diff()==p.set_domain(q.domain)
assert pickle.loads(pickle.dumps(q)).domain==q.domain
k=s.QQ.algebraic_field(s.sqrt(2),alias='alpha');p=s.Poly(s.sqrt(2)*x,x,domain=k)
assert p.integrate(auto=False).domain==k
assert 'sympy' not in sys.modules
'''
    subprocess.run([sys.executable,'-c',code],check=True)


def test_exact_quotient_exception_is_exported_and_pickleable():
    from symox.polys.domains import ExactQuotientFailed
    from symox.polys.polyerrors import ExactQuotientFailed as FromErrors
    assert ExactQuotientFailed is FromErrors is sm.ExactQuotientFailed
    original=ExactQuotientFailed('nonexact coefficient division')
    restored=pickle.loads(pickle.dumps(original))
    assert type(restored) is type(original) and restored.args==original.args
