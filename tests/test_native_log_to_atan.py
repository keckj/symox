"""Exact Euclidean representatives and principal-atan branch compatibility."""
import importlib
import random

import pytest
import sympy as sp
import symox as sx
from symox import Poly,QQ,ZZ

rt=importlib.import_module('symox.integrals.rationaltools')
ref=importlib.import_module('sympy.integrals.rationaltools')
x=sx.Symbol('x',real=True)
xp=sp.Symbol('x',real=True)


def pair(a,b,domain=QQ):
    f=Poly.from_list([sx.Rational(str(c)) for c in reversed(a)],x,domain=domain)
    g=Poly.from_list([sx.Rational(str(c)) for c in reversed(b)],x,domain=domain)
    return f,g


def expected(a,b):
    f=sp.Poly.from_list(list(reversed(a)),xp,domain=sp.QQ)
    g=sp.Poly.from_list(list(reversed(b)),xp,domain=sp.QQ)
    return ref.log_to_atan(f,g),f,g


@pytest.mark.parametrize('domain',[ZZ,QQ])
@pytest.mark.parametrize('seed',range(32))
def test_exact_sum_matches_reference_recurrence_including_branch_constants(domain,seed,monkeypatch):
    rng=random.Random(seed+48109)
    m,n=seed%5,(seed*3+1)%5
    a=[sp.Rational(rng.randrange(-5,6),rng.randrange(1,6) if domain==QQ else 1) for _ in range(m+1)]
    b=[sp.Rational(rng.randrange(-5,6),rng.randrange(1,6) if domain==QQ else 1) for _ in range(n+1)]
    a[-1]=a[-1] or 1;b[-1]=b[-1] or -1
    if seed%8==0:
        # A shared factor must preserve the same Bezout representative.
        a=[0]+a;b=[0]+b
    f,g=pair(a,b,domain);value,fp,gp=expected(a,b)
    native=sx._core._log_to_atan(f,g)
    assert native is not None
    reference=sx.sympify(str(value),locals={'x':x})
    assert native==reference
    calls=[];old=sx._core._log_to_atan
    def selected(*args):
        result=old(*args);calls.append(result);return result
    monkeypatch.setattr(sx._core,'_log_to_atan',selected)
    assert rt.log_to_atan(f,g)==reference and len(calls)==1 and calls[0] is not None
    # Disable the route to compare the entire previously supported algorithm.
    monkeypatch.setattr(sx._core,'_log_to_atan',lambda *args:None)
    assert rt.log_to_atan(f,g)==native
    # Exact rational derivative is an independent certificate, in addition to
    # expression equality (derivative equality alone permits branch constants).
    target=2*(gp.as_expr()*sp.diff(fp.as_expr(),xp)-fp.as_expr()*sp.diff(gp.as_expr(),xp))/(fp.as_expr()**2+gp.as_expr()**2)
    assert sp.cancel(sp.diff(value,xp)-target)==0
    if seed%8==0:
        for point in (-sp.Rational(7,3),-1,sp.Rational(1,5),2):
            p=native.subs(x,sx.Rational(str(point))).evalf(40)
            q=value.subs(xp,point).evalf(40)
            assert abs(sp.N(sp.sympify(str(p))-q,30))<sp.Rational(1,10**30)


@pytest.mark.parametrize('a,b',[
    ([0,1],[1,1]),([1,0,1],[-1,1]),([1,1],[1,0,1]),
    ([0,0,1],[0,1]),([-1,-1],[-1,0,-1]),([2],[3]),([1,2],[2,4]),
])
def test_branch_sensitive_examples_match_at_all_real_sides(a,b):
    f,g=pair(a,b);value,_,_=expected(a,b)
    actual=rt.log_to_atan(f,g)
    assert actual==sx.sympify(str(value),locals={'x':x})
    for point in (-4,-1,0,1,4):
        q=value.subs(xp,point).evalf(45)
        p=actual.subs(x,point).evalf(45)
        assert abs(sp.N(sp.sympify(str(p))-q,35))<sp.Rational(1,10**35)


@pytest.mark.parametrize('a,b',[([],[]),([],[1]),([1],[])])
def test_zero_inputs_keep_the_reference_division_error(a,b):
    f,g=pair(a,b)
    assert sx._core._log_to_atan(f,g) is None
    with pytest.raises(ZeroDivisionError,match='^polynomial division by zero$'):rt.log_to_atan(f,g)


@pytest.mark.parametrize('domain',['EX','RR'])
def test_nonrational_domains_keep_the_complete_original_path(domain,monkeypatch):
    if domain=='EX':f,g=Poly(x+sx.Rational(1,2),x,domain=domain),Poly(sx.sqrt(3)/2,x,domain=domain)
    else:f,g=Poly(x+0.5,x,domain=domain),Poly(2.0,x,domain=domain)
    assert sx._core._log_to_atan(f,g) is None
    actual=rt.log_to_atan(f,g)
    monkeypatch.setattr(sx._core,'_log_to_atan',lambda *args:None)
    assert rt.log_to_atan(f,g)==actual


def test_different_generators_and_multivariate_polynomials_decline_without_conversion():
    y=sx.Symbol('y')
    assert sx._core._log_to_atan(Poly(x,x),Poly(y,y)) is None
    assert sx._core._log_to_atan(Poly(x+y,x,y),Poly(1,x,y)) is None
    class Protocol:
        def __getattribute__(self,name):raise AssertionError('eligibility invoked a protocol')
    assert sx._core._log_to_atan(Protocol(),Protocol()) is None


@pytest.mark.parametrize('method',['degree','to_field','div','gcdex','quo','as_expr','__neg__','__mul__','__add__'])
def test_real_poly_overrides_are_observed_and_refuse_native(method,monkeypatch):
    f,g=pair([1,0,1],[-1,1]);original=getattr(Poly,method);calls=[]
    def callback(*args,**kwargs):calls.append(method);return original(*args,**kwargs)
    monkeypatch.setattr(Poly,method,callback)
    assert sx._core._log_to_atan(f,g) is None
    value=rt.log_to_atan(f,g)
    assert calls and value==2*sx.atan(x+1)-2*sx.atan(x+x**2/2+x**3/2)


def test_custom_atan_and_recursive_callbacks_are_preserved(monkeypatch):
    f,g=pair([1,0,1],[-1,1]);calls=[];original=rt.atan
    def atan(arg):calls.append(arg);return original(arg)
    monkeypatch.setattr(rt,'atan',atan)
    rt.log_to_atan(f,g)
    assert len(calls)==2
    monkeypatch.setattr(rt,'atan',original)
    calls.clear();original_recursion=rt.log_to_atan
    def recurse(*args):calls.append(args);return original_recursion(*args)
    monkeypatch.setattr(rt,'log_to_atan',recurse)
    original_recursion(f,g)
    assert len(calls)==1


def test_held_evaluation_keeps_the_old_unevaluated_construction(monkeypatch):
    f,g=pair([0,1],[1,1]);old=sx._core._log_to_atan
    def forbidden(*args):raise AssertionError('held context entered evaluated kernel')
    monkeypatch.setattr(sx._core,'_log_to_atan',forbidden)
    with sx.evaluate(False):
        actual=rt.log_to_atan(f,g)
    monkeypatch.setattr(sx._core,'_log_to_atan',lambda *args:None)
    with sx.evaluate(False):
        expected=rt.log_to_atan(f,g)
    assert actual==expected


def test_interruption_from_real_poly_callback_is_not_a_successful_refusal(monkeypatch):
    f,g=pair([1,0,1],[-1,1])
    def interrupted(*args,**kwargs):raise KeyboardInterrupt('gcdex')
    monkeypatch.setattr(Poly,'gcdex',interrupted)
    with pytest.raises(KeyboardInterrupt,match='gcdex'):rt.log_to_atan(f,g)


def test_reloading_the_adapter_keeps_stock_dispatch_and_rejects_real_overrides(monkeypatch):
    import functools
    f,g=pair([1,0,1],[-1,1]);original=Poly.as_expr
    importlib.reload(rt)
    assert rt._stock_atan_protocol()
    calls=[];helper=sx._core._log_to_atan
    def selected(*args):calls.append(True);return helper(*args)
    monkeypatch.setattr(sx._core,'_log_to_atan',selected)
    rt.log_to_atan(f,g)
    assert calls
    @functools.wraps(original)
    def callback(*args,**kwargs):raise KeyboardInterrupt('wrapped as_expr')
    monkeypatch.setattr(Poly,'as_expr',callback)
    importlib.reload(rt)
    with pytest.raises(KeyboardInterrupt,match='wrapped as_expr'):rt.log_to_atan(f,g)


def test_native_as_expr_delegate_override_keeps_the_old_callback(monkeypatch):
    f,g=pair([0,1],[1]);evaluation=importlib.import_module('symox.polys._evaluation')
    def callback(*args):raise KeyboardInterrupt('as_expr delegate')
    monkeypatch.setattr(evaluation,'_native_as_expr',callback)
    assert sx._core._log_to_atan(f,g) is None
    with pytest.raises(KeyboardInterrupt,match='as_expr delegate'):rt.log_to_atan(f,g)


def test_atan_factory_and_result_arithmetic_hooks_are_not_bypassed(monkeypatch):
    f,g=pair([0,1],[1]);calls=[];original=sx.atan.__dict__['_symox_factory']
    def factory(*args,**kwargs):calls.append(args);return original(*args,**kwargs)
    with monkeypatch.context() as patch:
        patch.setattr(sx.atan,'_symox_factory',factory)
        rt.log_to_atan(f,g)
        assert len(calls)==1
    def rmul(*args):raise KeyboardInterrupt('atan rmul')
    monkeypatch.setattr(type(sx.atan(x)),'__rmul__',rmul)
    with pytest.raises(KeyboardInterrupt,match='atan rmul'):rt.log_to_atan(f,g)


def test_domain_property_override_remains_observable(monkeypatch):
    f,g=pair([0,1],[1]);T=type(QQ)
    monkeypatch.setattr(T,'is_EX',property(lambda *args:(_ for _ in ()).throw(KeyboardInterrupt('is_EX'))))
    assert sx._core._log_to_atan(f,g) is None
    with pytest.raises(KeyboardInterrupt,match='is_EX'):rt.log_to_atan(f,g)


def test_degree_selection_declines_without_losing_the_general_algorithm(monkeypatch):
    f,g=Poly(x**129+1,x,domain=QQ),Poly(1,x,domain=QQ)
    assert sx._core._log_to_atan(f,g) is None
    assert rt.log_to_atan(f,g)==2*sx.atan(x**129+1)


def test_native_route_never_calls_the_python_euclidean_recurrence(monkeypatch):
    f,g=pair([1,0,1],[-1,1]);calls=[]
    # A profile checks actual execution without changing stock method identities.
    import sys
    def profile(frame,event,arg):
        if event=='call' and frame.f_code is rt.log_to_atan.__code__:calls.append(1)
    sys.setprofile(profile)
    try:actual=rt.log_to_atan(f,g)
    finally:sys.setprofile(None)
    assert calls==[1] and actual==2*sx.atan(x+1)-2*sx.atan(x+x**2/2+x**3/2)


def test_native_interruptions_and_independence_from_runtime_sympy():
    import os,subprocess,sys
    script="""import sys,signal,random
sys.modules['sympy']=None
import symox as s
from symox.integrals.rationaltools import log_to_atan
x=s.Symbol('x')
f=s.Poly(x,x,domain=s.QQ);g=s.Poly(x+1,x,domain=s.QQ)
assert log_to_atan(f,g)==2*s.atan(2*x+1)+s.pi/2
r=random.Random(7801)
f=s.Poly.from_list([r.randrange(1,8) for _ in range(80)],x,domain=s.QQ)
g=s.Poly.from_list([r.randrange(1,8) for _ in range(75)],x,domain=s.QQ)
signal.signal(signal.SIGALRM,lambda *args:(_ for _ in ()).throw(KeyboardInterrupt()))
signal.setitimer(signal.ITIMER_REAL,0.01)
try:s._core._log_to_atan(f,g)
except KeyboardInterrupt:pass
else:raise AssertionError('native signal was swallowed')
finally:signal.setitimer(signal.ITIMER_REAL,0)
"""
    subprocess.run([sys.executable,'-c',script],env=os.environ.copy(),check=True,timeout=30)


@pytest.mark.parametrize('left,right',[(ZZ,QQ),(QQ,ZZ)])
def test_mixed_rational_ground_domains_keep_the_reference_normalization(left,right):
    f=Poly(x**2+1,x,domain=left);g=Poly(x-1,x,domain=right)
    assert rt.log_to_atan(f,g)==2*sx.atan(x+1)-2*sx.atan(x+x**2/2+x**3/2)
    assert sx._core._log_to_atan(f,g) is not None


def test_complex_variable_values_preserve_principal_atan_sides():
    z=sx.Symbol('z');zp=sp.Symbol('z')
    f,g=Poly(z**2+1,z,domain=QQ),Poly(z-1,z,domain=QQ)
    actual=rt.log_to_atan(f,g)
    expected=ref.log_to_atan(sp.Poly(zp**2+1,zp),sp.Poly(zp-1,zp))
    assert actual==sx.sympify(str(expected),locals={'z':z})
    for a,b in [(sp.Rational(1,1000),2),(-sp.Rational(1,1000),2),(2,-1),(-3,2)]:
        value=sx.Rational(str(a))+sx.I*sx.Rational(str(b))
        left=actual.subs(z,value).evalf(45)
        right=expected.subs(zp,a+sp.I*b).evalf(45)
        assert abs(sp.N(sp.sympify(str(left))-right,35))<sp.Rational(1,10**34)


def test_unused_poly_global_override_does_not_create_an_eligibility_error(monkeypatch):
    f,g=pair([0,1],[1]);calls=[]
    def callback(*args):calls.append(args);raise AssertionError('unused Poly converter')
    monkeypatch.setattr(rt,'Poly',callback)
    assert rt.log_to_atan(f,g)==2*sx.atan(x) and calls==[]


def test_atan_native_factory_delegate_callback_is_preserved(monkeypatch):
    f,g=pair([0,1],[1]);original=sx._core._apply;calls=[]
    def callback(*args,**kwargs):calls.append(args);return original(*args,**kwargs)
    monkeypatch.setattr(sx._core,'_apply',callback)
    assert rt.log_to_atan(f,g)==2*sx.atan(x)
    assert sum(args[0]=='atan' for args in calls)==2  # computation and expected expression
