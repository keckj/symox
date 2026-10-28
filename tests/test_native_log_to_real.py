"""Exact linear/quadratic residues, including their log/atan representatives."""
import importlib
import sys

import pytest
import sympy as sp
import symox as sx

rt=importlib.import_module('symox.integrals.rationaltools')
ref=importlib.import_module('sympy.integrals.rationaltools')

CASES=[(h,q) for h in ('x+t','x**2+t*x+1','x**3+(t**2+1)*x+t',
                        '(t+1)*x**2+(t**2-2)*x+2*t')
       for q in ('t-2','2*t+3','t','t**2+1','3*t**2+1',
                 't**2+2*t+3','-2*t**2+4*t-7','t**2-2',
                 't**2+3*t+1','2*t**2-5*t+2')]


def pair(h='x+t',q='3*t**2+1',domain='QQ'):
    x,t=sx.symbols('x t')
    return x,t,sx.Poly(sx.sympify(h,locals={'x':x,'t':t}),x,domain=domain+'[t]'),sx.Poly(sx.sympify(q,locals={'t':t}),t,domain=domain)


def normalize_representative(e):
    # Normalize polynomial arguments and algebraic coefficients only. Never
    # combine logarithms or apply an atan addition/reciprocal identity.
    return sp.expand(e.replace(lambda f:f.func in (sp.log,sp.atan),
        lambda f:f.func(sp.cancel(f.args[0],extension=True))))


@pytest.mark.parametrize('h,q',CASES)
@pytest.mark.parametrize('domain',['ZZ','QQ'])
def test_exact_representatives_match_reference_and_complete_old_path(h,q,domain,monkeypatch):
    x,t,H,Q=pair(h,q,domain);X,T=sp.symbols('x t')
    expected=ref.log_to_real(sp.Poly(h,X,domain=domain+'[t]'),sp.Poly(q,T,domain=domain),X,T)
    direct=sx._core._log_to_real(H,Q,x,t)
    assert direct is not None
    assert rt.log_to_real(H,Q,x,t)==direct
    with monkeypatch.context() as patch:
        patch.setattr(sx._core,'_log_to_real',lambda *args:None)
        old=rt.log_to_real(H,Q,x,t)
    assert normalize_representative(direct.to_sympy())==normalize_representative(expected)
    assert normalize_representative(old.to_sympy())==normalize_representative(expected)
    # Independent rational trace certificate, without computing any roots.
    K=sp.QQ.frac_field(X);qpoly=sp.Poly(q,T,domain=K);hpoly=sp.Poly(h,T,domain=K)
    r=(sp.Poly(T*sp.diff(sp.sympify(h),X),T,domain=K)*sp.invert(hpoly,qpoly)).rem(qpoly)
    trace=sum((r*sp.Poly(T**j,T,domain=K)).rem(qpoly).nth(j) for j in range(qpoly.degree()))
    assert sp.cancel(sp.diff(direct.to_sympy(),X)-trace,extension=True)==0


@pytest.mark.parametrize('q',['t-2','t**2-2','3*t**2+1','t**2+2*t+3'])
def test_real_and_complex_substitutions_preserve_the_same_branches(q,monkeypatch):
    x,t,h,q=pair(q=q);value=rt.log_to_real(h,q,x,t)
    with monkeypatch.context() as patch:
        patch.setattr(sx._core,'_log_to_real',lambda *args:None)
        old=rt.log_to_real(h,q,x,t)
    for point in (-3,-1,0,1,3,sx.I,1+sx.I,-1-sx.I,sx.sqrt(2)):
        a,b=value.subs(x,point),old.subs(x,point)
        assert normalize_representative(a.to_sympy())==normalize_representative(b.to_sympy())
        assert a.doit()==a and b.doit()==b


@pytest.mark.parametrize('source',['1/(x*x+1)','1/(x*x+x+1)','1/(x*x+2*x+3)'])
@pytest.mark.parametrize('options',[{}, {'real':True}, {'real':True,'symbol':'residue'}])
def test_public_ratint_uses_native_conversion_without_python_root_enumeration(source,options,monkeypatch):
    x=sx.Symbol('x');f=sx.sympify(source,locals={'x':x});calls=[];old=sx._core._log_to_real
    def selected(*args):
        value=old(*args);calls.append(value is not None);return value
    monkeypatch.setattr(sx._core,'_log_to_real',selected)
    host=[]
    def profile(frame,event,arg):
        if event=='call' and frame.f_code is rt._get_real_roots.__code__:host.append(1)
    previous=sys.getprofile();sys.setprofile(profile)
    try:value=rt.ratint(f,x,**options)
    finally:sys.setprofile(previous)
    assert calls and all(calls) and host==[]
    assert sx.cancel(sx.diff(value,x)-f)==0


def test_real_false_option_keeps_the_rootsum_path(monkeypatch):
    x=sx.Symbol('x')
    def forbidden(*args):raise AssertionError('real=False entered real conversion')
    monkeypatch.setattr(sx._core,'_log_to_real',forbidden)
    result=rt.ratint(1/(x*x+x+1),x,real=False)
    assert sx.cancel(sx.diff(result,x)-1/(x*x+x+1))==0


@pytest.mark.parametrize('hs,qs,domain',[
    ('x+t','(t+1)**2','QQ'),('x+t','t**3+t+1','QQ'),
    ('x+t','0','QQ'),('x+t','1','QQ'),('x+t','t*t+1','EX'),
    ('x+t','t*t+1','RR'),('0','t*t+1','QQ'),
    ('t*(x+1)','t*t+1','QQ'),('x**33+t','t+1','QQ')])
def test_declined_domains_and_degenerate_cases_keep_complete_errors_or_results(hs,qs,domain,monkeypatch):
    x,t=sx.symbols('x t');h=sx.Poly(sx.sympify(hs,locals={'x':x,'t':t}),x,domain=domain+'[t]' if domain in ('QQ','RR') else domain)
    q=sx.Poly(sx.sympify(qs,locals={'t':t}),t,domain=domain)
    assert sx._core._log_to_real(h,q,x,t) is None
    def outcome():
        try:return ('value',rt.log_to_real(h,q,x,t))
        except Exception as error:return ('error',type(error),str(error))
    actual=outcome();monkeypatch.setattr(sx._core,'_log_to_real',lambda *args:None)
    assert outcome()==actual


@pytest.mark.parametrize('target',['roots','root_helper','symbols','dummy','collect','log',
                                    'as_expr','count_roots','xreplace','expand','evalf','subs','priority',
                                    'xreplace_helper','accum_helper','dummy_new','core_dummy','args'])
def test_callbacks_and_interruptions_remain_observable_in_the_complete_path(target,monkeypatch):
    x,t,h,q=pair(q='t**2-2' if target in ('subs','accum_helper') else 't**2+2*t+3');calls=[]
    def observed(*args,**kwargs):calls.append(1);raise KeyboardInterrupt('real conversion callback')
    if target in ('roots','symbols','log'):
        monkeypatch.setattr(rt,target,observed)
    elif target=='dummy':monkeypatch.setattr(rt,'Dummy',observed)
    elif target=='root_helper':monkeypatch.setattr(rt,'_get_real_roots',observed)
    elif target=='collect':monkeypatch.setattr(rt._real_radsimp_module,'collect',observed)
    elif target=='xreplace_helper':monkeypatch.setattr(rt._real_subs_module,'xreplace_rule',observed)
    elif target=='accum_helper':
        method=sx._core.Expr.subs
        cells=dict(zip(method.__code__.co_freevars,method.__closure__))
        monkeypatch.setattr(cells['_CREATED'],'cell_contents',[False])
        monkeypatch.setattr(rt._real_accum_module,'_substitutes_illegal',observed)
    elif target=='dummy_new':monkeypatch.setattr(rt.Dummy,'__new__',staticmethod(observed))
    elif target=='core_dummy':monkeypatch.setattr(sx._core,'Dummy',observed)
    elif target=='args':monkeypatch.setattr(sx._core.Expr,'args',property(observed))
    elif target in ('as_expr','count_roots'):monkeypatch.setattr(sx.Poly,target,observed)
    elif target=='priority':monkeypatch.setattr(sx._core.Expr,'_op_priority',property(observed))
    else:monkeypatch.setattr(sx._core.Expr,target,observed)
    original=sx._core._log_to_real
    def selected(*args):
        value=original(*args)
        assert value is None, 'callback bypassed by native conversion'
        return value
    monkeypatch.setattr(sx._core,'_log_to_real',selected)
    with pytest.raises(KeyboardInterrupt,match='real conversion callback'):rt.log_to_real(h,q,x,t)
    trace=list(calls);calls.clear()
    monkeypatch.setattr(sx._core,'_log_to_real',lambda *args:None)
    with pytest.raises(KeyboardInterrupt,match='real conversion callback'):rt.log_to_real(h,q,x,t)
    assert trace==calls and calls


def test_evaluate_false_retains_held_constructors(monkeypatch):
    x,t,h,q=pair()
    def outcome():
        try:return ('value',sx.srepr(rt.log_to_real(h,q,x,t)))
        except Exception as error:return ('error',type(error),str(error))
    with sx.evaluate(False):
        assert not rt._stock_real_protocol(h,q,rt._real_radsimp_module.collect)
        actual=outcome()
        monkeypatch.setattr(sx._core,'_log_to_real',lambda *args:None)
        assert outcome()==actual


def test_import_orders_reload_and_noref_runtime():
    import os,subprocess
    script="""import sys,importlib
sys.modules['sympy']=None
import symox as S
import symox.vector,symox.holonomic.recurrence
r=importlib.import_module('symox.integrals.rationaltools')
x,t=S.symbols('x t');h=S.Poly(x+3*t/2+S.S.Half,x,domain='QQ[t]');q=S.Poly(3*t*t+1,t,domain='ZZ')
for i in range(3):
 importlib.reload(r)
 assert r._stock_real_protocol(h,q,r._real_radsimp_module.collect)
 result=r.log_to_real(h,q,x,t)
 assert result==2*S.sqrt(3)*S.atan(2*S.sqrt(3)*x/3+S.sqrt(3)/3)/3
assert sys.modules['sympy'] is None
"""
    subprocess.run([sys.executable,'-c',script],env=os.environ.copy(),check=True)


@pytest.mark.parametrize('hs,qs',[
    ('x**2/7+(t/3+2/5)*x+t**2/11','2*t**2/3+5*t/7+11/13'),
    ('x+t/7','-t**2+2*t-17'),
    ('x+t+2**80/3','3*t**2+1'),
    ('x+t/3','t/7+2/3')])
def test_rational_and_large_integer_coefficients_are_exact(hs,qs):
    x,t,h,q=pair(hs,qs);X,T=sp.symbols('x t')
    actual=sx._core._log_to_real(h,q,x,t)
    assert actual is not None
    expected=ref.log_to_real(sp.Poly(hs,X,domain='QQ[t]'),sp.Poly(qs,T,domain='QQ'),X,T)
    assert normalize_representative(actual.to_sympy())==normalize_representative(expected)


@pytest.mark.parametrize('target',['as_expr_code','symbols_delegate','xreplace_code'])
def test_inplace_code_and_forwarder_global_changes_are_not_bypassed(target,monkeypatch):
    x,t,h,q=pair();calls=[]
    def observed(*args,**kwargs):calls.append(1);raise KeyboardInterrupt('modified real protocol')
    if target=='symbols_delegate':monkeypatch.setitem(rt.symbols.__globals__,'_symbols',observed)
    else:
        function=sx.Poly.as_expr if target=='as_expr_code' else sx._core.Expr.xreplace
        names=','.join('v%d'%i for i in range(len(function.__closure__ or ())))
        namespace={'observed':observed}
        use='\n        if ('+names+',):pass' if names else ''
        exec('def factory('+names+'):\n    def changed(*args,**kwargs):'+use+'\n        raise KeyboardInterrupt("modified real protocol")\n    return changed',namespace)
        monkeypatch.setattr(function,'__code__',namespace['factory'](*([None]*len(function.__closure__ or ()))).__code__)
    original=sx._core._log_to_real
    def selected(*args):
        value=original(*args);assert value is None;return value
    monkeypatch.setattr(sx._core,'_log_to_real',selected)
    with pytest.raises(KeyboardInterrupt,match='modified real protocol'):rt.log_to_real(h,q,x,t)
    monkeypatch.setattr(sx._core,'_log_to_real',lambda *args:None)
    with pytest.raises(KeyboardInterrupt,match='modified real protocol'):rt.log_to_real(h,q,x,t)


def test_custom_inputs_and_multivariate_or_fraction_domains_keep_the_full_path(monkeypatch):
    x,t,h,q=pair()
    class Protocol:
        def as_expr(self):raise KeyboardInterrupt('custom polynomial protocol')
    obj=Protocol()
    assert sx._core._log_to_real(obj,q,x,t) is None
    with pytest.raises(KeyboardInterrupt,match='custom polynomial protocol'):rt.log_to_real(obj,q,x,t)
    for h in (sx.Poly(x+t,x,t,domain='QQ'),sx.Poly(x+t,x,domain='QQ(t)')):
        assert sx._core._log_to_real(h,q,x,t) is None
        result=rt.log_to_real(h,q,x,t)
        with monkeypatch.context() as patch:
            patch.setattr(sx._core,'_log_to_real',lambda *args:None)
            assert rt.log_to_real(h,q,x,t)==result


def test_sparse_huge_degree_is_declined_before_any_dense_allocation():
    x,t,h,q=pair()
    q=sx.Poly(t**1000000+1,t,domain='QQ')
    assert sx._core._log_to_real(h,q,x,t) is None
