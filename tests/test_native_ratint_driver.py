"""QQ/ZZ public driver orchestration and certified branch representatives."""
import importlib
import sys

import pytest
import sympy as sp
import symox as sx

rt=importlib.import_module('symox.integrals.rationaltools')
reference=importlib.import_module('sympy.integrals.rationaltools')

CASES=['0','7/3','x**2+3/2','(x*x-1)/(x-1)','1/x','-3/(2*x+2)',
       '1/x**2','3/(x-1)**4','1/(x*x+1)','1/(x*x+x+1)',
       '1/(x*x+2*x+3)','2*x/(x*x-1)','x/(x*x+1)',
       '(x**3+1)/(x*x+1)','(x+1)/(x*x+1)**2',
       '36/(x**5-2*x**4-2*x**3+4*x*x+x-2)',
       '(x*x+1)/(3*x*x-6*x+3)','(x+2)/(2*x*x-2)',
       '((x+1)**4-1)/(x+1)**2','1/(x-1)+1/(x+1)']


def normalize(e):
    return sp.expand(e.replace(lambda a:a.func in (sp.log,sp.atan),
        lambda a:a.func(sp.cancel(a.args[0],extension=True))))


@pytest.mark.parametrize('source',CASES)
@pytest.mark.parametrize('pair',[False,True])
@pytest.mark.parametrize('options',[{}, {'real':True}, {'real':None,'symbol':'residue'}])
def test_complete_native_value_matches_old_driver_reference_and_exact_derivative(source,pair,options,monkeypatch):
    x=sx.Symbol('x');X=sp.Symbol('x');f=sx.sympify(source,locals={'x':x});F=sp.sympify(source,locals={'x':X})
    value=f.as_numer_denom() if pair else f
    expected=reference.ratint(F,X,**options)
    native=sx._core._ratint_driver(value,x,options)
    assert native is not None
    result=rt.ratint(value,x,**options)
    assert result==native[0]
    with monkeypatch.context() as patch:
        patch.setattr(sx._core,'_ratint_driver',lambda *args:None)
        old=rt.ratint(value,x,**options)
    assert normalize(result.to_sympy())==normalize(old.to_sympy())==normalize(expected)
    assert sp.cancel(sp.diff(result.to_sympy(),X)-F,extension=True)==0


@pytest.mark.parametrize('symbol_kind',['string','symbol','dummy'])
@pytest.mark.parametrize('dummy_variable',[False,True])
def test_auxiliary_symbols_and_variable_assumptions_keep_the_same_forms(symbol_kind,dummy_variable,monkeypatch):
    x=sx.Dummy('x',real=True) if dummy_variable else sx.Symbol('x',real=True)
    symbol={'string':'u','symbol':sx.Symbol('u',positive=True),'dummy':sx.Dummy('u',integer=True)}[symbol_kind]
    f=1/(x*x+x+1);flags={'symbol':symbol}
    native=sx._core._ratint_driver(f,x,flags)
    assert native is not None
    result=rt.ratint(f,x,**flags)
    monkeypatch.setattr(sx._core,'_ratint_driver',lambda *args:None)
    assert normalize(result.to_sympy())==normalize(rt.ratint(f,x,**flags).to_sympy())


def test_driver_has_no_python_poly_or_subkernel_recurrence(monkeypatch):
    x=sx.Symbol('x');f=36/(x**5-2*x**4-2*x**3+4*x*x+x-2)
    calls=[];old=sx._core._ratint_driver
    def selected(*args):
        result=old(*args);calls.append(result is not None);return result
    monkeypatch.setattr(sx._core,'_ratint_driver',selected)
    frames=[];codes={fn.__code__ for fn in (rt.ratint_ratpart,rt.ratint_logpart,rt.log_to_real,rt.log_to_atan,
        rt._lrt_constructor.construct,rt._lrt_constructor._options)}
    def profile(frame,event,arg):
        if event=='call' and frame.f_code in codes:frames.append(frame.f_code.co_name)
    previous=sys.getprofile();sys.setprofile(profile)
    try:value=rt.ratint(f,x)
    finally:sys.setprofile(previous)
    assert calls==[True] and frames==[]
    assert sx.cancel(sx.diff(value,x)-f)==0


def test_polynomial_and_rational_only_branches_do_not_require_unobservable_log_protocol(monkeypatch):
    x=sx.Symbol('x');calls=[]
    def observed(*args,**kwargs):calls.append(1);raise KeyboardInterrupt('log protocol')
    monkeypatch.setattr(rt,'log',observed)
    for f in (x*x,sx.S.Zero,1/x**2):
        native=sx._core._ratint_driver(f,x,{})
        assert native is not None and not native[1]
        assert rt.ratint(f,x)==native[0]
    assert calls==[]
    # A genuine log path must restart the old public function and run its hook.
    with pytest.raises(KeyboardInterrupt,match='log protocol'):rt.ratint(1/(x*x+1),x)
    assert calls


@pytest.mark.parametrize('real',[False,0,1,sx.S.true])
def test_unported_real_options_keep_the_complete_public_branch(real,monkeypatch):
    x=sx.Symbol('x');f=1/(x*x+x+1)
    assert sx._core._ratint_driver(f,x,{'real':real}) is None
    actual=rt.ratint(f,x,real=real)
    monkeypatch.setattr(sx._core,'_ratint_driver',lambda *args:None)
    old=rt.ratint(f,x,real=real)
    assert normalize(actual.to_sympy())==normalize(old.to_sympy())


@pytest.mark.parametrize('source',['1/(x**3+x+1)','sqrt(2)/(x*x+1)','1/(x*x+a)','1.0/(x*x+1)'])
def test_uncertified_or_other_domains_restart_the_full_driver(source,monkeypatch):
    x,a=sx.symbols('x a');f=sx.sympify(source,locals={'x':x,'a':a})
    assert sx._core._ratint_driver(f,x,{}) is None
    def run():
        try:return ('value',rt.ratint(f,x))
        except Exception as error:return ('error',type(error),str(error))
    actual=run();monkeypatch.setattr(sx._core,'_ratint_driver',lambda *args:None)
    old=run()
    if actual[0]=='value':assert actual[1].to_sympy().dummy_eq(old[1].to_sympy())
    else:assert actual==old


@pytest.mark.parametrize('f',[1,True,'x',None,(1,0),(0,0),(1,), (1,2,3)])
def test_input_errors_are_identical_to_the_old_driver(f,monkeypatch):
    x=sx.Symbol('x')
    def run():
        try:return ('value',rt.ratint(f,x))
        except Exception as error:return ('error',type(error),str(error))
    actual=run();monkeypatch.setattr(sx._core,'_ratint_driver',lambda *args:None)
    assert run()==actual


@pytest.mark.parametrize('target',['cancel','integrate','div','to_field','primitive','as_expr',
                                 'as_numer_denom','atoms','is_extended_real','as_dummy',
                                 'ratint_ratpart','ratint_logpart','log_to_real','log_to_atan',
                                 'core_ratpart','core_lrt','constructor','constructor_options','iterable','resultant','as_poly'])
def test_real_callbacks_and_interruptions_keep_their_complete_order(target,monkeypatch):
    x=sx.Symbol('x');f=1/(x*x+x+1);flags={'symbol':sx.Symbol('u')};calls=[]
    def observed(*args,**kwargs):calls.append(target);raise KeyboardInterrupt('driver protocol callback')
    if target in ('cancel','integrate','div','to_field','primitive','as_expr'):
        monkeypatch.setattr(sx.Poly,target,observed)
    elif target in ('ratint_ratpart','ratint_logpart','log_to_real','log_to_atan'):
        monkeypatch.setattr(rt,target,observed)
    elif target=='resultant':monkeypatch.setattr(rt,'resultant',observed)
    elif target in ('as_numer_denom','atoms','as_dummy','as_poly'):
        monkeypatch.setattr(sx._core.Expr,target,observed)
    elif target=='is_extended_real':monkeypatch.setattr(type(sx.S.One),target,property(observed))
    elif target=='core_ratpart':monkeypatch.setattr(sx._core,'_ratint_ratpart',observed)
    elif target=='core_lrt':monkeypatch.setattr(sx._core,'_ratint_logpart',observed)
    elif target=='constructor':monkeypatch.setattr(rt._lrt_constructor,'construct',observed)
    elif target=='constructor_options':monkeypatch.setattr(rt._lrt_constructor,'_options',observed)
    else:monkeypatch.setattr(rt._ratint_iterables,'iterable',observed)
    with pytest.raises(KeyboardInterrupt,match='driver protocol callback'):rt.ratint(f,x,**flags)
    trace=list(calls);calls.clear();assert trace
    monkeypatch.setattr(sx._core,'_ratint_driver',lambda *args:None)
    with pytest.raises(KeyboardInterrupt,match='driver protocol callback'):rt.ratint(f,x,**flags)
    assert calls==trace


def test_unrelated_converter_keeps_native_but_applicable_integer_converter_is_observed(monkeypatch):
    x=sx.Symbol('x');calls=[];old=sx._core._ratint_driver
    def selected(*args):
        value=old(*args);calls.append(value is not None);return value
    monkeypatch.setattr(sx._core,'_ratint_driver',selected)
    class Unrelated:pass
    monkeypatch.setitem(rt._ratint_sympify.converter,Unrelated,lambda value:None)
    value=rt.ratint((1,x*x+1),x)
    assert calls==[True]
    events=[]
    def converter(value):events.append(1);raise KeyboardInterrupt('integer conversion')
    monkeypatch.setitem(rt._ratint_sympify.converter,int,converter)
    calls.clear()
    with pytest.raises(KeyboardInterrupt,match='integer conversion'):rt.ratint((1,x*x+1),x)
    assert not calls and events
    # Native Expr inputs bypass the registry in public sympify.
    calls.clear();events.clear()
    assert rt.ratint((sx.S.One,x*x+1),x)==value
    assert calls==[True] and not events


def test_protocol_and_tuple_subclass_eligibility_has_no_side_effects(monkeypatch):
    x=sx.Symbol('x');calls=[]
    class Pair(tuple):
        def __iter__(self):calls.append('iter');raise KeyboardInterrupt('pair iteration')
    pair=Pair((1,x*x+1))
    assert sx._core._ratint_driver(pair,x,{}) is None and not calls
    with pytest.raises(KeyboardInterrupt,match='pair iteration'):rt.ratint(pair,x)
    assert calls==['iter']
    class Term:
        def as_numer_denom(self):calls.append('numerator');return sx.S.One,x*x+1
    term=Term();calls.clear()
    assert sx._core._ratint_driver(term,x,{}) is None and not calls
    value=rt.ratint(term,x,real=True)
    assert calls==['numerator'] and sx.cancel(sx.diff(value,x)-1/(x*x+1))==0


def test_unknown_flags_are_ignored_without_observing_their_values():
    x=sx.Symbol('x')
    class Unknown:
        def __bool__(self):raise AssertionError('unused flag was inspected')
    result=rt.ratint(1/(x*x+1),x,unrecognized=Unknown())
    assert result==sx.atan(x)


def test_sparse_huge_degrees_decline_before_expansion_or_dense_allocation():
    x=sx.Symbol('x')
    assert sx._core._ratint_driver(x**1000000+1,x,{}) is None
    assert sx._core._ratint_driver((x**1000000+1,x+1),x,{}) is None


def test_evaluate_false_keeps_the_full_held_path(monkeypatch):
    x=sx.Symbol('x');f=1/(x*x+1)
    with sx.evaluate(False):
        assert not rt._stock_ratint_protocol(f)
        def run():
            try:return ('value',sx.srepr(rt.ratint(f,x)))
            except Exception as error:return ('error',type(error),str(error))
        result=run();monkeypatch.setattr(sx._core,'_ratint_driver',lambda *args:None)
        assert run()==result


def test_noref_native_driver_with_vector_ore_and_reload():
    import os,subprocess
    script="""import sys,importlib
sys.modules['sympy']=None
import symox as S
import symox.vector,symox.holonomic.recurrence
r=importlib.import_module('symox.integrals.rationaltools')
x=S.Symbol('x');f=1/(x*x+x+1)
for i in range(3):
 importlib.reload(r)
 assert r._stock_ratint_protocol(f) and r._stock_ratint_log_protocol()
 native=S._core._ratint_driver(f,x,{})
 assert native is not None and native[1]
 result=r.ratint(f,x)
 assert result==native[0]
 assert S.cancel(S.diff(result,x)-f)==0
assert sys.modules['sympy'] is None
"""
    subprocess.run([sys.executable,'-c',script],env=os.environ.copy(),check=True)


@pytest.mark.parametrize('source',['-3/(2*x+2)','2*x/(x*x-1)','1/(x*x+x+1)',
    '36/(x**5-2*x**4-2*x**3+4*x*x+x-2)'])
def test_poles_and_complex_specializations_keep_old_principal_branches(source,monkeypatch):
    x=sx.Symbol('x');f=sx.sympify(source,locals={'x':x})
    native=sx._core._ratint_driver(f,x,{})
    assert native is not None
    result=rt.ratint(f,x)
    monkeypatch.setattr(sx._core,'_ratint_driver',lambda *args:None)
    old=rt.ratint(f,x)
    for point in (-3,-1,0,1,2,3,sx.I,1+sx.I,-1-sx.I):
        assert normalize(result.subs(x,point).to_sympy())==normalize(old.subs(x,point).to_sympy())


@pytest.mark.parametrize('shape',['degree','numeric'])
def test_compounded_power_refusal_preserves_the_complete_public_result(shape,monkeypatch):
    x=sx.Symbol('x')
    if shape=='degree':f=sx.Pow(sx.Pow(x+1,8,evaluate=False)+1,8,evaluate=False)
    else:f=sx.Pow(sx.Integer(2)**4095,2,evaluate=False)*x
    assert sx._core._ratint_driver(f,x,{}) is None
    actual=rt.ratint(f,x)
    monkeypatch.setattr(sx._core,'_ratint_driver',lambda *args:None)
    old=rt.ratint(f,x)
    assert actual==old
    assert sx.expand(sx.diff(actual,x)-f)==0


def test_changed_lrt_guard_is_observed_once_in_the_complete_path(monkeypatch):
    x=sx.Symbol('x');f=1/(x*x+1);calls=[]
    def observed():calls.append(1);return False
    monkeypatch.setattr(rt,'_stock_lrt_resultant',observed)
    actual=rt.ratint(f,x)
    first=list(calls);calls.clear()
    monkeypatch.setattr(sx._core,'_ratint_driver',lambda *args:None)
    assert rt.ratint(f,x)==actual and calls==first and calls


def test_changed_registry_zero_keeps_the_old_log_accumulation(monkeypatch):
    x=sx.Symbol('x');f=1/(x*x+1)
    monkeypatch.setattr(rt.S,'Zero',sx.Integer(7))
    actual=rt.ratint(f,x)
    monkeypatch.setattr(sx._core,'_ratint_driver',lambda *args:None)
    assert rt.ratint(f,x)==actual
