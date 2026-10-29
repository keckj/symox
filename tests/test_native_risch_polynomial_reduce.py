"""Bounded exact Risch polynomial reduction; unchanged complete old fallback."""
from fractions import Fraction
import importlib
import random

import pytest
import symox as sx
import sympy as sp

risch=importlib.import_module('symox.integrals.risch')
protocol=importlib.import_module('symox.integrals._risch_native')
reference=importlib.import_module('sympy.integrals.risch')


def input_case(seed,package,risch_module):
    rng=random.Random(884311+seed)
    x,t,u=package.symbols('risch_x risch_t risch_u')
    variables=[x,t,u][:1+seed%3]
    t=variables[-1]
    choices=[package.S.One,t,t*t+1,3*t**3+2*t-1,(t*t+1)/2]
    dt=choices[seed%len(choices)]
    coefficients=[package.Rational(rng.randrange(-19,20),1 if seed%2==0 else rng.randrange(1,11))
                  for _ in range(seed%21)]
    p=package.Poly(sum(c*t**i for i,c in enumerate(coefficients)),t,domain='ZZ' if seed%2==0 else 'QQ')
    ds=[package.Poly(1,v)for v in variables[:-1]]+[package.Poly(dt,t)]
    de=risch_module.DifferentialExtension(extension={'D':ds})
    return p,de


def clear_copy_slot_caches(monkeypatch):
    """Undo copyreg's ``__slotnames__`` cache left by an earlier copy/pickle.

    ``copy.copy`` of a DifferentialExtension anywhere earlier in the process
    stores ``__slotnames__`` in the class dictionary; the exact class proof of
    the admission then (correctly) declines. Remove only that cache, restored
    at teardown, so the stock assertions below do not depend on test order.
    """
    from symox import _node_types as nodes
    for cls in (risch.DifferentialExtension,sx._core.Poly,nodes.node_class('NegativeInfinity')):
        if '__slotnames__' in vars(cls):monkeypatch.delattr(cls,'__slotnames__')


def signature(pair):
    return tuple((p.gens[0].name,str(p.domain),tuple(sorted(
        (monomial, Fraction(int(value.p),int(value.q)))for monomial,value in p.as_dict().items())))for p in pair)


@pytest.mark.parametrize('seed',range(120))
def test_selected_complete_loop_matches_old_and_reference_with_exact_certificate(monkeypatch,seed):
    clear_copy_slot_caches(monkeypatch)
    p,de=input_case(seed,sx,risch);rp,rde=input_case(seed,sp,reference)
    original=(p,de.level,de.t,de.d,tuple(de.T),tuple(de.D));calls=[]
    native=sx._core._risch_polynomial_reduce
    def spy(*args):
        out=native(*args);calls.append(out is not None);return out
    monkeypatch.setattr(sx._core,'_risch_polynomial_reduce',spy)
    actual=risch.polynomial_reduce(p,de)
    no_step=p.is_zero or p.degree(de.t)<de.d.degree(de.t)
    assert calls==([] if no_step else [True])
    with monkeypatch.context() as patch:
        patch.setattr(risch,'_native_polynomial_reduce',lambda *args:None)
        old=risch.polynomial_reduce(p,de)
    expected=reference.polynomial_reduce(rp,rde)
    assert signature(actual)==signature(old)==signature(expected)
    q,r=actual
    assert p.as_dict()==(risch.derivation(q,de)+r).as_dict()
    assert r.is_zero or r.degree(de.t)<de.d.degree(de.t)
    assert original==(p,de.level,de.t,de.d,tuple(de.T),tuple(de.D))
    assert all(a is b for a,b in zip(original[-1],de.D))
    if p.is_zero or p.degree(de.t)<de.d.degree(de.t):assert r is p


def old_and_new(monkeypatch,p,de):
    calls=[];native=sx._core._risch_polynomial_reduce
    def spy(*args):
        out=native(*args);calls.append(out is not None);return out
    monkeypatch.setattr(sx._core,'_risch_polynomial_reduce',spy)
    actual=risch.polynomial_reduce(p,de)
    with monkeypatch.context() as patch:
        patch.setattr(risch,'_native_polynomial_reduce',lambda *args:None)
        expected=risch.polynomial_reduce(p,de)
    assert actual==expected
    assert all(a.gens==b.gens and str(a.domain)==str(b.domain) and a.as_dict()==b.as_dict()
               for a,b in zip(actual,expected))
    return actual,calls


@pytest.mark.parametrize('domain',['RR','EX','GF(7)','QQ[y]'])
def test_unsupported_domains_keep_the_whole_original_algorithm(monkeypatch,domain):
    t,y=sx.symbols('risch_t y')
    if domain=='QQ[y]':
        p=sx.Poly(t**4+y*t+1,t,domain=domain)
    elif domain=='RR':
        p=sx.Poly(t**4+sx.Float('1.25')*t+1,t,domain=domain)
    else:p=sx.Poly(t**4+2*t+1,t,domain=domain)
    de=risch.DifferentialExtension(extension={'D':[sx.Poly(t*t+1,t)]})
    # Domain signatures above intentionally use exact p/q coefficients only.
    with monkeypatch.context() as patch:
        patch.setattr(risch,'_native_polynomial_reduce',lambda *args:None);expected=risch.polynomial_reduce(p,de)
    calls=[];native=sx._core._risch_polynomial_reduce
    def spy(*args):
        out=native(*args);calls.append(out is not None);return out
    monkeypatch.setattr(sx._core,'_risch_polynomial_reduce',spy)
    actual=risch.polynomial_reduce(p,de)
    assert actual==expected and not any(calls)


@pytest.mark.parametrize('mode',['derivation_alias','derivation_code','degree','LC','mul','as_poly','as_expr','module_Poly','module_len','module_zip'])
def test_consumed_callbacks_and_changed_function_bodies_are_preserved(monkeypatch,mode):
    p,de=input_case(19,sx,risch);events=[]
    if mode=='derivation_alias':
        original=risch.derivation
        def changed(*args,**kwargs):events.append('derivation');return original(*args,**kwargs)
        monkeypatch.setattr(risch,'derivation',changed)
    elif mode=='derivation_code':
        from types import FunctionType
        original=risch.derivation
        delegate=FunctionType(original.__code__,original.__globals__,argdefs=original.__defaults__)
        monkeypatch.setattr(risch,'_native_test_derivation_events',events,raising=False)
        monkeypatch.setattr(risch,'_native_test_derivation_delegate',delegate,raising=False)
        def changed(p,DE,coefficientD=False,basic=False):
            _native_test_derivation_events.append('changed-code')
            return _native_test_derivation_delegate(p,DE,coefficientD,basic)
        monkeypatch.setattr(original,'__code__',changed.__code__)
    elif mode=='constructor':
        original=sx.Poly.__new__
        def changed(cls,*args,**kwargs):events.append('constructor');return original(cls,*args,**kwargs)
        monkeypatch.setattr(sx.Poly,'__new__',staticmethod(changed))
    elif mode=='module_Poly':
        original=risch.Poly
        def changed(*args,**kwargs):events.append('Poly');return original(*args,**kwargs)
        monkeypatch.setattr(risch,'Poly',changed)
    elif mode in ('module_len','module_zip'):
        name=mode[len('module_'):];original={'len':len,'zip':zip}[name]
        def changed(*args):events.append(name);return original(*args)
        monkeypatch.setattr(risch,name,changed,raising=False)
    else:
        original=getattr(sx.Poly,mode)
        def changed(self,*args,**kwargs):events.append(mode);return original(self,*args,**kwargs)
        monkeypatch.setattr(sx.Poly,mode,changed)
    if mode=='constructor':
        with pytest.raises(TypeError) as actual:risch.polynomial_reduce(p,de)
        observed=events[:];events.clear()
        with monkeypatch.context() as patch:
            patch.setattr(risch,'_native_polynomial_reduce',lambda *args:None)
            with pytest.raises(TypeError) as expected:risch.polynomial_reduce(p,de)
        assert str(actual.value)==str(expected.value) and observed==events and events
        return
    native=sx._core._risch_polynomial_reduce;calls=[]
    def spy(*args):
        out=native(*args);calls.append(out is not None);return out
    monkeypatch.setattr(sx._core,'_risch_polynomial_reduce',spy)
    actual=risch.polynomial_reduce(p,de);observed=events[:];events.clear()
    with monkeypatch.context() as patch:
        patch.setattr(risch,'_native_polynomial_reduce',lambda *args:None);expected=risch.polynomial_reduce(p,de)
    assert actual==expected and events==observed and not any(calls)
    assert events


@pytest.mark.parametrize('mode',['subclass','class_lookup','slot_property'])
def test_extension_attribute_protocols_keep_real_observable_reads(monkeypatch,mode):
    p,de=input_case(19,sx,risch);events=[]
    if mode=='subclass':
        class Custom(risch.DifferentialExtension):
            def __getattribute__(self,name):events.append(name);return super().__getattribute__(name)
        de=Custom(extension={'D':de.D});events.clear()
    elif mode=='class_lookup':
        def changed(self,name):events.append(name);return object.__getattribute__(self,name)
        monkeypatch.setattr(risch.DifferentialExtension,'__getattribute__',changed)
    else:
        descriptor=risch.DifferentialExtension.t
        def changed(self):events.append('t');return descriptor.__get__(self,type(self))
        monkeypatch.setattr(risch.DifferentialExtension,'t',property(changed))
    actual=risch.polynomial_reduce(p,de);observed=events[:];events.clear()
    with monkeypatch.context() as patch:
        patch.setattr(risch,'_native_polynomial_reduce',lambda *args:None);expected=risch.polynomial_reduce(p,de)
    assert actual==expected and observed==events and observed


def test_uninitialized_unused_extension_slots_do_not_create_a_new_error(monkeypatch):
    t=sx.Symbol('risch_t');de=risch.DifferentialExtension.__new__(risch.DifferentialExtension)
    de.t=t;de.d=sx.Poly(t*t+1,t);p=sx.Poly(t,t)
    actual,calls=old_and_new(monkeypatch,p,de)
    assert actual[1] is p and not any(calls)


@pytest.mark.parametrize('mode',['degree','bits','fraction_field','future_variables','duplicate_generators'])
def test_nonadmitted_inputs_keep_available_original_computations(monkeypatch,mode):
    t,x=sx.symbols('risch_t risch_x')
    p=sx.Poly(t**5+2*t+1,t);de=risch.DifferentialExtension(extension={'D':[sx.Poly(t*t+1,t)]})
    if mode=='degree':p=sx.Poly(t**257+1,t)
    elif mode=='bits':p=sx.Poly((2**4100)*t**3+1,t)
    elif mode=='fraction_field':
        de=risch.DifferentialExtension(extension={'D':[sx.Poly(1,x),sx.Poly(t*t/x+1,t)]})
    elif mode=='future_variables':
        de=risch.DifferentialExtension(extension={'D':[sx.Poly(t*t+1,t),sx.Poly(x,x)]});de.decrement_level()
    else:
        de.T=[t,t];de.D=[sx.Poly(1,t),sx.Poly(t*t+1,t)];de.level=-1
    actual,calls=old_and_new(monkeypatch,p,de)
    if mode!='future_variables':assert not any(calls)
    assert p.as_dict()==(risch.derivation(actual[0],de)+actual[1]).as_dict()


def test_late_import_hook_sees_only_the_existing_algorithm_imports():
    import builtins
    p,de=input_case(19,sx,risch);events=[]
    original=builtins.__import__;native=risch._native_polynomial_reduce
    def changed(name,globals=None,locals=None,fromlist=(),level=0):
        events.append((name,fromlist,level))
        return original(name,globals,locals,fromlist,level)
    try:
        builtins.__import__=changed
        actual=risch.polynomial_reduce(p,de);observed=events[:];events.clear()
        risch._native_polynomial_reduce=lambda *args:None
        expected=risch.polynomial_reduce(p,de)
    finally:
        builtins.__import__=original
        risch._native_polynomial_reduce=native
    assert actual==expected and events==observed and events


def test_stock_native_admission_is_restored_after_protocol_mutations(monkeypatch):
    clear_copy_slot_caches(monkeypatch)
    p,de=input_case(18,sx,risch)
    calls=[];native=sx._core._risch_polynomial_reduce
    def spy(*args):
        out=native(*args);calls.append(out is not None);return out
    monkeypatch.setattr(sx._core,'_risch_polynomial_reduce',spy)
    q,r=risch.polynomial_reduce(p,de)
    assert calls==[True]
    assert p.as_dict()==(risch.derivation(q,de)+r).as_dict()


def test_stock_assertions_are_independent_of_an_earlier_extension_copy(monkeypatch):
    import copy
    p,de=input_case(18,sx,risch)
    had='__slotnames__' in vars(risch.DifferentialExtension)
    copy.copy(de)
    assert '__slotnames__' in vars(risch.DifferentialExtension)
    calls=[];native=sx._core._risch_polynomial_reduce
    def spy(*args):
        out=native(*args);calls.append(out is not None);return out
    monkeypatch.setattr(sx._core,'_risch_polynomial_reduce',spy)
    with monkeypatch.context() as patch:
        clear_copy_slot_caches(patch)
        q,r=risch.polynomial_reduce(p,de)
        assert calls==[True]
    assert p.as_dict()==(risch.derivation(q,de)+r).as_dict()
    if not had:
        # Leave the process as an untouched import would have it.
        del risch.DifferentialExtension.__slotnames__


def test_modified_constructor_error_is_preserved_in_an_isolated_type_process():
    import os,subprocess,sys
    source="""
import symox as s
from symox.integrals import risch
x=s.Symbol('risch_t');p=s.Poly(x**4+1,x)
de=risch.DifferentialExtension(extension={'D':[s.Poly(x*x+1,x)]})
original=s.Poly.__new__;events=[]
def changed(cls,*args,**kwargs):
 events.append('constructor');return original(cls,*args,**kwargs)
s.Poly.__new__=staticmethod(changed)
errors=[]
for mode in range(2):
 if mode:risch._native_polynomial_reduce=lambda *args:None
 events.clear()
 try:risch.polynomial_reduce(p,de)
 except TypeError as error:errors.append(str(error))
 else:raise AssertionError('original constructor error omitted')
 assert events==['constructor']
assert errors[0]==errors[1]
"""
    subprocess.run([sys.executable,'-c',source],check=True,env=os.environ.copy())


def test_restored_constructor_dictionary_does_not_hide_a_changed_native_new_slot():
    import os,subprocess,sys
    source="""
import symox as s
from symox.integrals import risch
x=s.Symbol('risch_t');p=s.Poly(x**4+1,x);zero=s.Poly(0,x)
de=risch.DifferentialExtension(extension={'D':[s.Poly(x*x+1,x)]})
original=s.Poly.__new__
def changed(cls,*args,**kwargs):return original(cls,*args,**kwargs)
s.Poly.__new__=staticmethod(changed)
s.Poly.__new__=original
calls=[];native=s._core._risch_polynomial_reduce
native_original=risch._native_polynomial_reduce
def spy(*args):
 out=native(*args);calls.append(out is not None);return out
s._core._risch_polynomial_reduce=spy
def prefix(*args,**kwargs):
 risch.Poly=s.Poly
 return zero
errors=[]
for mode in range(2):
 risch.Poly=prefix
 if mode:risch._native_polynomial_reduce=lambda *args:None
 try:risch.polynomial_reduce(p,de)
 except TypeError as error:errors.append(str(error))
 else:raise AssertionError('restored dictionary hid unsafe constructor slot')
assert errors[0]==errors[1] and not any(calls)
"""
    subprocess.run([sys.executable,'-c',source],check=True,env=os.environ.copy())


@pytest.mark.parametrize('scenario',['stock','foreign_globals','defaults','kwdefaults','closure','nested_constants','guard_callback'])
def test_rationaltools_reload_authenticates_code_namespace_defaults_and_closure(monkeypatch,scenario):
    from types import FunctionType,CodeType
    rt=importlib.import_module('symox.integrals.rationaltools')
    importlib.reload(rt)
    clear_copy_slot_caches(monkeypatch)
    p,de=input_case(19,sx,risch);events=[]
    current=rt._stock_ratint_protocol
    if scenario=='foreign_globals':
        monkeypatch.setattr(rt,'_stock_ratint_protocol',FunctionType(current.__code__,dict(rt.__dict__)))
    elif scenario=='defaults':
        replacement=FunctionType(current.__code__,rt.__dict__,argdefs=(None,))
        monkeypatch.setattr(rt,'_stock_ratint_protocol',replacement)
    elif scenario=='kwdefaults':
        replacement=FunctionType(current.__code__,rt.__dict__);replacement.__kwdefaults__={'foreign':None}
        monkeypatch.setattr(rt,'_stock_ratint_protocol',replacement)
    elif scenario=='closure':
        def cell(value):return (lambda:value).__closure__[0]
        code=current.__code__.replace(co_freevars=('foreign_cell',))
        replacement=FunctionType(code,rt.__dict__,closure=(cell(events),))
        monkeypatch.setattr(rt,'_stock_ratint_protocol',replacement)
    elif scenario=='nested_constants':
        constants=list(current.__code__.co_consts)
        index=next(i for i,c in enumerate(constants)if type(c)is CodeType)
        constants[index]=constants[index].replace(co_consts=constants[index].co_consts+('foreign nested constant',))
        replacement=FunctionType(current.__code__.replace(co_consts=tuple(constants)),rt.__dict__)
        monkeypatch.setattr(rt,'_stock_ratint_protocol',replacement)
    elif scenario=='guard_callback':
        def foreign(*args,**kwargs):events.append('foreign guard');return True
        monkeypatch.setattr(rt,'_atan_stock_operators',foreign)
    calls=[];native=sx._core._risch_polynomial_reduce
    def spy(*args):
        value=native(*args);calls.append(value is not None);return value
    monkeypatch.setattr(sx._core,'_risch_polynomial_reduce',spy)
    actual=risch.polynomial_reduce(p,de)
    with monkeypatch.context() as patch:
        patch.setattr(risch,'_native_polynomial_reduce',lambda *args:None)
        old=risch.polynomial_reduce(p,de)
    assert actual==old and events==[]
    if scenario=='stock':assert calls==[True]
    else:assert not any(calls)


@pytest.mark.parametrize('scenario',['before_import','alias_after_reload','alias_raising_after_reload','module_slot_after_reload','alias_foreign_owner_after_reload','ground_alias_after_reload','ground_alias_raising_after_reload'])
@pytest.mark.parametrize('input_kind',['valid','invalid'])
def test_serializer_callbacks_are_never_captured_or_called(scenario,input_kind):
    import os,subprocess,sys
    source="""
import marshal,sys
scenario,input_kind=sys.argv[1:]
original=marshal.dumps;events=[]
def wrapper(*args,**kwargs):
 events.append('serializer')
 if scenario in ('alias_raising_after_reload','ground_alias_raising_after_reload'):raise RuntimeError('foreign serializer')
 return original(*args,**kwargs)
if scenario=='before_import':marshal.dumps=wrapper
import symox as s
from symox.integrals import risch,rationaltools
from symox.integrals import _risch_native as selector
import importlib
x=s.Symbol('risch_t');p=s.Poly(x**4+2*x+1,x)
de=risch.DifferentialExtension(extension={'D':[s.Poly(x*x+1,x)]})
if input_kind=='invalid':p=None
importlib.reload(rationaltools)
if scenario in ('alias_after_reload','alias_raising_after_reload'):selector._CODE=wrapper
elif scenario=='module_slot_after_reload':marshal.dumps=wrapper
elif scenario in ('ground_alias_after_reload','ground_alias_raising_after_reload'):
 from symox.polys import _polymatrix_native as ground
 ground._CODE=wrapper
elif scenario=='alias_foreign_owner_after_reload':
 import _pickle
 selector._CODE=_pickle.dumps
 def trace(frame,event,arg):
  if event=='c_call' and arg is _pickle.dumps:events.append('foreign builtin')
 sys.setprofile(trace)
calls=[];native=s._core._risch_polynomial_reduce
s._core._risch_polynomial_reduce=lambda *args:(lambda value:(calls.append(value is not None),value)[1])(native(*args))
def outcome():
 try:
  pair=risch.polynomial_reduce(p,de)
  return ('value',[(part.as_dict(),str(part.domain),part.gens)for part in pair])
 except Exception as error:return ('error',type(error).__name__,str(error))
events.clear();actual=outcome();observed=events[:]
events.clear();risch._native_polynomial_reduce=lambda *args:None
expected=outcome()
assert actual==expected and observed==events==[] and not any(calls)
if input_kind=='valid':assert actual[0]=='value'
else:assert actual[0:2]==('error','AttributeError')
"""
    subprocess.run([sys.executable,'-c',source,scenario,input_kind],check=True,env=os.environ.copy())


@pytest.mark.parametrize('zero',[True,False])
def test_no_step_declines_before_heavy_proof_and_runs_the_original_body(monkeypatch,zero):
    t=sx.Symbol('risch_fast_t');p=sx.Poly(0 if zero else sx.Rational(7,3),t,domain='QQ')
    de=risch.DifferentialExtension(extension={'D':[sx.Poly(t*t+1,t)]})
    def forbidden(*args,**kwargs):raise AssertionError('heavy proof entered')
    monkeypatch.setattr(protocol._ground,'_class_unchanged',forbidden)
    q,r=risch.polynomial_reduce(p,de)
    assert q.is_zero and str(q.domain)=='ZZ' and r is p


@pytest.mark.parametrize('zero',[True,False])
@pytest.mark.parametrize('mode',['t_getter','d_getter','lookup','degree_error','constructor','constructor_token'])
def test_fast_refusal_preserves_real_getters_errors_and_constructor_prefix(monkeypatch,zero,mode):
    t=sx.Symbol('risch_fast_t');p=sx.Poly(0 if zero else 7,t)
    de=risch.DifferentialExtension(extension={'D':[sx.Poly(t*t+1,t)]});events=[]
    if mode in ('t_getter','d_getter'):
        name=mode[0];member=risch.DifferentialExtension.__dict__[name]
        def getter(self):events.append(name);return member.__get__(self,type(self))
        monkeypatch.setattr(risch.DifferentialExtension,name,property(getter))
    elif mode=='lookup':
        def getter(self,name):events.append(name);return object.__getattribute__(self,name)
        monkeypatch.setattr(risch.DifferentialExtension,'__getattribute__',getter)
    elif mode=='degree_error':
        original=sx.Poly.degree
        def degree(self,*args,**kwargs):
            events.append('degree')
            if self is p:raise RuntimeError('real degree error')
            return original(self,*args,**kwargs)
        monkeypatch.setattr(sx.Poly,'degree',degree)
    else:
        original=risch.Poly;token=object()
        def constructor(*args,**kwargs):
            events.append('constructor')
            return token if mode=='constructor_token' else original(*args,**kwargs)
        monkeypatch.setattr(risch,'Poly',constructor)
    def outcome():
        try:return ('value',risch.polynomial_reduce(p,de))
        except Exception as error:return ('error',type(error),str(error))
    actual=outcome();observed=events[:];events.clear()
    with monkeypatch.context() as patch:
        patch.setattr(risch,'_native_polynomial_reduce',lambda *args:None)
        old=outcome()
    assert actual==old and events==observed and events
    if mode in ('constructor','constructor_token'):assert events==['constructor']
    if mode=='degree_error':assert actual[0:2]==('error',RuntimeError)
    else:assert actual[1][1]is p


def test_payload_suggested_no_step_never_returns_a_partial_mathematical_answer(monkeypatch):
    t=sx.Symbol('risch_fast_t');p=sx.Poly(7,t)
    de=risch.DifferentialExtension(extension={'D':[sx.Poly(t*t+1,t)]})
    actual_d=sx.Poly(1,t);de.D[0]=actual_d;events=[]
    def getter(self):events.append('real d');return actual_d
    monkeypatch.setattr(risch.DifferentialExtension,'d',property(getter))
    # The original member still stores t^2+1, so the read-only hint declines.
    # The REAL getter makes the original body integrate7 instead of returning7.
    q,r=risch.polynomial_reduce(p,de);observed=events[:];events.clear()
    with monkeypatch.context()as patch:
        patch.setattr(risch,'_native_polynomial_reduce',lambda *args:None)
        old=risch.polynomial_reduce(p,de)
    assert (q,r)==old and events==observed and events
    assert q==sx.Poly(7*t,t)and r.is_zero


@pytest.mark.parametrize('registration',[False,True])
@pytest.mark.parametrize('kind',['python','python_raising','wrong_native_name'])
def test_new_preflight_origin_is_not_an_arbitrary_captured_alias(monkeypatch,registration,kind):
    p,de=input_case(19,sx,risch);events=[]
    if kind=='wrong_native_name':foreign=sx._core._risch_polynomial_reduce
    else:
        def foreign(*args):
            events.append('foreign preflight')
            if kind=='python_raising':raise AssertionError('foreign preflight called')
            return True
        foreign.__name__='_risch_decline_no_step';foreign.__self__=sx._core
    monkeypatch.setattr(sx._core,'_risch_decline_no_step',foreign)
    if registration:
        monkeypatch.setattr(protocol,'_ROOT',None)
        protocol.register_root(sx.__dict__)
        assert protocol._ROOT is None
    actual=risch.polynomial_reduce(p,de)
    with monkeypatch.context()as patch:
        patch.setattr(risch,'_native_polynomial_reduce',lambda *args:None)
        old=risch.polynomial_reduce(p,de)
    assert actual==old and events==[]


def test_uninitialized_or_shadowed_member_remains_an_unknown_hint(monkeypatch):
    t=sx.Symbol('risch_fast_t');p=sx.Poly(7,t)
    de=risch.DifferentialExtension.__new__(risch.DifferentialExtension);de.t=t
    actual_d=sx.Poly(t*t+1,t);events=[]
    def getter(self):events.append('d');return actual_d
    monkeypatch.setattr(risch.DifferentialExtension,'d',property(getter))
    actual=risch.polynomial_reduce(p,de);observed=events[:];events.clear()
    with monkeypatch.context()as patch:
        patch.setattr(risch,'_native_polynomial_reduce',lambda *args:None)
        old=risch.polynomial_reduce(p,de)
    assert actual==old and actual[1]is p and events==observed==['d']
