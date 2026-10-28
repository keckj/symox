"""Exact AUTO heuristics, public FF/GJ/CD behavior and protocol fallbacks."""
import importlib
import itertools
import os
import random
import subprocess
import sys

import pytest
import symox as sx
import sympy as sp
from symox.polys.matrices import DomainMatrix
from symox.polys.matrices.ddm import DDM
from symox.polys.matrices.sdm import SDM
from sympy.polys.matrices import DomainMatrix as RefMatrix

rref=importlib.import_module('symox.polys.matrices.rref')
guard=importlib.import_module('symox.polys.matrices._rref_native')
ref_rref=importlib.import_module('sympy.polys.matrices.rref')


def make(seed,kind,fmt):
    rng=random.Random(67305+seed)
    m,n=[(0,0),(0,9),(4,0),(1,1),(3,7),(9,10),(10,10),(12,14)][seed%8]
    bits=[1,4,17,49,73,211,900][seed%7]
    domain,refdomain=(sx.ZZ,sp.ZZ) if kind=='ZZ' else (sx.QQ,sp.QQ)
    values=[];refs=[]
    for i in range(m):
        row=[];refrow=[]
        for j in range(n):
            active=(i+j+seed)%6<seed%6
            numerator=((1<<bits)+2*rng.randrange(1,15)+1)*(-1 if rng.randrange(2) else 1) if active else 0
            d=[1,2,3,7,1<<49,1<<50,1<<400][(i+j+seed)%7] if kind=='QQ' else 1
            row.append(domain(numerator,d) if kind=='QQ' else domain(numerator))
            refrow.append(refdomain(numerator,d) if kind=='QQ' else refdomain(numerator))
        values.append(row);refs.append(refrow)
    if fmt=='sparse':
        values={i:{j:e for j,e in enumerate(row) if e} for i,row in enumerate(values) if any(row)}
        refs={i:{j:e for j,e in enumerate(row) if e} for i,row in enumerate(refs) if any(row)}
    return DomainMatrix(values,(m,n),domain),RefMatrix(refs,(m,n),refdomain)


@pytest.mark.parametrize('seed',range(48))
@pytest.mark.parametrize('kind',['ZZ','QQ'])
@pytest.mark.parametrize('fmt',['dense','sparse'])
@pytest.mark.parametrize('denominator',[False,True])
def test_auto_decisions_match_complete_old_and_reference(monkeypatch,seed,kind,fmt,denominator):
    a,ref=make(seed,kind,fmt);calls=[];native=sx._core._dm_rref_auto
    def spy(*args):
        result=native(*args);calls.append(result);return result
    monkeypatch.setattr(sx._core,'_dm_rref_auto',spy)
    got=rref._dm_rref_choose_method(a,'auto',denominator=denominator)
    old=(rref._dm_rref_choose_method_ZZ if kind=='ZZ' else rref._dm_rref_choose_method_QQ)(a,denominator=denominator)
    expected=ref_rref._dm_rref_choose_method(ref,'auto',denominator=denominator)
    assert got==(old,'sparse')==expected
    assert calls==[old]


def public_case(method,fmt):
    kind='ZZ' if method in ('ZZ-GJ','ZZ-FF') else 'QQ'
    k,rk=(sx.ZZ,sp.ZZ) if kind=='ZZ' else (sx.QQ,sp.QQ)
    m,n=(3,9) if method=='ZZ-GJ' else (3,4)
    rows=[];refs=[]
    for i in range(m):
        row=[];refrow=[]
        for j in range(n):
            p=(i+1)**(j+1)+i+j+3
            if method=='ZZ-GJ' and j not in (i,i+3):p=0
            d={'QQ-CD':6,'QQ-FF':1<<55,'QQ-GJ':1<<400}.get(method,1)
            if method=='QQ-FF':p+=(1<<60)+1
            row.append(k(p,d) if kind=='QQ' else k(p));refrow.append(rk(p,d) if kind=='QQ' else rk(p))
        rows.append(row);refs.append(refrow)
    if fmt=='sparse':
        rows={i:{j:e for j,e in enumerate(row) if e} for i,row in enumerate(rows) if any(row)}
        refs={i:{j:e for j,e in enumerate(row) if e} for i,row in enumerate(refs) if any(row)}
    return DomainMatrix(rows,(m,n),k),RefMatrix(refs,(m,n),rk)


def rational(element):
    return sx.to_sympy(element.as_expr()) if type(element) is sx._core.DomainElement else sp.Rational(element)


@pytest.mark.parametrize('case',['ZZ-GJ','ZZ-FF','QQ-GJ','QQ-FF','QQ-CD'])
@pytest.mark.parametrize('fmt',['dense','sparse'])
@pytest.mark.parametrize('operation',['rref','rref_den','rref_den_field'])
def test_public_auto_uses_native_decision_and_preserves_every_return(monkeypatch,case,fmt,operation):
    a,reference=public_case(case,fmt);native=sx._core._dm_rref_auto;calls=[]
    before=a.to_list();row_ids=[id(row) for row in a.rep] if fmt=='dense' else [id(row) for row in a.rep.values()]
    def spy(*args):
        result=native(*args);calls.append(result);return result
    monkeypatch.setattr(sx._core,'_dm_rref_auto',spy)
    options={} if operation=='rref' else {'keep_domain':operation=='rref_den'}
    name='rref' if operation=='rref' else 'rref_den'
    actual=getattr(a,name)(**options);expected=getattr(reference,name)(**options)
    assert calls==[case.split('-')[1]]
    reduced=actual[0];refreduced=expected[0]
    assert type(reduced) is DomainMatrix and reduced.rep.fmt==fmt
    assert str(reduced.domain)==str(refreduced.domain)
    assert reduced.shape==refreduced.shape and actual[-1]==expected[-1]
    got=[[rational(e) for e in row] for row in reduced.to_list()]
    assert got==[[sp.Rational(e) for e in row] for row in refreduced.to_list()]
    if name=='rref_den':
        assert rational(actual[1])==sp.Rational(expected[1])
        assert actual[1].parent()==reduced.domain
    # Exact row-space certificate; avoids depending only on another RREF.
    original=sp.Matrix([[rational(e) for e in row] for row in before])
    final=sp.Matrix(got)
    assert original.col_join(final).rank()==original.rank()==len(actual[-1])
    assert a.to_list()==before
    assert row_ids==([id(row) for row in a.rep] if fmt=='dense' else [id(row) for row in a.rep.values()])
    monkeypatch.setattr(rref,'_native_choose_method',lambda *args:None)
    legacy=getattr(a,name)(**options)
    assert actual==legacy and reduced.domain==legacy[0].domain


@pytest.mark.parametrize('kind',['ZZ','QQ'])
@pytest.mark.parametrize('cols',[0,1,5,13])
def test_literal_sparse_empty_rows_are_not_removed_from_statistics(kind,cols):
    k,rk=(sx.ZZ,sp.ZZ) if kind=='ZZ' else (sx.QQ,sp.QQ)
    data={i:{} for i in range(12)}
    a=DomainMatrix(data,(12,cols),k);b=RefMatrix(data,(12,cols),rk)
    assert rref._dm_row_density(a)==ref_rref._dm_row_density(b)
    assert guard.native_choose_method(a)==ref_rref._dm_rref_choose_method(b,'auto')[0]


@pytest.mark.parametrize('kind',['ZZ','QQ'])
def test_literal_sparse_zero_coefficients_preserve_choice_or_exact_error(kind):
    k,rk=(sx.ZZ,sp.ZZ) if kind=='ZZ' else (sx.QQ,sp.QQ)
    a=DomainMatrix({i:{j:k.zero for j in range(6)} for i in range(10)},(10,10),k)
    b=RefMatrix({i:{j:rk.zero for j in range(6)} for i in range(10)},(10,10),rk)
    if kind=='ZZ':
        assert guard.native_choose_method(a) is None
        with pytest.raises(ZeroDivisionError):rref._dm_rref_choose_method(a,'auto')
        with pytest.raises(ZeroDivisionError):ref_rref._dm_rref_choose_method(b,'auto')
    else:
        assert guard.native_choose_method(a)=='GJ'
        assert rref._dm_rref_choose_method(a,'auto')==ref_rref._dm_rref_choose_method(b,'auto')


@pytest.mark.parametrize('method',['GJ','FF','CD','GJ_dense','FF_dense','CD_dense','unknown'])
def test_explicit_method_options_do_not_touch_native_auto(monkeypatch,method):
    a,_=public_case('QQ-CD','dense')
    def forbidden(*args):raise AssertionError('explicit selection invoked AUTO')
    monkeypatch.setattr(sx._core,'_dm_rref_auto',forbidden)
    expected=(method[:-6],'dense') if method.endswith('_dense') else (method,'sparse')
    assert rref._dm_rref_choose_method(a,method)==expected


@pytest.mark.parametrize('kind',['ZZ','QQ'])
def test_sparse_shortcuts_do_not_probe_foreign_coefficients(kind):
    k=sx.ZZ if kind=='ZZ' else sx.QQ
    a=DomainMatrix({0:{0:k.one}},(1,20),k);calls=[]
    class Foreign:
        def __int__(self):calls.append('int');raise AssertionError('early integer probe')
        @property
        def numerator(self):calls.append('numerator');raise AssertionError('early numerator probe')
        def __bool__(self):calls.append('bool');raise AssertionError('early boolean probe')
    a.rep[0][0]=Foreign()
    assert guard.native_choose_method(a)=='GJ'
    assert rref._dm_rref_choose_method(a,'auto')==('GJ','sparse') and calls==[]


def observe_fallback(monkeypatch,a,events):
    assert guard.native_choose_method(a) is None and events==[]
    actual=rref._dm_rref_choose_method(a,'auto');observed=events[:];events.clear()
    with monkeypatch.context() as patch:
        patch.setattr(rref,'_native_choose_method',lambda *args:None)
        expected=rref._dm_rref_choose_method(a,'auto')
    assert actual==expected and events==observed
    return observed


@pytest.mark.parametrize('name',['_dm_row_density','_dm_elements','_dm_QQ_numers_denoms',
    '_dm_rref_choose_method_QQ','_dm_rref_choose_method_ZZ'])
def test_heuristic_overrides_keep_complete_callback_order(monkeypatch,name):
    a,_=make(14,'ZZ' if name.endswith('_ZZ') else 'QQ','dense')
    if name=='_dm_rref_choose_method_ZZ':a,_=make(14,'ZZ','dense')
    # Dense enough to read the coefficient statistics.
    if name!='_dm_rref_choose_method_ZZ':a,_=public_case('QQ-CD','dense')
    events=[];original=getattr(rref,name)
    def callback(*args,**kwargs):events.append(name);return original(*args,**kwargs)
    monkeypatch.setattr(rref,name,callback)
    assert observe_fallback(monkeypatch,a,events)


@pytest.mark.parametrize('name',['sum','map','len','max','min','int'])
def test_heuristic_builtin_overrides_are_not_called_speculatively(monkeypatch,name):
    import builtins
    a,_=public_case('QQ-CD','dense');events=[];original=getattr(builtins,name)
    def callback(*args,**kwargs):events.append(name);return original(*args,**kwargs)
    monkeypatch.setattr(rref,name,callback,raising=False)
    assert observe_fallback(monkeypatch,a,events)


@pytest.mark.parametrize('target,name',[(DDM,'to_sdm'),(DDM,'to_flat_nz'),(SDM,'to_sdm'),
    (SDM,'to_flat_nz'),(SDM,'to_dok'),(DomainMatrix,'to_flat_nz')])
def test_matrix_method_overrides_preserve_complete_callback_order(monkeypatch,target,name):
    a,_=public_case('QQ-CD','sparse' if target is SDM else 'dense');events=[];original=getattr(target,name)
    def callback(*args,**kwargs):events.append(name);return original(*args,**kwargs)
    monkeypatch.setattr(target,name,callback)
    assert observe_fallback(monkeypatch,a,events)


@pytest.mark.parametrize('name',['__bool__','numerator','denominator'])
def test_numeric_descriptor_hooks_preserve_order_without_preflight(monkeypatch,name):
    a,_=public_case('QQ-CD','dense');events=[];cls=sx._core.DomainElement;original=cls.__dict__[name]
    def callback(self):
        events.append(name)
        return original.__get__(self,cls)() if name=='__bool__' else original.__get__(self,cls)
    monkeypatch.setattr(cls,name,callback if name=='__bool__' else property(callback))
    assert observe_fallback(monkeypatch,a,events)


def test_integer_conversion_hook_is_never_skipped(monkeypatch):
    # Ten rows avoid the small-matrix shortcut and force the integer bit count.
    a,_=make(46,'ZZ','dense');events=[];cls=sx._core.DomainElement;original=cls.__int__
    def callback(self):events.append('int');return original(self)
    monkeypatch.setattr(cls,'__int__',callback)
    assert observe_fallback(monkeypatch,a,events)


@pytest.mark.parametrize('name',['is_ZZ','is_QQ','of_type'])
def test_domain_predicate_hooks_keep_old_order(monkeypatch,name):
    a,_=public_case('QQ-CD','dense');events=[];cls=sx._core.Domain;original=cls.__dict__[name]
    if name=='of_type':
        def callback(self,element):events.append(name);return original.__get__(self,cls)(element)
        changed=callback
    else:
        def callback(self):events.append(name);return original.__get__(self,cls)
        changed=property(callback)
    monkeypatch.setattr(cls,name,changed)
    assert observe_fallback(monkeypatch,a,events)


def test_dense_constructor_native_elements_hook_is_not_bypassed(monkeypatch):
    module=importlib.import_module('symox.polys.matrices.sdm');a,_=public_case('QQ-CD','dense')
    original=module._native_elements;events=[]
    def callback(*args):events.append('native_elements');return original(*args)
    monkeypatch.setattr(module,'_native_elements',callback)
    assert observe_fallback(monkeypatch,a,events)


def test_dense_convert_descriptor_getter_retains_its_callbacks(monkeypatch):
    a,_=public_case('QQ-CD','dense');cls=sx._core.Domain;original=cls.__dict__['convert'];events=[]
    def getter(self):events.append('convert-get');return original.__get__(self,cls)
    monkeypatch.setattr(cls,'convert',property(getter))
    assert observe_fallback(monkeypatch,a,events)


def test_unused_convert_function_is_not_called_for_already_typed_dense_entries(monkeypatch):
    a,_=public_case('QQ-CD','dense');events=[]
    def forbidden(*args):events.append('convert');raise AssertionError('conversion of typed native input')
    monkeypatch.setattr(sx._core.Domain,'convert',forbidden)
    assert guard.native_choose_method(a)=='CD' and events==[]
    assert rref._dm_rref_choose_method_QQ(a)=='CD' and events==[]


def test_dense_untyped_entries_keep_the_conversion_callback(monkeypatch):
    a,_=public_case('QQ-CD','dense');a.rep[0][0]=3;cls=sx._core.Domain;original=cls.convert;events=[]
    def callback(self,*args,**kwargs):events.append('convert');return original(self,*args,**kwargs)
    monkeypatch.setattr(cls,'convert',callback)
    assert observe_fallback(monkeypatch,a,events)


@pytest.mark.parametrize('mode',['row','key'])
def test_sparse_foreign_container_and_key_protocols_remain_observable(monkeypatch,mode):
    a,_=public_case('QQ-CD','sparse');events=[]
    if mode=='row':
        class Row(dict):
            def __len__(self):events.append('row-len');return super().__len__()
        a.rep[0]=Row(a.rep[0])
    else:
        class Key(int):
            def __hash__(self):events.append('key-hash');return super().__hash__()
        a.rep[0]={Key(j):value for j,value in a.rep[0].items()};events.clear()
    assert observe_fallback(monkeypatch,a,events)


def test_host_integer_callbacks_are_not_probed_before_full_statistics(monkeypatch):
    a,_=make(46,'ZZ','sparse');events=[]
    class Scalar:
        def __int__(self):events.append('int');return 3
    a.rep[0][0]=Scalar()
    assert observe_fallback(monkeypatch,a,events)


@pytest.mark.parametrize('target', ['heuristic','numeric'])
def test_interruptions_are_not_transformed_into_a_method_or_missing_solution(monkeypatch,target):
    a,_=public_case('QQ-CD','dense')
    def interrupt(*args,**kwargs):raise KeyboardInterrupt('AUTO callback')
    if target=='heuristic':monkeypatch.setattr(rref,'_dm_row_density',interrupt)
    else:monkeypatch.setattr(sx._core.DomainElement,'numerator',property(interrupt))
    assert guard.native_choose_method(a) is None
    with pytest.raises(KeyboardInterrupt,match='AUTO callback'):rref._dm_rref_choose_method(a,'auto')


@pytest.mark.parametrize('kind',['ZZ','QQ'])
@pytest.mark.parametrize('cols',[2**53+1,2**63+7,2**64-1,10**400])
def test_large_sparse_metadata_keeps_exact_decisions_or_original_overflow(kind,cols):
    k,rk=(sx.ZZ,sp.ZZ) if kind=='ZZ' else (sx.QQ,sp.QQ)
    a=DomainMatrix({0:{0:k.one}},(1,cols),k);b=RefMatrix({0:{0:rk.one}},(1,cols),rk)
    if cols==10**400:
        assert guard.native_choose_method(a) is None
        with pytest.raises(OverflowError):rref._dm_rref_choose_method(a,'auto')
        with pytest.raises(OverflowError):ref_rref._dm_rref_choose_method(b,'auto')
    else:
        assert guard.native_choose_method(a)=='GJ'
        assert rref._dm_rref_choose_method(a,'auto')==ref_rref._dm_rref_choose_method(b,'auto')


@pytest.mark.parametrize('order',list(itertools.permutations(('stats','vector','holonomic.recurrence'))))
def test_optional_collectors_and_stock_rref_reloads_keep_native_selection(order):
    code='''
import importlib,sys
import symox as s
from symox.polys.matrices import DomainMatrix
from symox.polys.matrices._rref_native import native_choose_method
r=importlib.import_module('symox.polys.matrices.rref')
a=DomainMatrix([[s.QQ(1,2),s.QQ(1,3)],[s.QQ(3,7),s.QQ(4,7)]],(2,2),s.QQ)
assert native_choose_method(a)=='CD'
for name in sys.argv[1:]:importlib.import_module('symox.'+name)
for repeat in range(3):
 importlib.reload(r)
 assert native_choose_method(a)=='CD'
 assert r._dm_rref_choose_method(a,'auto')==('CD','sparse')
'''
    subprocess.run([sys.executable,'-c',code,*order],check=True,env=os.environ.copy())


@pytest.mark.parametrize('mode',['builtin','helper'])
def test_rref_reload_does_not_bless_user_callbacks(mode):
    code='''
import importlib
import symox as s
from symox.polys.matrices import DomainMatrix
from symox.polys.matrices._rref_native import native_choose_method
r=importlib.import_module('symox.polys.matrices.rref')
a=DomainMatrix([[s.QQ(1,2),s.QQ(1,3)],[s.QQ(3,7),s.QQ(4,7)]],(2,2),s.QQ)
def callback(*args,**kwargs):raise KeyboardInterrupt('kept user callback')
'''
    if mode=='builtin':code+='r.int=callback\nimportlib.reload(r)\n'
    else:code+='importlib.reload(r)\nr._dm_row_density=callback\n'
    code+='''
assert native_choose_method(a) is None
try:r._dm_rref_choose_method(a,'auto')
except KeyboardInterrupt as error:assert str(error)=='kept user callback'
else:raise AssertionError('AUTO skipped callback')
'''
    subprocess.run([sys.executable,'-c',code],check=True,env=os.environ.copy())


@pytest.mark.parametrize('name',['lcm','one'])
def test_lcm_and_initial_unit_hooks_remain_observable(monkeypatch,name):
    a,_=public_case('QQ-CD','sparse');cls=sx._core.Domain;original=cls.__dict__[name];events=[]
    if name=='one':
        def callback(self):events.append(name);return original.__get__(self,cls)
        changed=property(callback)
    else:
        def callback(self,*args):events.append(name);return original.__get__(self,cls)(*args)
        changed=callback
    monkeypatch.setattr(cls,name,changed)
    assert observe_fallback(monkeypatch,a,events)


def test_dense_sdm_super_constructor_hook_is_not_bypassed(monkeypatch):
    import builtins
    module=importlib.import_module('symox.polys.matrices.sdm');a,_=public_case('QQ-CD','dense');events=[]
    # Invoke the exact superclass with the original initializer's instance.
    # Only its global super lookup changes; no method is replaced here.
    def callback():
        events.append('super')
        return builtins.super(SDM,sys._getframe(1).f_locals['self'])
    monkeypatch.setattr(module,'super',callback,raising=False)
    assert observe_fallback(monkeypatch,a,events)


@pytest.mark.parametrize('fmt',['dense','sparse'])
def test_unrelated_matrix_attributes_do_not_disable_native_selection(monkeypatch,fmt):
    a,_=public_case('QQ-CD',fmt)
    monkeypatch.setattr(DomainMatrix,'unrelated_plugin_marker',True,raising=False)
    monkeypatch.setattr(DDM,'unrelated_plugin_marker',True,raising=False)
    monkeypatch.setattr(SDM,'unrelated_plugin_marker',True,raising=False)
    a.unrelated_payload=object()
    assert guard.native_choose_method(a)=='CD'


def test_instance_method_override_retains_its_call_order(monkeypatch):
    a,_=public_case('QQ-CD','dense');original=a.to_flat_nz;events=[]
    def callback():events.append('instance-flat');return original()
    monkeypatch.setattr(a,'to_flat_nz',callback)
    assert observe_fallback(monkeypatch,a,events)


@pytest.mark.parametrize('mode',['row_length','outer_length','shape','domain'])
def test_invalid_or_distinct_metadata_keeps_full_old_behavior(monkeypatch,mode):
    a,_=public_case('QQ-CD','dense')
    if mode=='row_length':a.rep[0].pop()
    elif mode=='outer_length':a.rep.pop()
    elif mode=='shape':a.shape=tuple(list(a.shape))
    else:a.rep.domain=sx.ZZ
    assert guard.native_choose_method(a) is None
    def outcome():
        try:return ('value',rref._dm_rref_choose_method(a,'auto'))
        except Exception as error:return ('error',type(error),str(error))
    actual=outcome()
    monkeypatch.setattr(rref,'_native_choose_method',lambda *args:None)
    assert outcome()==actual


def test_reloading_domainmatrix_does_not_bless_replaced_numeric_protocols():
    code='''
import importlib
import symox as s
r=importlib.import_module('symox.polys.matrices.rref')
dm=importlib.import_module('symox.polys.matrices.domainmatrix')
def callback(*args):raise KeyboardInterrupt('numeric callback preserved')
s._core.DomainElement.numerator=property(callback)
importlib.reload(dm)
a=dm.DomainMatrix([[s.QQ(1,2),s.QQ(1,3)],[s.QQ(3,7),s.QQ(4,7)]],(2,2),s.QQ)
try:r._dm_rref_choose_method(a,'auto')
except KeyboardInterrupt as error:assert str(error)=='numeric callback preserved'
else:raise AssertionError('class reload blessed callback')
'''
    subprocess.run([sys.executable,'-c',code],check=True,env=os.environ.copy())


@pytest.mark.parametrize('kind',['RR','CC','EX','GF','QQ[y]','QQ(y)'])
@pytest.mark.parametrize('denominator',[False,True])
def test_other_domains_keep_complete_selection_and_original_payloads(kind,denominator,monkeypatch):
    t=sx.Symbol('t')
    k={'RR':lambda:sx.RealField(211),'CC':lambda:sx.ComplexField(140),'EX':lambda:sx.EX,
       'GF':lambda:sx.GF(5),'QQ[y]':lambda:sx.QQ.poly_ring(t),'QQ(y)':lambda:sx.QQ.frac_field(t)}[kind]()
    a=DomainMatrix([[k.one,k.zero],[k.zero,k.one]],(2,2),k)
    values=[value for row in a.rep for value in row]
    assert guard.native_choose_method(a) is None
    got=rref._dm_rref_choose_method(a,'auto',denominator=denominator)
    monkeypatch.setattr(rref,'_native_choose_method',lambda *args:None)
    assert got==rref._dm_rref_choose_method(a,'auto',denominator=denominator)
    assert all(old is new for old,new in zip(values,(value for row in a.rep for value in row)))
