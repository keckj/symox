"""Exact QQ sparse Jordan elimination, ordered metadata and full fallbacks."""
import importlib
import itertools
import os
import random
import subprocess
import sys
from fractions import Fraction

import pytest
import symox as sx
import sympy as sp
from symox.polys.matrices import DomainMatrix
from symox.polys.matrices.sdm import SDM
from sympy.polys.matrices.sdm import sdm_irref as reference_irref

module=importlib.import_module('symox.polys.matrices.sdm')
protocol=importlib.import_module('symox.polys.matrices._sdm_irref_native')


def rational(value):
    return Fraction(int(value.numerator),int(value.denominator))


def signature(result):
    rows,pivots,columns=result
    return ([(i,[(j,rational(v))for j,v in row.items()])for i,row in rows.items()],
            list(pivots),[(j,list(indices))for j,indices in columns.items()])


def aliases(a,result):
    objects=[]
    for i,row in a.items():
        objects.extend((i,row))
        for j,value in row.items():objects.extend((j,value))
    rows,pivots,columns=result
    for i,row in rows.items():
        objects.extend((i,row))
        for j,value in row.items():objects.extend((j,value))
    objects.extend(pivots)
    for j,indices in columns.items():
        objects.extend((j,indices));objects.extend(sorted(indices))
    identities={};return_value=[]
    for obj in objects:
        identity=id(obj)
        if identity not in identities:identities[identity]=len(identities)
        return_value.append(identities[identity])
    return return_value


def generated(seed):
    rng=random.Random(512097+seed)
    width=seed%10;count=(seed*7)%9
    pool=[0,1,8,16,31,32,63,64,65,127,128,256,257,1023,1024,1025,
          (1<<31)+1,(1<<53)+7,(1<<61)-1,2*((1<<61)-1),3*((1<<61)-1),(1<<64)-1]
    rng.shuffle(pool);keys=sorted(pool[:width]);a={};ref={};dense=[]
    for i in range(count):
        order=keys[:];rng.shuffle(order);row={};reference={};values={}
        for j in order:
            if rng.randrange(5)>seed%4:
                p=rng.randrange(-17,18)
                if p:
                    p*=((1<<211)+1) if seed%13==0 else 1
                    q=rng.randrange(1,17)
                    # Distinct equal key objects expose representative aliases.
                    key=int(str(j));row[key]=sx.QQ(p,q);reference[key]=sp.QQ(p,q);values[j]=sp.Rational(p,q)
        if row:
            outer=10000+17*i if seed%2 else -i
            a[outer]=row;ref[outer]=reference;dense.append([values.get(j,sp.S.Zero)for j in keys])
    if not dense:dense=sp.zeros(0,len(keys))
    else:dense=sp.Matrix(dense)
    return a,ref,keys,dense


@pytest.mark.parametrize('seed',range(120))
def test_complete_reduction_matches_old_reference_metadata_aliases_and_rank(monkeypatch,seed):
    a,reference,keys,dense=generated(seed);before=[(id(row),[(id(j),id(value))for j,value in row.items()])for row in a.values()]
    native=sx._core._sdm_irref;calls=[]
    def spy(*args):
        result=native(*args);calls.append(result is not None);return result
    monkeypatch.setattr(sx._core,'_sdm_irref',spy)
    got=module.sdm_irref(a)
    assert calls==[bool(a)]
    expected=reference_irref(reference)
    assert signature(got)==signature(expected)
    with monkeypatch.context() as patch:
        patch.setattr(module,'_native_irref',lambda *args:None)
        old=module.sdm_irref(a)
    assert signature(got)==signature(old) and aliases(a,got)==aliases(a,old)
    assert before==[(id(row),[(id(j),id(value))for j,value in row.items()])for row in a.values()]
    reduced,pivots,_=got
    compressed=sp.Matrix([[sp.Rational(rational(row.get(j,sx.QQ.zero)).numerator,rational(row.get(j,sx.QQ.zero)).denominator)for j in keys]for row in reduced.values()]) if reduced else sp.zeros(0,len(keys))
    assert len(pivots)==dense.rank()==compressed.rank()
    assert dense.col_join(compressed).rank()==len(pivots)
    assert [keys[i]for i in dense.rref()[1]]==pivots


@pytest.mark.parametrize('values',[[(0,1,1)],[(1000,1,1),(2000,2,3)],[(2**61-1,-7,11),(2*(2**61-1),13,17)]])
def test_single_row_normalizes_every_scalar_without_reusing_input_values(monkeypatch,values):
    a={77:{int(str(j)):sx.QQ(p,q)for j,p,q in values}}
    actual=module.sdm_irref(a)
    with monkeypatch.context() as patch:
        patch.setattr(module,'_native_irref',lambda *args:None);old=module.sdm_irref(a)
    assert signature(actual)==signature(old) and aliases(a,actual)==aliases(a,old)
    assert all(new is not original for row in actual[0].values()for new in row.values()for original in a[77].values())


def test_reverse_stable_ties_keep_original_large_key_objects():
    p0,p1=int('1000'),int('1000');q0,q1=int('2000'),int('2000')
    a={17:{p0:sx.QQ(1),q0:sx.QQ(2)},-5:{p1:sx.QQ(1),q1:sx.QQ(2)}}
    rows,pivots,columns=module.sdm_irref(a)
    assert pivots[0] is p1 and list(rows[0])[0] is p1 and list(rows[0])[1] is q1
    assert next(iter(columns)) is q1


def test_large_row_number_metadata_aliases_match_old_without_dense_storage(monkeypatch):
    a={i:{i:sx.QQ(1),10000:sx.QQ(i+1),20000:sx.QQ(i+2)}for i in range(270)}
    result=module.sdm_irref(a)
    with monkeypatch.context() as patch:
        patch.setattr(module,'_native_irref',lambda *args:None);old=module.sdm_irref(a)
    assert signature(result)==signature(old) and aliases(a,result)==aliases(a,old)
    one=next(value for value in result[2][10000]if value==269)
    two=next(value for value in result[2][20000]if value==269)
    outer=next(value for value in result[0]if value==269)
    assert one is two and one is not outer and one is not result[1][269]


@pytest.mark.parametrize('a',[{}, {0:{}}, {0:{0:0}}, {0:{0:1,1:0}}, {0:{0:1},1:{}},
    {0:{-1:1,2:3}}, {0:{10**100:1,10**101:2}}, {0:{True:1}}])
def test_unhandled_storage_preserves_exact_old_result_or_error(monkeypatch,a):
    a={i:{j:sx.QQ(v)for j,v in row.items()}for i,row in a.items()}
    assert module._native_irref(a) is None
    def outcome():
        try:return ('result',module.sdm_irref(a))
        except Exception as error:return ('error',type(error),str(error))
    actual=outcome()
    monkeypatch.setattr(module,'_native_irref',lambda *args:None)
    old=outcome()
    if actual[0]=='result':assert signature(actual[1])==signature(old[1])
    else:assert actual==old


@pytest.mark.parametrize('kind',['ZZ','RR','CC','GF','QQ_I','QQ(y)','EX'])
def test_other_domains_retain_the_complete_general_algorithm(monkeypatch,kind):
    y=sx.Symbol('y')
    k={'ZZ':lambda:sx.ZZ,'RR':lambda:sx.RealField(140),'CC':lambda:sx.ComplexField(211),
       'GF':lambda:sx.GF(5),'QQ_I':lambda:sx.QQ_I,'QQ(y)':lambda:sx.QQ.frac_field(y),'EX':lambda:sx.EX}[kind]()
    a={0:{0:k(1),1:k(2)},1:{0:k(3),1:k(4)}}
    # GF(p) and QQ_I are exact fields with native arithmetic since wave 18.
    native=module._native_irref(a)
    assert (native is not None) if kind in ('GF','QQ_I') else (native is None)
    def outcome():
        try:return ('result',module.sdm_irref(a))
        except Exception as error:return ('error',type(error),str(error))
    actual=outcome();monkeypatch.setattr(module,'_native_irref',lambda *args:None);old=outcome()
    if actual[0]=='result':
        assert actual[1]==old[1]
    else:assert actual==old


def coefficients():
    return {13:{0:sx.QQ(1,2),2:sx.QQ(2,3),5:sx.QQ(5,7)},
            7:{1:sx.QQ(3,5),2:sx.QQ(7,11),5:sx.QQ(-2,7)},
            41:{0:sx.QQ(2,3),1:sx.QQ(1,7),4:sx.QQ(4,9)}}


def fallback_matches(monkeypatch,a,events):
    assert module._native_irref(a) is None and events==[]
    actual=module.sdm_irref(a);observed=events[:];events.clear()
    with monkeypatch.context() as patch:
        patch.setattr(module,'_native_irref',lambda *args:None)
        old=module.sdm_irref(a)
    legacy=events[:]
    assert observed==legacy
    return actual,old,observed


@pytest.mark.parametrize('name',['__pow__','__mul__','__sub__','__neg__','__bool__'])
def test_consumed_coefficient_hooks_keep_the_complete_old_sequence(monkeypatch,name):
    a=coefficients();events=[];cls=sx._core.DomainElement;original=getattr(cls,name)
    def callback(*args):events.append(name);return original(*args)
    monkeypatch.setattr(cls,name,callback)
    actual,old,observed=fallback_matches(monkeypatch,a,events)
    assert observed and signature(actual)==signature(old)


def test_added_inplace_multiplication_is_observed(monkeypatch):
    a=coefficients();events=[];cls=sx._core.DomainElement;original=cls.__mul__
    def callback(self,other):events.append('imul');return original(self,other)
    monkeypatch.setattr(cls,'__imul__',callback,raising=False)
    actual,old,observed=fallback_matches(monkeypatch,a,events)
    assert observed and signature(actual)==signature(old)


@pytest.mark.parametrize('name',['__getattribute__','__getattr__'])
def test_element_attribute_hooks_are_not_probed_during_selection(monkeypatch,name):
    a=coefficients();events=[];cls=sx._core.DomainElement;original=cls.__getattribute__
    def callback(self,key):
        events.append(key)
        if name=='__getattr__':raise AttributeError(key)
        return original(self,key)
    monkeypatch.setattr(cls,name,callback,raising=False)
    actual,old,observed=fallback_matches(monkeypatch,a,events)
    assert observed
    # Restore before examining coefficients: those examinations have their own
    # attribute lookups and do not belong to the elimination contract.
    monkeypatch.undo()
    assert signature(actual)==signature(old)


@pytest.mark.parametrize('marker',['_symox_domainmatrix','_symox_quotient_element','is_Matrix','_symox_sparse','_symox_ring_element'])
def test_descriptor_temporaries_and_restoration_do_not_become_trusted(monkeypatch,marker):
    a=coefficients();events=[];cls=sx._core.DomainElement;temporary='_sdm_hostile_temporary'
    class Flag:
        def __bool__(self):events.append('bool');return False
        def __del__(self):
            events.append('restore')
            if temporary in cls.__dict__:delattr(cls,temporary)
    def getter(self):
        events.append('get');setattr(cls,temporary,True);return Flag()
    monkeypatch.setattr(cls,marker,property(getter),raising=False)
    actual,old,observed=fallback_matches(monkeypatch,a,events)
    assert observed and 'get' in observed and 'restore' in observed
    assert temporary not in cls.__dict__ and signature(actual)==signature(old)


def test_dynamic_native_element_finalizer_is_not_consumed_by_cpp_allocations(monkeypatch):
    a=coefficients();events=[]
    def finalizer(self):events.append('del')
    with monkeypatch.context() as patch:
        patch.setattr(sx._core.DomainElement,'__del__',finalizer,raising=False)
        assert module._native_irref(a) is not None
        actual=module.sdm_irref(a)
        saved=module._native_irref
        try:
            module._native_irref=lambda *args:None
            old=module.sdm_irref(a)
        finally:module._native_irref=saved
        assert events==[]
    assert signature(actual)==signature(old)


@pytest.mark.parametrize('name',['sorted','min','set','dict','enumerate','len','defaultdict'])
def test_used_factory_and_builtin_callbacks_are_not_bypassed(monkeypatch,name):
    import builtins,collections
    a=coefficients();events=[];original=collections.defaultdict if name=='defaultdict' else getattr(builtins,name)
    def callback(*args,**kwargs):events.append(name);return original(*args,**kwargs)
    monkeypatch.setattr(module,name,callback,raising=False)
    actual,old,observed=fallback_matches(monkeypatch,a,events)
    assert observed and signature(actual)==signature(old)


@pytest.mark.parametrize('kind',['values','copy','key','scalar'])
def test_custom_containers_keys_and_scalars_keep_the_old_callbacks(monkeypatch,kind):
    a=coefficients();events=[]
    if kind=='values':
        class Outer(dict):
            def values(self):events.append('values');return super().values()
        a=Outer(a)
    elif kind=='copy':
        class Row(dict):
            def copy(self):events.append('copy');return super().copy()
        a[13]=Row(a[13])
    elif kind=='key':
        class Key(int):
            def __hash__(self):events.append('hash');return super().__hash__()
        a[13]={Key(j):v for j,v in a[13].items()};events.clear()
    else:
        class Scalar:
            def __bool__(self):events.append('bool');return True
            def __pow__(self,n):events.append('pow');return sx.QQ(2)**n
        # One pivot only avoids requiring a general host algebra in the fixture.
        a={0:{0:Scalar()}}
    assert module._native_irref(a) is None and events==[]
    if kind=='scalar':
        def outcome():
            try:return ('value',module.sdm_irref(a))
            except Exception as error:return ('error',type(error),str(error))
        actual=outcome();observed=events[:];events.clear()
        monkeypatch.setattr(module,'_native_irref',lambda *args:None)
        assert outcome()==actual and events==observed and observed
    else:
        actual,old,observed=fallback_matches(monkeypatch,a,events)
        assert observed and signature(actual)==signature(old)


@pytest.mark.parametrize('kind',['division','add','domain_getattribute','domain_getattr','basic_instancecheck','import'])
def test_unconsumed_protocols_do_not_cause_false_refusal_or_callbacks(monkeypatch,kind):
    import builtins
    a=coefficients();events=[]
    def forbidden(*args,**kwargs):events.append(kind);raise AssertionError('unconsumed callback')
    with monkeypatch.context() as patch:
        if kind=='division':patch.setattr(sx._core.DomainElement,'__truediv__',forbidden)
        elif kind=='add':patch.setattr(sx._core.DomainElement,'__add__',forbidden)
        elif kind=='domain_getattribute':patch.setattr(sx._core.Domain,'__getattribute__',forbidden,raising=False)
        elif kind=='domain_getattr':patch.setattr(sx._core.Domain,'__getattr__',forbidden,raising=False)
        elif kind=='basic_instancecheck':
            basic=importlib.import_module('symox._basic');patch.setattr(basic._BasicMeta,'__instancecheck__',forbidden)
        else:patch.setattr(builtins,'__import__',forbidden)
        assert module._native_irref(a) is not None
        actual=module.sdm_irref(a)
        saved=module._native_irref
        try:
            module._native_irref=lambda *args:None
            old=module.sdm_irref(a)
        finally:module._native_irref=saved
        assert events==[]
    assert signature(actual)==signature(old)


@pytest.mark.parametrize('mode',['cold_cache','stale_token'])
def test_basic_abc_cache_state_is_not_consumed_by_native_qq_arithmetic(mode):
    source="""
import _abc,abc
import symox as s
from symox.polys.matrices.sdm import _native_irref,sdm_irref
_abc._reset_caches(s.Basic)
"""
    if mode=='stale_token':source+="""
class Unrelated(metaclass=abc.ABCMeta):pass
Unrelated.register(type('Registered',(),{}))
"""
    source+="""
a={0:{0:s.QQ(1),2:s.QQ(2)},1:{1:s.QQ(3),2:s.QQ(4)}}
assert _native_irref(a) is not None
assert len(sdm_irref(a)[1])==2
"""
    subprocess.run([sys.executable,'-c',source],check=True,env=os.environ.copy())


@pytest.mark.parametrize('fmt',['dense','sparse'])
@pytest.mark.parametrize('operation',['rref','rref_den','rref_den_field'])
@pytest.mark.parametrize('seed',range(5))
def test_public_auto_gj_uses_native_elimination_with_full_exact_certificates(monkeypatch,fmt,operation,seed):
    rows,cols=5,12;rng=random.Random(7801+seed);data={}
    for i in range(rows):
        data[i]={i:sx.QQ(rng.randrange(1,10),rng.randrange(1,10)),i+5:sx.QQ(rng.randrange(1,10),rng.randrange(1,10))}
    if seed%2:data[4]=data[1].copy()
    rep=data if fmt=='sparse' else [[data[i].get(j,sx.QQ.zero)for j in range(cols)]for i in range(rows)]
    a=DomainMatrix(rep,(rows,cols),sx.QQ);original=a.to_list();calls=[];native=sx._core._sdm_irref
    def spy(*args):
        result=native(*args);calls.append(result is not None);return result
    monkeypatch.setattr(sx._core,'_sdm_irref',spy)
    rref=importlib.import_module('symox.polys.matrices.rref')
    assert rref._dm_rref_choose_method(a,'auto')==('GJ','sparse')
    name='rref' if operation=='rref' else 'rref_den';options={} if name=='rref' else {'keep_domain':operation=='rref_den'}
    actual=getattr(a,name)(**options)
    assert calls and all(calls)
    ref=sp.polys.matrices.DomainMatrix([[sp.QQ(v.numerator.__int__(),v.denominator.__int__())for v in row]for row in original],(rows,cols),sp.QQ)
    expected=getattr(ref,name)(**options)
    got=actual[0];assert type(got) is DomainMatrix and got.rep.fmt==fmt and got.domain==sx.QQ
    assert actual[-1]==expected[-1]
    matrix=sp.Matrix([[sp.Rational(rational(v).numerator,rational(v).denominator)for v in row]for row in original])
    reduced=sp.Matrix([[sp.Rational(rational(v).numerator,rational(v).denominator)for v in row]for row in got.to_list()])
    assert reduced==expected[0].to_Matrix()
    assert matrix.col_join(reduced).rank()==matrix.rank()==len(actual[-1])
    for vector in reduced.nullspace():assert matrix*vector==sp.zeros(rows,1)
    if name=='rref_den':assert rational(actual[1])==rational(expected[1])
    with monkeypatch.context() as patch:
        patch.setattr(module,'_native_irref',lambda *args:None);old=getattr(a,name)(**options)
    assert actual==old and a.to_list()==original


@pytest.mark.parametrize('order',list(itertools.permutations(('stats','vector','holonomic.recurrence'))))
def test_official_imports_gaussian_wrappers_and_reloads_keep_native_qq(order):
    source='''
import importlib,sys
import symox as s
m=importlib.import_module('symox.polys.matrices.sdm')
a={0:{0:s.QQ(1,2),2:s.QQ(2,3)},1:{1:s.QQ(3,5),2:s.QQ(4,7)}}
assert m._native_irref(a) is not None
expected=m.sdm_irref(a)
for name in sys.argv[1:]:importlib.import_module('symox.'+name)
import symox.physics.units
for repeat in range(3):
 importlib.reload(m)
 assert m._native_irref(a) is not None and m.sdm_irref(a)==expected
'''
    subprocess.run([sys.executable,'-c',source,*order],check=True,env=os.environ.copy())


def test_signal_interruptions_propagate_without_mutating_original_sparse_rows():
    source='''
import signal
import symox as s
from symox.polys.matrices.sdm import sdm_irref
n=90
a={i:{j:s.QQ((i*13+j*17)%29+1,31+i+j)for j in range(n)}for i in range(n)}
original=[(id(row),[(id(j),id(v))for j,v in row.items()])for row in a.values()]
def interrupt(*args):raise KeyboardInterrupt('sparse signal')
signal.signal(signal.SIGALRM,interrupt);signal.setitimer(signal.ITIMER_REAL,0.005)
try:
 try:sdm_irref(a)
 except KeyboardInterrupt as error:assert str(error)=='sparse signal'
 else:raise AssertionError('signal not observed')
finally:signal.setitimer(signal.ITIMER_REAL,0)
assert original==[(id(row),[(id(j),id(v))for j,v in row.items()])for row in a.values()]
'''
    subprocess.run([sys.executable,'-c',source],check=True,env=os.environ.copy(),timeout=20)


@pytest.mark.parametrize('mode',['values','getattribute','instance_values'])
def test_stock_sdm_consumed_lookup_overrides_keep_old_calls(monkeypatch,mode):
    a=SDM(coefficients(),(50,6),sx.QQ);events=[]
    if mode=='values':
        def changed(self):events.append('values');return dict.values(self)
        monkeypatch.setattr(SDM,'values',changed,raising=False)
    elif mode=='getattribute':
        def changed(self,key):events.append(key);return dict.__getattribute__(self,key)
        monkeypatch.setattr(SDM,'__getattribute__',changed,raising=False)
    else:
        def changed():events.append('instance-values');return dict.values(a)
        a.values=changed
    actual,old,observed=fallback_matches(monkeypatch,a,events)
    assert observed
    monkeypatch.undo()
    assert signature(actual)==signature(old)


def test_unused_sdm_metadata_and_unrelated_class_fields_are_not_probed(monkeypatch):
    a=SDM(coefficients(),(50,6),sx.QQ);events=[]
    class Metadata:
        def __getattribute__(self,key):events.append(key);raise AssertionError('unused metadata')
    a.domain=Metadata();a.unrelated_payload=Metadata()
    monkeypatch.setattr(SDM,'unrelated_marker',Metadata(),raising=False)
    assert module._native_irref(a) is not None and events==[]
    actual=module.sdm_irref(a)
    with monkeypatch.context() as patch:
        patch.setattr(module,'_native_irref',lambda *args:None);old=module.sdm_irref(a)
    assert signature(actual)==signature(old) and events==[]


@pytest.mark.parametrize('mode',['builtin_before','pow_before','attribute_after'])
def test_reload_does_not_authenticate_callbacks_as_new_stock(mode):
    source='''
import importlib
import symox as s
m=importlib.import_module('symox.polys.matrices.sdm')
a={0:{0:s.QQ(1,2),2:s.QQ(2,3)},1:{1:s.QQ(3,5),2:s.QQ(4,7)}}
events=[]
'''
    if mode=='builtin_before':source+='''
original=min
def changed(*args,**kwargs):events.append('min');return original(*args,**kwargs)
m.min=changed
importlib.reload(m)
'''
    elif mode=='pow_before':source+='''
original=s._core.DomainElement.__pow__
def changed(*args,**kwargs):events.append('pow');return original(*args,**kwargs)
s._core.DomainElement.__pow__=changed
importlib.reload(m)
'''
    else:source+='''
importlib.reload(m)
def changed(self,name):events.append(name);raise AttributeError(name)
s._core.DomainElement.__getattr__=changed
'''
    source+='''
assert m._native_irref(a) is None and events==[]
result=m.sdm_irref(a)
assert events and result[1]==[0,1]
'''
    subprocess.run([sys.executable,'-c',source],check=True,env=os.environ.copy())


@pytest.mark.parametrize('name',['min','sorted','len'])
def test_first_import_does_not_authenticate_replaced_builtins(name):
    source='''
PRELOAD
import builtins
name=NAME
original=getattr(builtins,name)
events=[]
def changed(*args,**kwargs):
    events.append(name)
    return original(*args,**kwargs)
setattr(builtins,name,changed)
import importlib,symox as s
m=importlib.import_module('symox.polys.matrices.sdm')
REGISTER
a={0:{0:s.QQ(1),2:s.QQ(2)},1:{1:s.QQ(1),2:s.QQ(3)}}
events.clear()
assert m._native_irref(a) is None and events==[]
actual=m.sdm_irref(a)
observed=events[:]
m._native_irref=lambda *args:None
events.clear()
old=m.sdm_irref(a)
legacy=events[:]
setattr(builtins,name,original)
assert actual==old and observed==legacy and observed
'''.replace('NAME',repr(name))
    # Replacing len before importing dataclasses already prevents the existing
    # package from importing. Exercise the initial registration independently.
    source=source.replace('PRELOAD','import symox' if name=='len' else '')
    registration="p=importlib.import_module('symox.polys.matrices._sdm_irref_native'); p._PROTOCOL=None; m._native_irref=p.register_irref(m.SDM,m.__dict__)"
    source=source.replace('REGISTER',registration if name=='len' else '')
    subprocess.run([sys.executable,'-c',source],check=True,env=os.environ.copy())


@pytest.mark.parametrize('name',['set','dict','enumerate','defaultdict'])
def test_first_registration_rejects_heap_factories_with_copied_builtin_names(monkeypatch,name):
    import builtins,collections
    a=coefficients();events=[]
    original=collections.defaultdict if name=='defaultdict' else getattr(builtins,name)
    def new(cls,*args,**kwargs):
        events.append(name)
        return original.__new__(cls,*args,**kwargs)
    factory=type(name,(original,),{'__module__':original.__module__,'__new__':new})
    monkeypatch.setattr(module,name,factory,raising=False)
    monkeypatch.setattr(protocol,'_PROTOCOL',None)
    protocol.register_irref(SDM,module.__dict__)
    actual,old,observed=fallback_matches(monkeypatch,a,events)
    assert observed and signature(actual)==signature(old)


@pytest.mark.parametrize('name',['_ALL','_ANY'])
def test_selection_never_invokes_a_custom_helper_captured_before_registration(monkeypatch,name):
    a=coefficients();events=[];original=getattr(protocol,name)
    def callback(*args,**kwargs):events.append(name);return original(*args,**kwargs)
    monkeypatch.setattr(protocol,name,callback)
    monkeypatch.setattr(protocol,'_PROTOCOL',None)
    protocol.register_irref(SDM,module.__dict__)
    actual,old,observed=fallback_matches(monkeypatch,a,events)
    assert observed==[] and signature(actual)==signature(old)
