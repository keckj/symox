"""Sparse exact Berkowitz: complete coefficients, storage and protocol parity."""
import abc
import builtins
import importlib
import json
import os
import random
import subprocess
import sys
from collections import Counter

import pytest
import sympy as sp
import symox as sx
from sympy.polys.matrices.sdm import sdm_berk as reference
from symox.polys.matrices.domainmatrix import DomainMatrix
from symox.polys.matrices.sdm import SDM

sparse=importlib.import_module('symox.polys.matrices.sdm')
methods=importlib.import_module('symox.polys.domains._methods')


def values(raw,k):
    return {i:{j:k(int(sp.numer(v)),int(sp.denom(v))) if k is sx.QQ else k(int(v))
               for j,v in row.items()} for i,row in raw.items()}


def warm(k):
    # Keep cold semantic work on its established path. Native selection is
    # optional until the exact coefficient class has a valid Basic=False cache.
    assert k.sum([k.one])==k.one


def alias_partition(result,source):
    classes={}
    for v in source:classes.setdefault(id(v),len(classes))
    return [classes.setdefault(id(v),len(classes)) for v in result.values()]


@pytest.mark.parametrize('k',[sx.ZZ,sx.QQ])
@pytest.mark.parametrize('seed',range(36))
def test_entire_sparse_coefficients_reference_cayley_hamilton_and_aliases(k,seed):
    rng=random.Random(719231+seed);n=seed%9
    raw={i:{j:sp.Rational(rng.randrange(-9,10),rng.randrange(1,11) if k is sx.QQ else 1)
            for j in range(n) if rng.randrange(4)==0} for i in range(n)}
    if seed%7==0:raw={}
    a=values(raw,k)
    if n>1 and seed%5==0:a[n-1]=a.get(0,{})
    before={i:dict(row) for i,row in a.items()}
    row_ids={i:id(row) for i,row in a.items()}
    cell_ids={i:{j:id(v) for j,v in row.items()} for i,row in a.items()}
    actual=sx._core._sdm_berk(a,n,k)
    expected=reference(a,n,k)
    assert actual is not None and type(actual) is dict and actual==expected
    assert list(actual)==sorted(expected) and all(bool(v) for v in actual.values())
    assert all(type(v) is type(k.one) and v.parent()==k for v in actual.values())
    source=[v for row in a.values() for v in row.values()]
    assert alias_partition(actual,source)==alias_partition(expected,source)
    assert a==before and row_ids=={i:id(row) for i,row in a.items()}
    assert cell_ids=={i:{j:id(v) for j,v in row.items()} for i,row in a.items()}
    mat=sp.zeros(n)
    for i,row in a.items():
        for j,v in row.items():mat[i,j]=sp.Rational(str(v))
    coefficients=[sp.Rational(str(actual.get(i,k.zero))) for i in range(n+1)]
    assert coefficients==mat.charpoly().all_coeffs()
    evaluated=sp.zeros(n)
    for coefficient in coefficients:evaluated=evaluated*mat+coefficient*sp.eye(n)
    assert evaluated==sp.zeros(n)


@pytest.mark.parametrize('k',[sx.ZZ,sx.QQ])
def test_large_exact_coefficients_and_coprime_denominators(k):
    raw={0:{0:2**211+3,2:sp.Rational(2,17) if k is sx.QQ else 17},
         1:{0:-3**81,1:sp.Rational(7**61,11**31) if k is sx.QQ else 7**61},
         2:{1:sp.Rational(1,5**47) if k is sx.QQ else 5**47,2:3}}
    a=values(raw,k)
    assert sx._core._sdm_berk(a,3,k)==reference(a,3,k)


@pytest.mark.parametrize('entry',['SDM','charpoly_berk','charpoly'])
@pytest.mark.parametrize('k',[sx.ZZ,sx.QQ])
def test_real_public_sparse_cycle_uses_native_without_dense_conversion(monkeypatch,entry,k):
    n=7
    a=SDM({i:{(i+1)%n:k(i+2)} for i in range(n)},(n,n),k)
    dm=DomainMatrix.from_rep(a);before=a.copy();warm(k)
    calls=[];native=sx._core._sdm_berk
    def spy(*args):
        result=native(*args);calls.append(result is not None);return result
    monkeypatch.setattr(sx._core,'_sdm_berk',spy)
    def dense(*args,**kwargs):raise AssertionError('unexpected dense conversion')
    monkeypatch.setattr(SDM,'to_ddm',dense)
    actual=a.charpoly() if entry=='SDM' else getattr(dm,entry)()
    expected=[k.one]+[k.zero]*(n-1)+[-k(sp.prod(range(2,n+2)))]
    assert actual==expected and calls and all(calls) and a==before


@pytest.mark.parametrize('k',[sx.ZZ,sx.QQ])
def test_isolated_zero_vertices_keep_sparse_leading_indices_for_huge_dimensions(k):
    n=sys.maxsize
    a={n-9:{n-9:k(2),n-2:k(3)},n-2:{n-9:k(5),n-2:k(7)}}
    # Simultaneous permutation gives a2x2 block direct-summed with zeros:
    # lambda**(n-2)*(lambda**2-9*lambda-1). No n-vector is needed here.
    assert sx._core._sdm_berk(a,n,k)=={0:k.one,1:k(-9),2:k(-1)}
    assert sx._core._sdm_berk({},n,k)=={0:k.one}


def test_zero_dimension_preserves_the_old_ignored_storage_convention():
    assert sx._core._sdm_berk({object():object()},0,sx.ZZ)=={0:sx.ZZ.one}


@pytest.mark.parametrize('hook',['one','zero','sum','convert','neg','mul','add','bool',
                               'dot','matvec','recursion','dict_constructor','range',
                               'sdm_items','element_attributes','basic_check','convert_helper'])
def test_custom_protocols_keep_the_complete_existing_callback_trace(monkeypatch,hook):
    k=sx.QQ;a=SDM(values({0:{0:1,1:2,2:3},1:{0:4,1:5,2:6},2:{0:7,1:8,2:10}},k),(3,3),k)
    original=sparse.sdm_berk;events=[];warm(k)
    with monkeypatch.context() as patch:
        if hook in ('one','zero'):
            cls=type(k);descriptor=cls.__dict__[hook]
            def get(self):events.append(hook);return descriptor.__get__(self,cls)
            patch.setattr(cls,hook,property(get))
        elif hook in ('sum','convert'):
            cls=type(k);old=getattr(cls,hook)
            def callback(self,*args,**kwargs):events.append(hook);return old(self,*args,**kwargs)
            patch.setattr(cls,hook,callback)
        elif hook in ('neg','mul','add','bool'):
            cls=type(k.one);name={'neg':'__neg__','mul':'__mul__','add':'__add__','bool':'__bool__'}[hook];old=getattr(cls,name)
            def callback(self,*args):events.append(hook);return old(self,*args)
            patch.setattr(cls,name,callback)
        elif hook in ('dot','matvec','recursion'):
            name={'dot':'sdm_dotvec','matvec':'sdm_matvecmul','recursion':'sdm_berk'}[hook];old=getattr(sparse,name)
            def callback(*args):events.append(hook);return old(*args)
            patch.setattr(sparse,name,callback)
        elif hook in ('dict_constructor','range'):
            name='dict' if hook=='dict_constructor' else 'range';old=getattr(builtins,name)
            def callback(*args,**kwargs):events.append(hook);return old(*args,**kwargs)
            patch.setattr(sparse,name,callback,raising=False)
        elif hook=='sdm_items':
            old=SDM.items
            def callback(self):events.append(hook);return old(self)
            patch.setattr(SDM,'items',callback)
        elif hook=='element_attributes':
            cls=type(k.one)
            def get(self,name):
                if name=='__class__':events.append(hook)
                return object.__getattribute__(self,name)
            patch.setattr(cls,'__getattribute__',get)
        elif hook=='basic_check':
            cls=type(sx.Basic);old=cls.__instancecheck__
            def callback(self,value):events.append(hook);return old(self,value)
            patch.setattr(cls,'__instancecheck__',callback)
        else:
            old=methods._inexact_value
            def callback(value):events.append(hook);return old(value)
            patch.setattr(methods,'_inexact_value',callback)
        actual=original(a,3,k);trace=events[:];events.clear()
        with monkeypatch.context() as disable:
            disable.setattr(sparse,'_sdm_stock_protocol',lambda *args:False)
            disable.setattr(sx._core,'_sdm_berk',lambda *args:None)
            expected=original(a,3,k)
        assert actual==expected and events==trace and trace


@pytest.mark.parametrize('storage',['outer','row','integers','mixed','domain_mismatch','bool_key',
                                  'negative_key','out_of_range','noninteger_dimension','bool_dimension'])
def test_noncanonical_representations_keep_old_value_or_error(monkeypatch,storage):
    events=[]
    class Custom(dict):
        def items(self):events.append('items');return super().items()
    def make():
        a=values({0:{0:1,1:2},1:{0:3,1:4}},sx.ZZ);n=2
        if storage=='outer':a=Custom(a)
        elif storage=='row':a[1]=Custom(a[1])
        elif storage=='integers':a={0:{0:1,1:2},1:{0:3,1:4}}
        elif storage=='mixed':a[1][0]=3
        elif storage=='domain_mismatch':a[1][0]=sx.QQ(1,2)
        elif storage=='bool_key':a={False:a[0],True:a[1]}
        elif storage=='negative_key':a[-1]=a.pop(1)
        elif storage=='out_of_range':a[0][2]=sx.ZZ.one
        elif storage=='noninteger_dimension':n=2.0
        elif storage=='bool_dimension':n=True
        return a,n
    def run():
        a,n=make()
        try:
            result=sparse.sdm_berk(a,n,sx.ZZ)
            return ('value',[(i,str(v),type(v).__name__) for i,v in result.items()])
        except Exception as error:return (type(error),str(error))
    actual=run();trace=events[:];events.clear()
    with monkeypatch.context() as patch:
        patch.setattr(sparse,'_sdm_stock_protocol',lambda *args:False)
        patch.setattr(sx._core,'_sdm_berk',lambda *args:None)
        expected=run()
    assert actual==expected and trace==events


@pytest.mark.parametrize('domain',[sx.GF(7),sx.ZZ.poly_ring(sx.Symbol('t')),sx.RR,sx.CC,sx.EXRAW])
def test_other_domains_retain_the_complete_algorithm(monkeypatch,domain):
    a={0:{0:domain(1),1:domain(2)},1:{0:domain(3),1:domain(4)}}
    actual=sparse.sdm_berk(a,2,domain)
    with monkeypatch.context() as patch:
        patch.setattr(sparse,'_sdm_stock_protocol',lambda *args:False)
        patch.setattr(sx._core,'_sdm_berk',lambda *args:None)
        expected=sparse.sdm_berk(a,2,domain)
    assert actual==expected and [(i,str(v)) for i,v in actual.items()]==[(i,str(v)) for i,v in expected.items()]


def test_invalidated_abc_cache_uses_full_fallback_then_stock_native(monkeypatch):
    a=values({0:{0:1,1:2},1:{0:3,1:4}},sx.QQ);warm(sx.QQ)
    class Temporary(abc.ABC):pass
    class Registered:pass
    Temporary.register(Registered)
    assert not sparse._sdm_stock_protocol(a,sx.QQ)
    first=sparse.sdm_berk(a,2,sx.QQ)
    assert sparse._sdm_stock_protocol(a,sx.QQ)
    calls=[];native=sx._core._sdm_berk
    def spy(*args):result=native(*args);calls.append(result is not None);return result
    monkeypatch.setattr(sx._core,'_sdm_berk',spy)
    assert sparse.sdm_berk(a,2,sx.QQ)==first and calls==[True]


def test_stock_public_route_does_not_delegate_matrix_arithmetic_to_python():
    k=sx.QQ;a=SDM(values({0:{0:1,1:2},1:{0:3,1:4}},k),(2,2),k);warm(k)
    forbidden={sparse.sdm_dotvec.__code__,sparse.sdm_matvecmul.__code__,methods.sum_.__code__,methods.convert.__code__}
    calls=[]
    def profile(frame,event,arg):
        if event=='call' and frame.f_code in forbidden:calls.append(frame.f_code.co_name)
    old=sys.getprofile()
    try:
        sys.setprofile(profile)
        result=a.charpoly()
    finally:sys.setprofile(old)
    assert result==[k.one,k(-5),k(-2)] and not calls


@pytest.mark.parametrize('hook',['code','native_convert','sum_code','basic_cell'])
def test_changed_code_and_delegates_decline_without_replacing_old_errors(monkeypatch,hook):
    k=sx.QQ;a=values({0:{0:1,1:2},1:{0:3,1:4}},k);warm(k);events=[]
    if hook=='basic_cell':
        # This cell is shared by methods of the metaclass. Keep its mutation
        # in a separate interpreter so even a failure cannot pollute pytest.
        script='''
import symox as s,importlib,json
m=importlib.import_module('symox.polys.matrices.sdm');b=importlib.import_module('symox._basic')
k=s.QQ;a={0:{0:k(1),1:k(2)},1:{0:k(3),1:k(4)}};k.sum([k.one])
cell=b._BasicMeta.__instancecheck__.__closure__[0];original=cell.cell_contents
try:
    del cell.cell_contents
    refused=not m._sdm_stock_protocol(a,k)
    try:m.sdm_berk(a,2,k)
    except Exception as e:actual=(type(e).__name__,str(e))
    native=s._core._sdm_berk;s._core._sdm_berk=lambda *args:None
    try:m.sdm_berk(a,2,k)
    except Exception as e:expected=(type(e).__name__,str(e))
finally:cell.cell_contents=original
assert refused and actual==expected
print(json.dumps(actual))
'''
        result=subprocess.run([sys.executable,'-c',script],capture_output=True,text=True,timeout=10)
        assert result.returncode==0,result.stdout+result.stderr
        assert 'RuntimeError' in result.stdout
        return
    with monkeypatch.context() as patch:
        if hook=='native_convert':
            old=methods._native_convert
            def convert(*args):events.append('native_convert');return old(*args)
            patch.setattr(methods,'_native_convert',convert)
        else:
            function=methods._inexact_value if hook=='code' else methods.sum_
            patch.setitem(function.__globals__,'_sdm_test_events',events)
            if hook=='code':
                def replacement(value):
                    _sdm_test_events.append('code')
                    return None
            else:
                def replacement(self,args):
                    _sdm_test_events.append('sum_code')
                    result=self.zero
                    for value in args:result+=self.convert(value)
                    return result
            patch.setattr(function,'__code__',replacement.__code__)
        assert not sparse._sdm_stock_protocol(a,k)
        actual=sparse.sdm_berk(a,2,k);trace=events[:];events.clear()
        with monkeypatch.context() as disable:
            disable.setattr(sparse,'_sdm_stock_protocol',lambda *args:False)
            disable.setattr(sx._core,'_sdm_berk',lambda *args:None)
            expected=sparse.sdm_berk(a,2,k)
        assert actual==expected and events==trace and trace


def test_module_reload_keeps_a_conservative_complete_fallback():
    script='''
import symox as s,importlib
m=importlib.import_module('symox.polys.matrices.sdm');k=s.QQ
a={0:{0:k(1),1:k(2)},1:{0:k(3),1:k(4)}};k.sum([k.one])
old=m.sdm_berk
m=importlib.reload(m)
assert not m._sdm_stock_protocol(a,k)
assert m.sdm_berk(a,2,k)==old(a,2,k)=={0:k.one,1:k(-5),2:k(-2)}
'''
    result=subprocess.run([sys.executable,'-c',script],capture_output=True,text=True,timeout=10)
    assert result.returncode==0,result.stdout+result.stderr


def test_native_interruptions_preserve_input_and_never_become_refusal():
    script='''
import symox as s,signal
n=256;k=s.ZZ;v=k(2**1024+1)
a={i:{(i+1)%n:v} for i in range(n)};before={i:dict(row) for i,row in a.items()}
class Stop(Exception):pass
def handler(*args):raise Stop('interrupt-token')
signal.signal(signal.SIGINT,handler)
try:
    print('ready',flush=True)
    s._core._sdm_berk(a,n,k)
except Stop as error:
    assert str(error)=='interrupt-token' and a==before
    print('interrupted',flush=True)
else:raise AssertionError('native computation ignored interrupt')
'''
    process=subprocess.Popen([sys.executable,'-u','-c',script],stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    try:
        assert process.stdout.readline().strip()=='ready'
        process.send_signal(__import__('signal').SIGINT)
        out,err=process.communicate(timeout=15)
        assert process.returncode==0 and 'interrupted' in out,out+err
    finally:
        if process.poll() is None:process.kill();process.wait()


def test_stock_sparse_path_works_with_reference_imports_blocked():
    script='''
import sys
class NoReference:
    def find_spec(self,fullname,path=None,target=None):
        if fullname=='sympy' or fullname.startswith('sympy.'):raise AssertionError(fullname)
sys.meta_path.insert(0,NoReference())
import symox as s,importlib
m=importlib.import_module('symox.polys.matrices.sdm');k=s.QQ
n=7;a=m.SDM({i:{(i+1)%n:k(i+2)} for i in range(n)},(n,n),k);k.sum([k.one])
assert m._sdm_stock_protocol(a,k)
assert a.charpoly()==[k.one]+[k.zero]*6+[k(-40320)]
assert not any(name=='sympy' or name.startswith('sympy.') for name in sys.modules)
'''
    result=subprocess.run([sys.executable,'-c',script],capture_output=True,text=True,timeout=10)
    assert result.returncode==0,result.stdout+result.stderr


def test_fresh_reimport_does_not_bless_a_custom_callable_without_code():
    script='''
import symox as s,sys,importlib
m=importlib.import_module('symox.polys.matrices.sdm');k=s.QQ
old=type(k).sum;events=[]
class Callback:
    def __get__(self,instance,owner):
        def call(args):events.append('sum');return old(instance,args)
        return call
type(k).sum=Callback()
del sys.modules['symox.polys.matrices.sdm']
m=importlib.import_module('symox.polys.matrices.sdm')
a={0:{0:k(1),1:k(2)},1:{0:k(3),1:k(4)}}
assert not m._sdm_stock_protocol(a,k)
assert m.sdm_berk(a,2,k)=={0:k.one,1:k(-5),2:k(-2)} and events
'''
    result=subprocess.run([sys.executable,'-c',script],capture_output=True,text=True,timeout=10)
    assert result.returncode==0,result.stdout+result.stderr


@pytest.mark.parametrize('hook',['sum','convert','zero'])
@pytest.mark.parametrize('failure',[ValueError,KeyboardInterrupt])
def test_callback_errors_propagate_at_the_established_point(monkeypatch,hook,failure):
    k=sx.QQ;a=values({0:{0:1,1:2},1:{0:3,1:4}},k);warm(k);events=[]
    cls=type(k)
    def callback(*args,**kwargs):events.append(hook);raise failure('callback-token')
    if hook=='zero':callback=property(callback)
    with monkeypatch.context() as patch:
        patch.setattr(cls,hook,callback)
        def run():
            try:sparse.sdm_berk(a,2,k)
            except BaseException as error:return type(error),str(error)
            raise AssertionError('callback was bypassed')
        actual=run();trace=events[:];events.clear()
        with monkeypatch.context() as disable:
            disable.setattr(sparse,'_sdm_stock_protocol',lambda *args:False)
            disable.setattr(sx._core,'_sdm_berk',lambda *args:None)
            expected=run()
        assert actual==expected==(failure,'callback-token') and events==trace==[hook]


@pytest.mark.parametrize('target',['basic','abc'])
def test_metaclass_descriptor_callbacks_are_not_probed_by_selection(monkeypatch,target):
    basic=importlib.import_module('symox._basic')
    cls=basic._BasicMeta if target=='basic' else abc.ABCMeta
    old=cls.__dict__['__instancecheck__']
    k=sx.QQ;a=values({0:{0:1,1:2},1:{0:3,1:4}},k);warm(k);events=[]
    class Descriptor:
        def __get__(self,instance,owner):
            events.append(('bind',instance is None))
            return lambda value:old(instance,value)
    with monkeypatch.context() as patch:
        patch.setattr(cls,'__instancecheck__',Descriptor())
        assert not sparse._sdm_stock_protocol(a,k)
        assert not events
        actual=sparse.sdm_berk(a,2,k);trace=events[:];events.clear()
        with monkeypatch.context() as disable:
            disable.setattr(sparse,'_sdm_stock_protocol',lambda *args:False)
            disable.setattr(sx._core,'_sdm_berk',lambda *args:None)
            expected=sparse.sdm_berk(a,2,k)
        assert actual==expected and events==trace and trace
        assert all(bound is False for _,bound in trace)


@pytest.mark.parametrize('target',['Domain','DomainElement','native','sum_alias','basic_alias','abc_dump'])
def test_missing_optional_aliases_keep_the_established_value_or_error(monkeypatch,target):
    basic=importlib.import_module('symox._basic');abc_native=importlib.import_module('_abc')
    k=sx.QQ;a=values({0:{0:1,1:2},1:{0:3,1:4}},k);warm(k)
    owner,name={'Domain':(sx._core,'Domain'),'DomainElement':(sx._core,'DomainElement'),
                'native':(sx._core,'_sdm_berk'),'sum_alias':(methods,'sum_'),
                'basic_alias':(basic,'Basic'),'abc_dump':(abc_native,'_get_dump')}[target]
    def run():
        try:return ('value',sparse.sdm_berk(a,2,k))
        except Exception as error:return (type(error),str(error))
    with monkeypatch.context() as patch:
        patch.delattr(owner,name)
        actual=run()
        with monkeypatch.context() as disable:
            disable.setattr(sx._core,'_sdm_berk',lambda *args:None,raising=False)
            expected=run()
        assert actual==expected


@pytest.mark.parametrize('method',['items','get'])
def test_instance_method_overrides_keep_existing_callbacks(monkeypatch,method):
    k=sx.QQ;a=SDM({0:{0:k(3)}},(1,1),k);warm(k);events=[];old=getattr(SDM,method)
    def callback(*args):events.append(method);return old(a,*args)
    a.__dict__[method]=callback
    assert not sparse._sdm_stock_protocol(a,k)
    actual=sparse.sdm_berk(a,1,k);trace=events[:];events.clear()
    with monkeypatch.context() as patch:
        patch.setattr(sparse,'_sdm_stock_protocol',lambda *args:False)
        patch.setattr(sx._core,'_sdm_berk',lambda *args:None)
        patch.setattr(sx._core,'_sdm_berk_combine',lambda *args:None)
        expected=sparse.sdm_berk(a,1,k)
    assert actual==expected and trace==events and trace


@pytest.mark.parametrize('shape',['dense','cycle','diagonal'])
def test_real_first_call_can_enter_native_at_the_original_conversion_boundary(monkeypatch,shape):
    # Invalidate only an unrelated ABC registry, without altering Basic's
    # membership rules. This recreates the cold protocol state honestly.
    k=sx.QQ;n=7
    if shape=='dense':raw={i:{j:(i*3+j)%9+1 for j in range(n)} for i in range(n)}
    elif shape=='cycle':raw={i:{(i+1)%n:i+2} for i in range(n)}
    else:raw={i:{i:i+1} for i in range(n)}
    a=SDM(values(raw,k),(n,n),k)
    class Temporary(abc.ABC):pass
    class Registered:pass
    Temporary.register(Registered)
    assert not sparse._sdm_stock_protocol(a,k)
    calls=[];native=sx._core._sdm_berk;combine=sx._core._sdm_berk_combine
    def spy(*args):result=native(*args);calls.append(('matrix',result is not None));return result
    def finish(*args):result=combine(*args);calls.append(('combine',result is not None));return result
    monkeypatch.setattr(sx._core,'_sdm_berk',spy);monkeypatch.setattr(sx._core,'_sdm_berk_combine',finish)
    actual=a.charpoly()
    assert any(ok for _,ok in calls)
    if shape!='diagonal':assert any(kind=='matrix' and ok for kind,ok in calls)
    else:assert any(kind=='combine' and ok for kind,ok in calls)
    mat=sp.zeros(n)
    for i,row in raw.items():
        for j,v in row.items():mat[i,j]=v
    assert [sp.Rational(str(v)) for v in actual]==mat.charpoly().all_coeffs()


@pytest.mark.parametrize('mutate',['input','prefix_result'])
def test_native_continuation_consumes_actual_prefix_results_without_replay(monkeypatch,mutate):
    k=sx.QQ;raw={0:{0:1,1:2,2:3},1:{0:4,1:5,2:6},2:{0:7,1:8,2:10}}
    cls=type(k);original_sum=cls.sum;native=sx._core._sdm_berk;combine=sx._core._sdm_berk_combine
    def run(enable):
        a=SDM(values(raw,k),(3,3),k);events=[];routes=[]
        def first_sum(self,args):
            events.append('first_sum');cls.sum=original_sum
            result=original_sum(self,args)
            if mutate=='input':a[2][2]=k(123)
            else:result+=k(1)
            return result
        with monkeypatch.context() as patch:
            patch.setattr(cls,'sum',first_sum)
            def matrix(*args):result=native(*args);routes.append(result is not None);return result
            def finish(*args):result=combine(*args);routes.append(result is not None);return result
            patch.setattr(sx._core,'_sdm_berk',matrix if enable else lambda *args:None)
            patch.setattr(sx._core,'_sdm_berk_combine',finish if enable else lambda *args:None)
            result=a.charpoly()
        return result,events,routes,a[2][2]
    actual,trace,routes,last=run(True);expected,reference_trace,_,reference_last=run(False)
    assert actual==expected and trace==reference_trace==['first_sum']
    assert last==reference_last and any(routes)


def test_empty_legacy_coefficient_map_keeps_the_original_explicit_error():
    assert sx._core._sdm_berk_combine({},[sx.QQ.one],2,sx.QQ) is None
    assert sx._core._sdm_berk_combine({0:sx.QQ.zero},[sx.QQ.one],2,sx.QQ)=={}


def test_nonstandard_instance_dictionary_keys_are_not_compared_by_selection(monkeypatch):
    k=sx.QQ;a=SDM(values({0:{0:1,1:2},1:{0:3,1:4},2:{2:5}},k),(3,3),k);warm(k);events=[]
    class Key:
        def __hash__(self):return hash('get')
        def __eq__(self,other):events.append(other);return False
    a.__dict__[Key()]=1
    actual=sparse.sdm_berk(a,3,k);trace=events[:];events.clear()
    with monkeypatch.context() as patch:
        patch.setattr(sparse,'_sdm_stock_protocol',lambda *args:False)
        patch.setattr(sx._core,'_sdm_berk',lambda *args:None)
        patch.setattr(sx._core,'_sdm_berk_combine',lambda *args:None)
        expected=sparse.sdm_berk(a,3,k)
    assert actual==expected and events==trace==[]


@pytest.mark.parametrize('k',[sx.ZZ,sx.QQ])
@pytest.mark.parametrize('route',['whole','prefix'])
@pytest.mark.parametrize('pattern',['leading','gapped','stored_zero','nilpotent'])
def test_zero_vertices_and_stored_zeros_on_whole_and_prefix_public_paths(monkeypatch,k,route,pattern):
    n=8
    raw={'leading':{3:{3:2,4:3},4:{3:5,4:7}},
         'gapped':{1:{1:2,6:3},6:{1:5,6:7}},
         'stored_zero':{0:{0:0,7:0},2:{2:2,5:3},5:{2:5,5:7},7:{}},
         'nilpotent':{2:{4:2},4:{7:3},7:{7:0}}}[pattern]
    a=SDM(values(raw,k),(n,n),k);before={i:dict(row) for i,row in a.items()};warm(k)
    calls=[];native=sx._core._sdm_berk;combine=sx._core._sdm_berk_combine
    def matrix(*args):
        result=native(*args) if route=='whole' else None
        calls.append(('matrix',result is not None));return result
    def finish(*args):
        result=combine(*args);calls.append(('combine',result is not None));return result
    monkeypatch.setattr(sx._core,'_sdm_berk',matrix)
    monkeypatch.setattr(sx._core,'_sdm_berk_combine',finish)
    actual=a.charpoly()
    expected_map=reference(before,n,k)
    assert actual==[expected_map.get(i,k.zero) for i in range(n+1)]
    assert a==before and any(kind==('matrix' if route=='whole' else 'combine') and ok for kind,ok in calls)
    if pattern=='nilpotent':assert actual==[k.one]+[k.zero]*n
    else:assert actual==[k.one,k(-9),k(-1)]+[k.zero]*(n-2)


@pytest.mark.parametrize('route',['whole','prefix'])
@pytest.mark.parametrize('location',['row_negative','column_negative','row_large','column_large'])
@pytest.mark.parametrize('zero',[False,True])
def test_invalid_indices_are_not_discarded_even_for_explicit_zeros(monkeypatch,route,location,zero):
    k=sx.QQ;n=3;a=values({0:{0:2,1:3},1:{0:5,1:7}},k)
    bad=-1 if location.endswith('negative') else n
    if location.startswith('row'):a[bad]={1:k.zero if zero else k.one}
    else:a.setdefault(1,{})[bad]=k.zero if zero else k.one
    assert sx._core._sdm_berk(a,n,k) is None
    before={i:dict(row) for i,row in a.items()};warm(k)
    native=sx._core._sdm_berk
    if route=='prefix':monkeypatch.setattr(sx._core,'_sdm_berk',lambda *args:None)
    def run():
        try:return ('value',sparse.sdm_berk(a,n,k))
        except Exception as error:return (type(error),str(error))
    actual=run()
    with monkeypatch.context() as patch:
        patch.setattr(sparse,'_sdm_stock_protocol',lambda *args:False)
        expected=run()
    assert actual==expected and a==before


@pytest.mark.parametrize('key',[-1,3,True,1.5])
def test_prefix_combiner_declines_invalid_coefficient_indices(key):
    assert sx._core._sdm_berk_combine({0:sx.ZZ.one,key:sx.ZZ.zero},[sx.ZZ.one],2,sx.ZZ) is None


def test_first_public_calls_in_a_fresh_process_report_the_actual_routes():
    script='''
import symox as s,importlib,json
from symox.polys.matrices.sdm import SDM
from symox.polys.matrices.domainmatrix import DomainMatrix
m=importlib.import_module('symox.polys.matrices.sdm')
# No warm-up, ABC manipulation or preliminary coefficient conversion.
n=7;k=s.QQ;a=SDM({i:{(i+1)%n:k(i+2)} for i in range(n)},(n,n),k)
routes=[];native=s._core._sdm_berk;combine=s._core._sdm_berk_combine
initial=m._sdm_stock_protocol(a,k)
def matrix(*args):
 r=native(*args);routes.append(['matrix',r is not None]);return r
def finish(*args):
 r=combine(*args);routes.append(['combine',r is not None]);return r
s._core._sdm_berk=matrix;s._core._sdm_berk_combine=finish
expected=[k.one]+[k.zero]*6+[-k(40320)]
assert a.charpoly()==expected and any(ok for _,ok in routes)
first=list(routes);routes.clear()
assert DomainMatrix.from_rep(a).charpoly()==expected and any(ok for _,ok in routes)
print(json.dumps({'initial':initial,'sdm':first,'domainmatrix':routes}))
'''
    result=subprocess.run([sys.executable,'-c',script],text=True,capture_output=True,check=True,timeout=30)
    routes=json.loads(result.stdout)
    assert routes['sdm'] and routes['domainmatrix']
