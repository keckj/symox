"""One authenticated late host hook must not disable exact ground arithmetic."""
import importlib
import types
import sys
import os
import subprocess

import pytest
import symox as sx

from symox import _host

pm=importlib.import_module('symox.polys.polymatrix')
guard=importlib.import_module('symox.polys._polymatrix_native')


def install():
    class Host(sx.Symbol):
        @property
        def free_symbols(self):return set()
    return Host


def matrix(kind):
    x,y,z=sx.symbols('host_ground_x host_ground_y host_ground_z')
    field=sx.QQ if kind=='QQ' else sx.QQ.frac_field(y,z)
    entries=[[1,2,3],[3,7,8]] if kind=='QQ' else [[y/(z+1),1,2],[2*y/(z+1),3,5]]
    return pm.PolyMatrix(entries,x,ring=field)


@pytest.mark.parametrize('kind',['QQ','RF'])
@pytest.mark.parametrize('operation',['rref','nullspace'])
def test_official_host_hook_keeps_native_ground_and_full_public_result(monkeypatch,kind,operation):
    install();a=matrix(kind);calls=[];native=sx._core._polymatrix_ground
    def spy(*args):
        result=native(*args);calls.append(result is not None);return result
    monkeypatch.setattr(sx._core,'_polymatrix_ground',spy)
    before=tuple(a._items);host_calls=[];getter_code=_host._tree_free_symbols.__code__
    def profile(frame,event,arg):
        if event=='call' and frame.f_code is getter_code:host_calls.append(event)
    previous=sys.getprofile();sys.setprofile(profile)
    try:result=getattr(a,operation)()
    finally:sys.setprofile(previous)
    assert calls and all(calls) and host_calls==[]
    monkeypatch.setattr(pm,'_native_ground',lambda *args:None)
    sys.setprofile(profile)
    try:expected=getattr(a,operation)()
    finally:sys.setprofile(previous)
    assert expected==result and host_calls==[]
    assert tuple(a._items)==before and all(a is b for a,b in zip(a._items,before))


@pytest.mark.parametrize('kind',['QQ','RF'])
@pytest.mark.parametrize('operation',['rref','nullspace'])
@pytest.mark.parametrize('raises',[False,True])
def test_foreign_host_getter_is_not_called_by_candidate_or_complete_old_body(monkeypatch,kind,operation,raises):
    install();a=matrix(kind);calls=[]
    def callback(self):
        calls.append('free_symbols')
        if raises:raise RuntimeError('foreign host getter')
        return set()
    monkeypatch.setattr(_host.HostTree,'free_symbols',property(callback))
    assert pm._native_ground(a) is None and calls==[]
    result=getattr(a,operation)();observed=calls[:];calls.clear()
    monkeypatch.setattr(pm,'_native_ground',lambda *args:None)
    assert getattr(a,operation)()==result
    assert calls==observed==[]


@pytest.mark.parametrize('mode',['getter_code','getter_defaults','getter_alias','getter_foreign_globals',
    'installer_code','installer_defaults','installer_alias','property_setter','unrelated',
    'builtin_shadow','core_alias','module','module_class'])
def test_only_exact_authenticated_official_host_property_is_admitted(monkeypatch,mode):
    install();a=matrix('QQ');calls=[]
    def callback(*args):calls.append('callback');raise RuntimeError('foreign host protocol')
    if mode=='getter_code':monkeypatch.setattr(_host._tree_free_symbols,'__code__',(lambda self:None).__code__)
    elif mode=='getter_defaults':monkeypatch.setattr(_host._tree_free_symbols,'__defaults__',(callback,))
    elif mode=='getter_alias':monkeypatch.setattr(_host,'_tree_free_symbols',callback)
    elif mode=='getter_foreign_globals':
        f=_host._tree_free_symbols;foreign=types.FunctionType(f.__code__,dict(f.__globals__))
        monkeypatch.setattr(_host.HostTree,'free_symbols',property(foreign))
    elif mode=='installer_code':monkeypatch.setattr(_host._install_free_symbols,'__code__',(lambda:None).__code__)
    elif mode=='installer_defaults':monkeypatch.setattr(_host._install_free_symbols,'__defaults__',(callback,))
    elif mode=='installer_alias':monkeypatch.setattr(_host,'_install_free_symbols',callback)
    elif mode=='property_setter':monkeypatch.setattr(_host.HostTree,'free_symbols',property(_host._tree_free_symbols,callback))
    elif mode=='unrelated':monkeypatch.setattr(_host.HostTree,'unrelated_callback',callback,raising=False)
    elif mode=='builtin_shadow':monkeypatch.setattr(_host,'property',callback,raising=False)
    elif mode=='core_alias':monkeypatch.setattr(_host,'_core',object())
    elif mode=='module_class':
        class Hooks(types.ModuleType):
            def __getattribute__(self,name):
                calls.append(name)
                if name=='__dict__':raise RuntimeError('foreign module lookup')
                return super().__getattribute__(name)
        monkeypatch.setattr(_host,'__class__',Hooks)
    else:
        fake=types.ModuleType(_host.__name__);fake.__dict__.update(_host.__dict__)
        monkeypatch.setitem(__import__('sys').modules,_host.__name__,fake)
    assert pm._native_ground(a) is None and calls==[]


@pytest.mark.parametrize('operation',['rref','nullspace'])
def test_real_host_coefficient_keeps_full_failure_and_its_free_symbol_protocol(monkeypatch,operation):
    calls=[]
    class Rate(sx.Symbol):
        @property
        def free_symbols(self):
            calls.append('host coefficient');raise RuntimeError('host coefficient free_symbols')
    x=sx.Symbol('host_ground_x');rate=Rate('opaque_host_rate')
    a=pm.PolyMatrix([[rate+1,1]],x,ring=sx.QQ)
    assert pm._native_ground(a) is None and calls==[]
    with pytest.raises(sx.CoercionFailed) as actual:getattr(a,operation)()
    observed=calls[:];calls.clear()
    monkeypatch.setattr(pm,'_native_ground',lambda *args:None)
    with pytest.raises(type(actual.value)) as expected:getattr(a,operation)()
    assert str(actual.value)==str(expected.value) and calls==observed==[]
    # The coefficient's real host protocol is still active and observable.
    with pytest.raises(RuntimeError,match='host coefficient free_symbols'):
        a._items[0].free_symbols
    assert calls==['host coefficient']


def test_shadowed_property_before_import_declines_without_descriptor_callbacks():
    code="""
import builtins
calls=[]
class ForeignProperty(builtins.property):
    def __getattribute__(self,name):
        if name in ('fget','fset','fdel'):
            calls.append(name)
        return super().__getattribute__(name)
builtins.property=ForeignProperty
import symox as s
from symox import _host
from symox.polys import polymatrix as pm
class Host(s.Symbol):
    @property
    def free_symbols(self):return set()
x=s.Symbol('foreign_property_x');a=pm.PolyMatrix([[1,2,3],[3,7,8]],x,ring=s.QQ)
calls.clear()
selected=pm._native_ground(a)
assert selected is None and calls==[],calls
result=a.rref();observed=calls[:];calls.clear()
pm._native_ground=lambda *args:None
assert a.rref()==result
assert calls==observed==[],(calls,observed)
"""
    subprocess.run([sys.executable,'-c',code],check=True,env=os.environ.copy())
