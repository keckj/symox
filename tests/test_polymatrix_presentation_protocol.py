"""Ground elimination consumes arithmetic, without presenting expressions."""
import importlib
import os
import subprocess
import sys

import pytest
import symox as sx

pm=importlib.import_module('symox.polys.polymatrix')
protocol=importlib.import_module('symox.polys._polymatrix_native')
x,y,z=sx.symbols('x y z')


def matrix(kind):
    if kind=='QQ':return pm.PolyMatrix([[1,2,3],[3,7,8]],x,ring=sx.QQ)
    k=sx.ZZ.frac_field(y) if kind=='ZRF' else (sx.QQ.frac_field(y,z) if kind=='multi' else sx.QQ.frac_field(y))
    a=y/(z+1) if kind=='multi' else y+1
    return pm.PolyMatrix([[a,1,2],[2*a,3,5]],x,ring=k)


@pytest.mark.parametrize('kind',['QQ','RF','multi','ZRF'])
@pytest.mark.parametrize('basis',[False,True])
@pytest.mark.parametrize('name',['__str__','__repr__'])
@pytest.mark.parametrize('mode',['count','raise'])
def test_native_and_complete_old_ground_never_call_presentation(monkeypatch,kind,basis,name,mode):
    a=matrix(kind);events=[];calls=[];native=sx._core._polymatrix_ground;original=getattr(sx._core.Expr,name)
    def spy(*args):
        result=native(*args);calls.append(result is not None);return result
    def presentation(self):
        events.append(name)
        if mode=='raise':raise RuntimeError('ground unexpectedly presented an expression')
        return original(self)
    monkeypatch.setattr(sx._core,'_polymatrix_ground',spy)
    with monkeypatch.context() as patch:
        patch.setattr(sx._core.Expr,name,presentation)
        assert pm._native_ground(a,basis) is not None
        actual=a.nullspace() if basis else a.rref()
        assert events==[]
    assert calls and all(calls)
    with monkeypatch.context() as patch:
        patch.setattr(pm,'_native_ground',lambda *args:None)
        patch.setattr(sx._core.Expr,name,presentation)
        expected=a.nullspace() if basis else a.rref()
        assert events==[]
    assert actual==expected


def test_official_units_printing_and_optional_imports_keep_native_after_reload():
    source='''
import importlib
import symox as s
from symox.polys.polymatrix import PolyMatrix,_native_ground
x,y=s.symbols('x y');a=PolyMatrix([[y,1,2],[2*y,3,5]],x,ring=s.QQ.frac_field(y))
before=_native_ground(a,True);assert before is not None
import symox.physics.units
import symox.stats,symox.vector,symox.holonomic.recurrence
r=importlib.import_module('symox.integrals.rationaltools')
for repeat in range(3):
 importlib.reload(r)
 assert _native_ground(a,True)==before
'''
    subprocess.run([sys.executable,'-c',source],check=True,env=os.environ.copy())


def test_class_snapshot_omits_print_functions_even_when_already_python(monkeypatch):
    events=[]
    def presentation(self):events.append(1);raise AssertionError('snapshot presentation')
    def replacement(self):events.append(2);raise AssertionError('changed presentation')
    monkeypatch.setattr(sx._core.Expr,'__str__',presentation)
    monkeypatch.setattr(sx._core.Expr,'__repr__',presentation)
    record=protocol._class(sx._core.Expr)
    assert all(state[0] is not presentation for state in record[3])
    monkeypatch.setattr(presentation,'__code__',replacement.__code__)
    assert protocol._class_unchanged(record,protocol._ROOT[4]) and events==[]


@pytest.mark.parametrize('mode',['count','raise'])
@pytest.mark.parametrize('context',['duplicate_gens','cancel_error','items_error'])
def test_live_error_presentation_is_retained_in_complete_fallback(monkeypatch,mode,context):
    class PresentationError(RuntimeError):pass
    class UserError(ValueError):pass
    name='__repr__' if context=='duplicate_gens' else '__str__'
    original=getattr(sx._core.Expr,name);events=[]
    a=pm.PolyMatrix([[1]],x,x,ring=sx.QQ) if context=='duplicate_gens' else matrix('QQ')
    if context=='cancel_error':
        def cancellation(expr,*args,**kwargs):raise UserError('cancel error: '+str(expr))
        monkeypatch.setattr(sx,'cancel',cancellation)
    elif context=='items_error':
        class Items(list):
            def __getitem__(self,key):raise UserError('items error: '+str(list.__getitem__(self,0)))
        a._items=Items(a._items)
    def presentation(self):
        events.append(name)
        if mode=='raise':raise PresentationError('actual error presentation')
        return original(self)
    def invoke():
        with monkeypatch.context() as patch:
            patch.setattr(sx._core.Expr,name,presentation)
            try:a.rref()
            except (PresentationError,UserError,sx.GeneratorsError) as error:
                return type(error),str(error)
        raise AssertionError('expected the established error')
    actual=invoke();seen=events[:];events.clear()
    with monkeypatch.context() as patch:
        patch.setattr(pm,'_native_ground',lambda *args:None)
        expected=invoke()
    assert actual==expected and events==seen and seen
    if mode=='raise':assert actual[0] is PresentationError


@pytest.mark.parametrize('name',['__truediv__','__getattribute__'])
def test_nonpresentation_arithmetic_and_attribute_changes_still_refuse(monkeypatch,name):
    a=matrix('RF');events=[];original=getattr(sx._core.Expr,name)
    def callback(*args):events.append(name);return original(*args)
    monkeypatch.setattr(sx._core.Expr,name,callback)
    assert pm._native_ground(a) is None and events==[]


@pytest.mark.parametrize('name',['__str__','__repr__'])
def test_removed_presentation_key_cannot_hide_an_added_attribute_hook(monkeypatch,name):
    a=matrix('QQ');events=[]
    def attribute(self,key):events.append(key);raise AttributeError(key)
    monkeypatch.delattr(sx._core.Expr,name)
    monkeypatch.setattr(sx._core.Expr,'__getattr__',attribute,raising=False)
    assert pm._native_ground(a) is None and events==[]
