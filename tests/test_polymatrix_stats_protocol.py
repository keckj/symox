"""Official late stats methods leave exact PolyMatrix ground arithmetic native."""
import importlib
import itertools
import os
import subprocess
import sys

import pytest
import symox as sx
import sympy as sp
from sympy.polys.polymatrix import PolyMatrix as Reference

pm=importlib.import_module('symox.polys.polymatrix')
x,y,z=sx.symbols('x y z')


def stats():return importlib.import_module('symox.stats._symbolic')


def matrix(kind='QQ'):
    ring=sx.QQ if kind=='QQ' else sx.QQ.frac_field(y,z)
    data=[[1,2,3],[3,7,8]] if kind=='QQ' else [[y/(z+1),1,2],[2*y/(z+1),3,5]]
    return pm.PolyMatrix(data,x,ring=ring)


@pytest.mark.parametrize('order',list(itertools.permutations(('stats','vector','holonomic.recurrence'))))
def test_stats_vector_ore_import_orders_preserve_native_exact_ground(order):
    code='''
import importlib,sys
import symox as s
from symox.polys.polymatrix import PolyMatrix,_native_ground
assert 'symox.stats._symbolic' not in sys.modules
x,y=s.symbols('x y');a=PolyMatrix([[y,1,2],[2*y,3,5]],x,ring=s.QQ.frac_field(y))
before=_native_ground(a,True)
assert before is not None
assert 'symox.stats._symbolic' not in sys.modules
for name in sys.argv[1:]:importlib.import_module('symox.'+name)
after=_native_ground(a,True)
assert after is not None and after==before
'''
    subprocess.run([sys.executable,'-c',code,*order],check=True,env=os.environ.copy())


@pytest.mark.parametrize('kind',['QQ','RF'])
@pytest.mark.parametrize('basis',[False,True])
def test_stats_wrappers_keep_native_route_reference_and_certificates(monkeypatch,kind,basis):
    stats();a=matrix(kind);calls=[];native=sx._core._polymatrix_ground
    def spy(*args):
        result=native(*args);calls.append(result is not None);return result
    monkeypatch.setattr(sx._core,'_polymatrix_ground',spy)
    result=a.nullspace() if basis else a.rref()
    assert calls and all(calls)
    X,Y,Z=sp.symbols('x y z');field=sp.QQ if kind=='QQ' else sp.QQ.frac_field(Y,Z)
    flat=list(map(sx.to_sympy,a.flat()))
    raw=[flat[i*a.cols:(i+1)*a.cols] for i in range(a.rows)]
    ref=Reference(a.rows,a.cols,flat,ring=field[X])
    expected=ref.nullspace() if basis else ref.rref()
    if basis:
        assert len(result)==len(expected)
        for actual,want in zip(result,expected):
            assert all(sp.cancel(got-wanted)==0 for got,wanted in zip(map(sx.to_sympy,actual.flat()),list(want.to_Matrix())))
            assert all(sp.cancel(sum(raw[i][j]*sx.to_sympy(actual.flat()[j]) for j in range(a.cols)))==0 for i in range(a.rows))
    else:
        actual,pivots=result;want,rpivots=expected
        assert pivots==rpivots
        assert all(sp.cancel(got-wanted)==0 for got,wanted in zip(map(sx.to_sympy,actual.flat()),list(want.to_Matrix())))
    monkeypatch.setattr(pm,'_native_ground',lambda *args:None)
    assert result==(a.nullspace() if basis else a.rref())


@pytest.mark.parametrize('name',['doit','rewrite'])
def test_public_stats_method_override_retains_full_legacy_protocol(monkeypatch,name):
    stats();a=matrix('RF');original=getattr(sx.Expr,name);calls=[]
    def callback(self,*args,**kwargs):calls.append(name);return original(self,*args,**kwargs)
    monkeypatch.setattr(sx.Expr,name,callback)
    assert pm._native_ground(a) is None and calls==[]
    result=a.rref();observed=calls[:];calls.clear()
    monkeypatch.setattr(pm,'_native_ground',lambda *args:None)
    expected=a.rref();legacy=calls[:]
    assert legacy==observed==[]  # Neither method occurs on exact ground path.
    assert result==expected


@pytest.mark.parametrize('name',['doit','rewrite'])
@pytest.mark.parametrize('mode',['delegate','code','defaults','empty_cell'])
def test_mutated_stats_wrappers_decline_without_speculative_callbacks(monkeypatch,name,mode):
    module=stats();a=matrix();function=getattr(sx.Expr,name);calls=[]
    def callback(*args,**kwargs):calls.append(args);raise AssertionError('speculative stats callback')
    if mode=='code':
        factory=lambda original: lambda self,*args,**kwargs:original(self,*args,**kwargs)
        monkeypatch.setattr(function,'__code__',factory(callback).__code__)
    elif mode=='defaults':monkeypatch.setattr(function,'__defaults__',(callback,))
    elif mode=='delegate':monkeypatch.setattr(function.__closure__[0],'cell_contents',callback)
    else:
        cell=function.__closure__[0];old=cell.cell_contents;del cell.cell_contents
        try:assert pm._native_ground(a) is None and calls==[]
        finally:cell.cell_contents=old
        return
    assert pm._native_ground(a) is None and calls==[]


@pytest.mark.parametrize('mode',['installer','doit_template','rewrite_template','marker','unrelated','module'])
def test_other_stats_protocol_mutations_keep_conservative_fallback(monkeypatch,mode):
    module=stats();a=matrix();native_expr=importlib.import_module('symox._node_types').node_class('Expr')
    if mode=='installer':monkeypatch.setattr(module.install,'__code__',(lambda:None).__code__)
    elif mode.endswith('_template'):
        monkeypatch.setattr(module,'_STATS_'+mode.split('_')[0].upper()+'_CODE',(lambda:None).__code__)
    elif mode=='marker':monkeypatch.setattr(native_expr,'_stats_doit_installed',object())
    elif mode=='unrelated':monkeypatch.setattr(native_expr,'unrelated_callback',lambda self:None,raising=False)
    else:
        from types import ModuleType
        fake=ModuleType(module.__name__);fake.__dict__.update(module.__dict__)
        monkeypatch.setitem(sys.modules,module.__name__,fake)
    assert pm._native_ground(a) is None


def test_stats_registry_is_not_probed_by_commutative_ground(monkeypatch):
    module=stats();a=matrix('RF');calls=[]
    class Registry(dict):
        def __contains__(self,key):calls.append('contains');raise AssertionError('stats registry probe')
    monkeypatch.setattr(module,'HELD',Registry())
    assert pm._native_ground(a) is not None and calls==[]
    result=a.rref();assert calls==[]
    monkeypatch.setattr(pm,'_native_ground',lambda *args:None)
    assert a.rref()==result and calls==[]


def test_held_stats_coefficient_retains_complete_non_ground_failure(monkeypatch):
    module=stats();held=module.Expectation(y)
    a=pm.PolyMatrix([[held,1]],x,ring=sx.QQ)
    assert pm._native_ground(a) is None
    with pytest.raises(sx.CoercionFailed) as actual:a.rref()
    monkeypatch.setattr(pm,'_native_ground',lambda *args:None)
    with pytest.raises(type(actual.value)) as expected:a.rref()
    assert str(actual.value)==str(expected.value)
