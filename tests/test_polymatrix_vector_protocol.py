"""Optional official vector forwarding and observable PolyMatrix fallbacks."""
import importlib
import os
import subprocess
import sys

import pytest
import symox as sx

pm=importlib.import_module('symox.polys.polymatrix')
x=sx.Symbol('x')


def matrix():
    return pm.PolyMatrix([[1,2,3],[3,7,8]],x,ring=sx.QQ)


def adapters():
    return importlib.import_module('symox.vector._adapters')


@pytest.mark.parametrize('imports',[('vector',),('holonomic.recurrence',),
    ('vector','holonomic.recurrence'),('holonomic.recurrence','vector')])
def test_optional_import_does_not_change_selection_or_native_results(imports):
    code='''
import sys,importlib
import symox as s
from symox.polys.polymatrix import PolyMatrix,_native_ground
assert 'symox.vector._adapters' not in sys.modules
x=s.Symbol('x'); a=PolyMatrix([[1,2,3],[3,7,8]],x,ring=s.QQ)
before=_native_ground(a,True)
assert before is not None
assert 'symox.vector._adapters' not in sys.modules
for name in sys.argv[1:]:importlib.import_module('symox.'+name)
after=_native_ground(a,True)
assert after is not None and after==before
'''
    subprocess.run([sys.executable,'-c',code,*imports],check=True,env=os.environ.copy())


@pytest.mark.parametrize('basis',[False,True])
def test_stock_vector_forwarders_use_native_ground_and_old_exact_results(monkeypatch,basis):
    adapters();a=matrix();calls=[];original=sx._core._polymatrix_ground
    def spy(*args):
        result=original(*args);calls.append(result is not None);return result
    monkeypatch.setattr(sx._core,'_polymatrix_ground',spy)
    result=a.nullspace() if basis else a.rref()
    assert calls and all(calls)
    monkeypatch.setattr(pm,'_native_ground',lambda *args:None)
    assert result==(a.nullspace() if basis else a.rref())


@pytest.mark.parametrize('name',['__add__','__sub__','__mul__','__truediv__'])
def test_user_operator_wrapping_preserves_complete_callback_order(monkeypatch,name):
    adapters();a=matrix();calls=[];original=getattr(sx._core.Expr,name)
    def callback(self,other):
        calls.append(name);return original(self,other)
    monkeypatch.setattr(sx._core.Expr,name,callback)
    assert pm._native_ground(a) is None and calls==[]
    result=a.rref();observed=calls[:];calls.clear()
    monkeypatch.setattr(pm,'_native_ground',lambda *args:None)
    old=a.rref();legacy_calls=calls[:]
    assert legacy_calls==observed
    assert old==result
    if name!='__add__':assert observed


@pytest.mark.parametrize('mode',['delegate','predicate','operation_name','basis','empty_cell'])
def test_mutated_forwarder_never_invokes_callbacks_during_selection(monkeypatch,mode):
    module=adapters();a=matrix();calls=[]
    operation=next(p[0] for p in module._SCALAR_DISPATCH_PROOFS if p[7]=='__truediv__')
    cells=dict(zip(operation.__code__.co_freevars,operation.__closure__))
    def observed(*args):calls.append(args);raise AssertionError('speculative callback')
    if mode=='delegate':monkeypatch.setattr(cells['original'],'cell_contents',observed)
    elif mode=='predicate':monkeypatch.setattr(module,'isinstance',observed,raising=False)
    elif mode=='operation_name':monkeypatch.setattr(cells['name'],'cell_contents','__mul__')
    elif mode=='basis':monkeypatch.setattr(cells['_BasisDependent'],'cell_contents',sx.Expr)
    else:
        previous=cells['original'].cell_contents
        del cells['original'].cell_contents
        try:
            assert pm._native_ground(a) is None and calls==[]
        finally:cells['original'].cell_contents=previous
        return
    assert pm._native_ground(a) is None and calls==[]


@pytest.mark.parametrize('mode',['registry','proof','states','mro','installer','template','marker','extra_attribute'])
def test_invalid_certificates_and_unrelated_class_changes_keep_fallback(monkeypatch,mode):
    module=adapters()
    operation=next(p[0] for p in module._SCALAR_DISPATCH_PROOFS if p[7]=='__add__')
    proof=operation._symox_vector_proof
    if mode=='registry':monkeypatch.setattr(module,'_SCALAR_DISPATCH_PROOFS',[])
    elif mode=='proof':monkeypatch.setattr(operation,'_symox_vector_proof',())
    elif mode in ('states','mro'):
        changed=list(proof);changed[2 if mode=='states' else 5]=[object()]
        changed=tuple(changed)
        proofs=tuple(changed if p is proof else p for p in module._SCALAR_DISPATCH_PROOFS)
        monkeypatch.setattr(module,'_SCALAR_DISPATCH_PROOFS',proofs)
        monkeypatch.setattr(module.install,'_stock_proofs',proofs)
        monkeypatch.setattr(operation,'_symox_vector_proof',changed)
    elif mode=='installer':monkeypatch.setattr(module.install,'__code__',(lambda:None).__code__)
    elif mode=='template':monkeypatch.setattr(module,'_SCALAR_OPERATION_CODE',(lambda:None).__code__)
    elif mode=='marker':monkeypatch.setattr(sx._core.Expr,'_coordinate_vector_dispatch_installed',object())
    else:monkeypatch.setattr(sx._core.Expr,'unrelated_callback',lambda self:None,raising=False)
    assert pm._native_ground(matrix()) is None



def test_stock_rationaltools_reload_retains_authenticated_forwarding(monkeypatch):
    adapters();rt=importlib.import_module('symox.integrals.rationaltools')
    importlib.reload(rt)
    a=matrix()
    assert pm._native_ground(a) is not None
    original=rt._atan_vector_proof;calls=[]
    def callback(*args):calls.append(args);return original(*args)
    monkeypatch.setattr(rt,'_atan_vector_proof',callback)
    assert pm._native_ground(a) is None and calls==[]


def test_reloaded_guard_code_mutation_keeps_complete_fallback(monkeypatch):
    adapters();rt=importlib.import_module('symox.integrals.rationaltools')
    importlib.reload(rt)
    assert pm._native_ground(matrix()) is not None
    def callback(*args):raise AssertionError('must not call changed guard')
    monkeypatch.setattr(rt._atan_vector_proof,'__code__',callback.__code__)
    assert pm._native_ground(matrix()) is None


def test_guard_code_constants_do_not_invoke_host_serialization_or_equality(monkeypatch):
    adapters();rt=importlib.import_module('symox.integrals.rationaltools')
    importlib.reload(rt);a=matrix();calls=[]
    class Host:
        def __eq__(self,other):calls.append('eq');raise AssertionError('host equality')
        def __reduce__(self):calls.append('reduce');raise AssertionError('host serialization')
    function=rt._atan_vector_proof
    changed=function.__code__.replace(co_consts=function.__code__.co_consts+(Host(),))
    monkeypatch.setattr(function,'__code__',changed)
    assert pm._native_ground(a) is None and calls==[]
