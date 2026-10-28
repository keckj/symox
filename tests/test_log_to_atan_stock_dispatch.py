"""Authenticate the official Ore forwarding layer without trusting wraps."""
import importlib
import functools
import subprocess
import sys

import pytest
import symox as sx

rt=importlib.import_module('symox.integrals.rationaltools')
rec=importlib.import_module('symox.holonomic.recurrence')
x=sx.Symbol('x')


def pair():return sx.Poly(x**2+1,x,domain=sx.QQ),sx.Poly(x-1,x,domain=sx.QQ)


def proof(name):
    # Vector arithmetic can be installed after Ore during test collection.
    # Mutate the same Ore layer even when it is no longer the outer wrapper.
    current=getattr(sx._core.Expr,name);seen=set()
    while id(current) not in seen:
        seen.add(id(current))
        record=next((p for p in rec._SCALAR_DISPATCH_PROOFS if p[0] is current),None)
        if record is not None:return record
        vector=sys.modules.get('symox.vector._adapters')
        records=() if vector is None else vector._SCALAR_DISPATCH_PROOFS
        record=next((p for p in records if p[0] is current),None)
        if record is None:
            geometry=sys.modules.get('symox.diffgeom._adapters')
            records=() if geometry is None else geometry._SCALAR_DISPATCH_PROOFS
            record=next((p for p in records if p[0] is current),None)
        assert record is not None, 'expected an installed official arithmetic wrapper'
        current=record[1]
    raise AssertionError('cyclic arithmetic wrappers in test setup')


def forbid_native(monkeypatch):
    def forbidden(*args):raise AssertionError('nonstock protocol reached native route')
    monkeypatch.setattr(sx._core,'_log_to_atan',forbidden)


def assert_old_fallback_and_hook(name,message,f,g,expected):
    if name=='__add__':
        with pytest.raises(KeyboardInterrupt,match=message):rt.log_to_atan(f,g)
    else:
        # The old recurrence uses Python-int * atan via __rmul__, so this
        # particular pair does not exercise Expr.__mul__. Preserve its full
        # result, and separately exercise a scalar multiplication that does.
        assert rt.log_to_atan(f,g)==expected
        with pytest.raises(KeyboardInterrupt,match=message):sx.Integer(2)*sx.atan(x)


def test_official_forwarding_selects_the_native_route_and_preserves_exact_result(monkeypatch):
    f,g=pair();calls=[];original=sx._core._log_to_atan
    expected=original(f,g)
    def selected(*args):calls.append(True);return original(*args)
    monkeypatch.setattr(sx._core,'_log_to_atan',selected)
    assert rt._stock_atan_protocol(f,g)
    assert rt.log_to_atan(f,g)==expected
    assert calls==[True]
    importlib.reload(rt)
    assert rt._stock_atan_protocol(f,g)
    assert rt.log_to_atan(f,g)==expected
    assert calls==[True,True]


@pytest.mark.parametrize('name',['__add__','__mul__'])
def test_wraps_metadata_cannot_authenticate_a_user_wrapper(name,monkeypatch):
    f,g=pair()
    expected=sx._core._log_to_atan(f,g)
    original=getattr(sx._core.Expr,name)
    @functools.wraps(original)
    def forged(*args):raise KeyboardInterrupt('user arithmetic wrapper')
    monkeypatch.setattr(sx._core.Expr,name,forged)
    assert not rt._stock_atan_protocol(f,g)
    forbid_native(monkeypatch)
    assert_old_fallback_and_hook(name,'user arithmetic wrapper',f,g,expected)


def interrupting_code(freevars,message):
    # A replacement code object with the exact closure arity of the official
    # function, allowing an in-place mutation without replacing its identity.
    names=','.join('v%d'%i for i in range(freevars));namespace={}
    use='\n        if ('+names+',):pass' if freevars else ''
    exec('def factory('+names+'):\n    def changed(self,other,*args):'+use+'\n        raise KeyboardInterrupt('+repr(message)+')\n    return changed',namespace)
    return namespace['factory'](*([None]*freevars)).__code__


@pytest.mark.parametrize('name',['__add__','__mul__'])
def test_official_wrapper_identity_does_not_hide_mutated_code(name,monkeypatch):
    f,g=pair()
    expected=sx._core._log_to_atan(f,g)
    operation=proof(name)[0]
    monkeypatch.setattr(operation,'__code__',interrupting_code(len(operation.__closure__ or ()),'mutated wrapper code'))
    assert not rt._stock_atan_protocol(f,g)
    forbid_native(monkeypatch)
    assert_old_fallback_and_hook(name,'mutated wrapper code',f,g,expected)


@pytest.mark.parametrize('name',['__add__','__mul__'])
def test_official_wrapper_identity_does_not_hide_a_replaced_delegate_cell(name,monkeypatch):
    f,g=pair()
    expected=sx._core._log_to_atan(f,g)
    operation=proof(name)[0]
    index=operation.__code__.co_freevars.index('original')
    def callback(*args):raise KeyboardInterrupt('mutated delegate cell')
    monkeypatch.setattr(operation.__closure__[index],'cell_contents',callback)
    assert not rt._stock_atan_protocol(f,g)
    forbid_native(monkeypatch)
    assert_old_fallback_and_hook(name,'mutated delegate cell',f,g,expected)


@pytest.mark.parametrize('name',['__add__','__mul__'])
def test_official_delegate_identity_does_not_hide_mutated_code(name,monkeypatch):
    f,g=pair()
    expected=sx._core._log_to_atan(f,g)
    delegate=proof(name)[1]
    monkeypatch.setattr(delegate,'__code__',interrupting_code(len(delegate.__closure__ or ()),'mutated delegate code'))
    assert not rt._stock_atan_protocol(f,g)
    forbid_native(monkeypatch)
    assert_old_fallback_and_hook(name,'mutated delegate code',f,g,expected)


@pytest.mark.parametrize('name',['isinstance','getattr'])
def test_wrapper_global_callbacks_remain_observable(name,monkeypatch):
    f,g=pair()
    def callback(*args):raise KeyboardInterrupt('dispatch global callback')
    monkeypatch.setattr(rec,name,callback,raising=False)
    assert not rt._stock_atan_protocol(f,g)
    forbid_native(monkeypatch)
    with pytest.raises(KeyboardInterrupt,match='dispatch global callback'):rt.log_to_atan(f,g)


def test_holonomic_marker_descriptors_are_not_bypassed(monkeypatch):
    f,g=pair()
    def callback(*args):raise KeyboardInterrupt('holonomic marker descriptor')
    monkeypatch.setattr(sx._core.Expr,'_symox_holonomic_function',property(callback),raising=False)
    assert not rt._stock_atan_protocol(f,g)
    forbid_native(monkeypatch)
    with pytest.raises(KeyboardInterrupt,match='holonomic marker descriptor'):rt.log_to_atan(f,g)


def test_native_scalar_attribute_hook_is_not_bypassed(monkeypatch):
    f,g=pair();expected=sx._core._log_to_atan(f,g)
    cls=type(sx.atan(x))
    def attribute(self,name):
        if name=='_symox_holonomic_function':raise KeyboardInterrupt('scalar attribute hook')
        return object.__getattribute__(self,name)
    monkeypatch.setattr(cls,'__getattribute__',attribute,raising=False)
    assert not rt._stock_atan_protocol(f,g)
    forbid_native(monkeypatch)
    assert rt.log_to_atan(f,g)==expected
    with pytest.raises(KeyboardInterrupt,match='scalar attribute hook'):sx.Integer(2)*sx.atan(x)


def test_host_generators_keep_the_general_route(monkeypatch):
    f,g=pair()
    class X(sx.Symbol):pass
    z=X('z');f=sx.Poly(z+1,z,domain=sx.QQ);g=sx.Poly(1,z,domain=sx.QQ)
    assert not rt._stock_atan_protocol(f,g)
    forbid_native(monkeypatch)
    assert rt.log_to_atan(f,g)==2*sx.atan(z+1)


def test_ore_and_holonomic_operands_keep_their_reflected_dispatch(monkeypatch):
    f,g=pair()
    from symox.holonomic import DifferentialOperators,HolonomicFunction
    _,D=DifferentialOperators(sx.QQ.old_poly_ring(x),'D')
    h=HolonomicFunction(D-1,x,0,[1]);sentinel=object();calls=[]
    def reflected(self,other):calls.append((self,other));return sentinel
    monkeypatch.setattr(type(D),'__rmul__',reflected)
    monkeypatch.setattr(type(h),'__rmul__',reflected)
    assert sx.Integer(2)*D is sentinel
    assert sx.Integer(2)*h is sentinel
    assert calls==[(D,sx.Integer(2)),(h,sx.Integer(2))]
    assert rt._stock_atan_protocol(f,g)


@pytest.mark.parametrize('order',['atan_first','ore_first'])
def test_fresh_process_import_orders_select_native_without_sympy(order):
    imports=("import symox.integrals.rationaltools as rt\nimport symox.holonomic" if order=='atan_first' else
             "import symox.holonomic\nimport symox.integrals.rationaltools as rt")
    script="import sys\nsys.modules['sympy']=None\nimport symox as s\n"+imports+'''
x=s.Symbol('x');f=s.Poly(x+1,x,domain=s.QQ);g=s.Poly(1,x,domain=s.QQ)
assert rt._stock_atan_protocol(f,g)
old=s._core._log_to_atan;calls=[]
def selected(*args):calls.append(True);return old(*args)
s._core._log_to_atan=selected
assert rt.log_to_atan(f,g)==2*s.atan(x+1)
assert calls==[True]
'''
    subprocess.run([sys.executable,'-c',script],check=True)


@pytest.mark.parametrize('forgery',[None,(),((),),(('invalid',),),'forged'])
def test_replaced_or_malformed_proof_registry_declines_without_selection_error(forgery,monkeypatch):
    f,g=pair();expected=sx._core._log_to_atan(f,g)
    monkeypatch.setattr(rec,'_SCALAR_DISPATCH_PROOFS',forgery)
    assert not rt._stock_atan_protocol(f,g)
    forbid_native(monkeypatch)
    assert rt.log_to_atan(f,g)==expected


def test_registry_tuple_copy_cannot_supply_a_new_certificate(monkeypatch):
    f,g=pair();expected=sx._core._log_to_atan(f,g)
    monkeypatch.setattr(rec,'_SCALAR_DISPATCH_PROOFS',tuple(list(rec._SCALAR_DISPATCH_PROOFS)))
    assert not rt._stock_atan_protocol(f,g)
    forbid_native(monkeypatch)
    assert rt.log_to_atan(f,g)==expected


def test_wrapper_certificate_metadata_replacement_is_not_trusted(monkeypatch):
    f,g=pair();operation=proof('__add__')[0]
    monkeypatch.setattr(operation,'_symox_ore_proof',None)
    assert not rt._stock_atan_protocol(f,g)
    forbid_native(monkeypatch)
    assert rt.log_to_atan(f,g)==2*sx.atan(x+1)-2*sx.atan(x+x**2/2+x**3/2)


@pytest.mark.parametrize('scenario',['recurrence_reload','rationaltools_reimport_with_callback'])
def test_module_reloads_or_removed_bootstrap_records_keep_a_conservative_fallback(scenario):
    script='''
import importlib,sys
import symox as s
import symox.integrals.rationaltools as rt
import symox.holonomic.recurrence as rec
x=s.Symbol('x');f=s.Poly(x**2+1,x,domain=s.QQ);g=s.Poly(x-1,x,domain=s.QQ)
assert rt._stock_atan_protocol(f,g)
'''
    if scenario=='recurrence_reload':
        script+='''
expected=s._core._log_to_atan(f,g)
importlib.reload(rec)
assert not rt._stock_atan_protocol(f,g)
def forbidden(*args):raise AssertionError('reload reached native route')
s._core._log_to_atan=forbidden
assert rt.log_to_atan(f,g)==expected
'''
    else:
        script+='''
def callback(*args):raise KeyboardInterrupt('late captured callback')
s._core.Expr.__add__=callback
sys.modules.pop('symox.integrals.rationaltools')
rt=importlib.import_module('symox.integrals.rationaltools')
assert not rt._stock_atan_protocol(f,g)
def forbidden(*args):raise AssertionError('reimport reached native route')
s._core._log_to_atan=forbidden
try:rt.log_to_atan(f,g)
except KeyboardInterrupt as error:assert str(error)=='late captured callback'
else:raise AssertionError('callback bypassed')
'''
    subprocess.run([sys.executable,'-c',script],check=True)


def test_a_fully_forged_registry_record_still_cannot_authenticate_user_code(monkeypatch):
    from types import FunctionType
    f,g=pair();old=proof('__add__')
    def callback(self,other):raise KeyboardInterrupt('forged complete record')
    forged=FunctionType(callback.__code__,rec.__dict__,old[0].__name__)
    functools.update_wrapper(forged,old[0])
    fake=(forged,old[1],rec._scalar_dispatch_states(forged),rec.__dict__,*old[4:])
    forged._symox_ore_proof=fake
    registry=tuple(fake if p is old else p for p in rec._SCALAR_DISPATCH_PROOFS)
    monkeypatch.setattr(rec,'_SCALAR_DISPATCH_PROOFS',registry)
    monkeypatch.setattr(rec._install_scalar_dispatch,'_stock_proofs',registry)
    monkeypatch.setattr(sx._core.Expr,'__add__',forged)
    assert not rt._stock_atan_protocol(f,g)
    forbid_native(monkeypatch)
    with pytest.raises(KeyboardInterrupt,match='forged complete record'):rt.log_to_atan(f,g)


@pytest.mark.parametrize('forgery',[None,(),((),),(('invalid',),),'forged'])
def test_malformed_registry_and_installer_anchor_cannot_raise_a_new_selection_error(forgery,monkeypatch):
    f,g=pair();expected=sx._core._log_to_atan(f,g)
    monkeypatch.setattr(rec,'_SCALAR_DISPATCH_PROOFS',forgery)
    monkeypatch.setattr(rec._install_scalar_dispatch,'_stock_proofs',forgery)
    assert not rt._stock_atan_protocol(f,g)
    forbid_native(monkeypatch)
    assert rt.log_to_atan(f,g)==expected
