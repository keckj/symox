"""Cold rejection inspects real ABC state without classification or warming."""
import importlib
import json
import os
from pathlib import Path
import subprocess
import sys
import types

import pytest
import symox as s

p=importlib.import_module('symox.polys.matrices._lll_native')


def run(script):
    r=subprocess.run([sys.executable,'-c',script],capture_output=True,text=True,timeout=40)
    assert r.returncode==0,r.stdout+r.stderr
    return r.stdout


@pytest.mark.parametrize('entry',['ddm','dense','sparse'])
def test_fresh_public_cold_path_does_not_scan_heavy_protocol_or_warm_cache(entry):
    script='''
import symox as s,importlib,_abc,json,os
from symox._basic import Basic
from symox.polys.matrices.ddm import DDM
from symox.polys.matrices.domainmatrix import DomainMatrix
p=importlib.import_module('symox.polys.matrices._lll_native')
# Preserve this original cold-rejection witness with only the new native
# certificate disabled. Separate tests require real first-call admission.
p.prefix,p.ready=p._LEGACY_PREFIX,p._LEGACY_READY
if os.environ.get('SYMOX_LLL_LEGACY_PREFIX'):
 exec('def original(x,delta,transform):\\n return (_CAPTURED and _TYPE(transform) is _builtins.bool and _TYPE(delta) is _ELEMENT and ready(x,x,None))',p.__dict__)
 p.prefix=p.original
if ENTRY=='ddm':a=DDM([[s.ZZ(2),s.ZZ(0)],[s.ZZ(3),s.ZZ(2)]],(2,2),s.ZZ)
else:
 a=DomainMatrix.from_list([[2,0],[3,2]],s.ZZ)
 if ENTRY=='sparse':a=a.to_sparse()
# Construction is real and there is no earlier LLL or QQ conversion.
assert not any(r() is s._core.DomainElement for r in _abc._get_dump(Basic)[2])
checks=[];native=[];old=p._protocol;call=s._core._ddm_lll_native
def protocol():checks.append(True);return old()
def spy(*args):r=call(*args);native.append(r is not None);return r
p._protocol=protocol;s._core._ddm_lll_native=spy
b,t=a.lll_transform()
assert checks==[],checks
assert native==[],native
assert any(r() is s._core.DomainElement for r in _abc._get_dump(Basic)[2])
a.lll_transform()
assert checks and native==[True],(checks,native)
rows=b.to_list() if hasattr(b,'to_list') else b
assert [[int(v) for v in row] for row in rows]==[[2,0],[-1,2]]
'''.replace('ENTRY',repr(entry))
    run(script)


def fingerprint(state):
    return tuple(tuple(sorted(id(v) for v in group)) for group in state[:3])+(state[3],)


def test_inspection_does_not_populate_or_invalidate_abc_caches():
    script='''
import symox as s,_abc
from symox._basic import Basic
f=s._core._lll_cache_is_cold
fp=lambda t:tuple(tuple(sorted(map(id,g))) for g in t[:3])+(t[3],)
a=fp(_abc._get_dump(Basic))
assert f(Basic,s._core.DomainElement,_abc._get_dump,_abc) is True
assert fp(_abc._get_dump(Basic))==a
# Only this actual classifier invocation warms the cache.
assert not isinstance(s.ZZ(2),Basic)
a=fp(_abc._get_dump(Basic))
assert f(Basic,s._core.DomainElement,_abc._get_dump,_abc) is False
assert fp(_abc._get_dump(Basic))==a
'''
    run(script)


@pytest.mark.parametrize('kind',['python_dump','foreign_c_dump','foreign_owner','owner_subclass','module_bad_key','class_bad_key','bad_data','nonclass'])
def test_untrusted_inspection_inputs_decline_without_callbacks(kind):
    import _abc
    from symox._basic import Basic
    from symox._basic import _BasicMeta
    f=s._core._lll_cache_is_cold;events=[]
    dump,owner,basic=_abc._get_dump,_abc,Basic
    if kind=='python_dump':
        def dump(value):events.append('dump');raise AssertionError('custom dump called')
    elif kind=='foreign_c_dump':dump=len
    elif kind=='foreign_owner':owner=types.ModuleType('_abc')
    elif kind=='owner_subclass':
        class Module(types.ModuleType):
            def __getattribute__(self,name):events.append(name);raise AssertionError(name)
        owner=Module('_abc')
    elif kind=='module_bad_key':
        # A fake owner is rejected even before reading arbitrary dictionary keys.
        class Key:
            def __hash__(self):return 1
            def __eq__(self,other):events.append('eq');raise AssertionError('key compared')
        owner=types.ModuleType('_abc');owner.__dict__[Key()]=None
    else:
        class Meta(type):
            def __getattribute__(cls,name):events.append(name);raise AssertionError(name)
        class Cold(metaclass=Meta):pass
        if kind=='bad_data':type.__setattr__(Cold,'_abc_impl',object())
        if kind=='class_bad_key':
            # setattr normalizes identifiers, but a string subclass must not
            # be used for a user equality comparison by the inspector.
            class Key(str):
                def __eq__(self,other):events.append('eq');raise AssertionError('key compared')
                __hash__=str.__hash__
            type.__setattr__(Cold,Key('extra'),None)
        basic=Cold if kind!='nonclass' else object()
    assert f(basic,s._core.DomainElement,dump,owner) is None
    assert events==[]


def test_basic_metaclass_attribute_hooks_are_not_invoked_by_the_inspector():
    import _abc
    from symox._basic import Basic
    class Meta(type):
        def __getattribute__(cls,name):raise AssertionError('inspector used a Python metaclass hook: '+name)
    class Class(metaclass=Meta):
        _abc_impl=Basic.__dict__['_abc_impl']
    before=fingerprint(_abc._get_dump(types.SimpleNamespace(_abc_impl=Basic.__dict__['_abc_impl'])))
    result=s._core._lll_cache_is_cold(Class,s._core.DomainElement,_abc._get_dump,_abc)
    assert result in (True,False)
    after=fingerprint(_abc._get_dump(types.SimpleNamespace(_abc_impl=Basic.__dict__['_abc_impl'])))
    assert before==after


def test_replaced_inspector_is_not_called_by_public_selection(monkeypatch):
    from symox.polys.matrices.ddm import DDM
    a=DDM([[s.ZZ(2),s.ZZ(0)],[s.ZZ(3),s.ZZ(2)]],(2,2),s.ZZ);events=[]
    def replaced(*args):events.append('replaced');raise AssertionError('replaced inspector')
    monkeypatch.setattr(s._core,'_lll_cache_is_cold',replaced)
    b,t=a.lll_transform()
    assert [[int(c) for c in row] for row in b]==[[2,0],[-1,2]]
    assert events==[]


def test_no_reference_imports_in_fresh_cold_rejection():
    run('''
import sys
class Block:
 def find_spec(self,name,*args):
  if name=='sympy' or name.startswith('sympy.'):raise AssertionError(name)
sys.meta_path.insert(0,Block())
import symox as s,_abc
from symox._basic import Basic
from symox.polys.matrices.domainmatrix import DomainMatrix
a=DomainMatrix.from_list([[2,0],[3,2]],s.ZZ);b,t=a.lll_transform()
assert [[int(c) for c in row] for row in b.to_list()]==[[2,0],[-1,2]]
assert not any(n=='sympy' or n.startswith('sympy.') for n in sys.modules)
''')


def test_removed_lll_module_cannot_reuse_the_old_namespace_proof():
    run('''
import symox as s,sys,importlib
from symox.polys.matrices.ddm import DDM
p=importlib.import_module('symox.polys.matrices._lll_native')
# An actual coefficient conversion makes this a warm-routing regression.
s.QQ.convert_from(s.ZZ(1),s.ZZ)
sys.modules.pop('symox.polys.matrices.lll')
lll=importlib.import_module('symox.polys.matrices.lll')
events=[];original=lll.mfloor
def mfloor(value):events.append(str(value));return original(value)
lll.mfloor=mfloor
a=DDM([[s.ZZ(2),s.ZZ(0)],[s.ZZ(3),s.ZZ(2)]],(2,2),s.ZZ)
actual=lll.ddm_lll_transform(a);trace=events[:];events.clear()
p.prefix=lambda *args:False
expected=lll.ddm_lll_transform(a)
assert actual==expected
assert trace==events and trace,(trace,events)
''')


def test_removed_ddm_module_keeps_its_real_constructor_callbacks():
    run('''
import symox as s,sys,importlib
from symox.polys.matrices.ddm import DDM
p=importlib.import_module('symox.polys.matrices._lll_native')
s.QQ.convert_from(s.ZZ(1),s.ZZ)
sys.modules.pop('symox.polys.matrices.ddm')
ddm=importlib.import_module('symox.polys.matrices.ddm')
events=[];original=ddm.DDM.copy
def copy(self):events.append('copy');return original(self)
ddm.DDM.copy=copy
a=ddm.DDM([[s.ZZ(2),s.ZZ(0)],[s.ZZ(3),s.ZZ(2)]],(2,2),s.ZZ)
actual=a.lll_transform();trace=events[:];events.clear()
p.prefix=lambda *args:False
expected=a.lll_transform()
assert actual==expected and trace==events and trace,(trace,events)
''')


@pytest.mark.parametrize('name',['max','sum','abs','all'])
def test_builtin_replaced_before_import_retains_the_actual_old_callbacks(name):
    run('''
import builtins
original=builtins.__dict__[NAME];events=[]
def wrapper(*args,**kwargs):events.append(NAME);return original(*args,**kwargs)
builtins.__dict__[NAME]=wrapper
try:
 import symox as s
 from symox.polys.matrices.ddm import DDM
 from symox.polys.matrices import _lll_native as p
 # This matrix actually swaps rows, so max is genuinely consumed.
 DDM([[s.ZZ(1)]],(1,1),s.ZZ).lll_transform()
 a=DDM([[s.ZZ(3),s.ZZ(2)],[s.ZZ(2),s.ZZ(0)]],(2,2),s.ZZ)
 events.clear();actual=a.lll_transform();trace=events[:];events.clear()
 p.prefix=lambda *args:False
 expected=a.lll_transform()
 assert actual==expected
 assert trace==events and events,(NAME,trace,events)
finally:builtins.__dict__[NAME]=original
'''.replace('NAME',repr(name)))


@pytest.mark.parametrize('name',['max','sum','abs','all','int','tuple','type','enumerate'])
def test_bootstrap_authentication_refuses_python_wrappers_and_heap_factories(name):
    import builtins
    events=[]
    if name in ('max','sum','abs','all'):
        def value(*args,**kwargs):events.append('call');return builtins.__dict__[name](*args,**kwargs)
    else:
        class Meta(type):
            def __getattribute__(cls,attribute):events.append(attribute);raise AssertionError(attribute)
        value=Meta(name,(),{'__module__':'builtins'})
    assert s._core._lll_bootstrap_stock((({},name,value),)) is False
    assert events==[]


def test_bootstrap_authentication_accepts_actual_immutable_primitives():
    import builtins
    names=('max','sum','abs','all','isinstance','getattr','hasattr','len','any',
           'int','tuple','type','enumerate','super','dict','bool','range')
    records=tuple(({},n,builtins.__dict__[n]) for n in names)
    assert s._core._lll_bootstrap_stock(records) is True


@pytest.mark.parametrize('alias',['_TYPE','_ALL','_ANY','_LEN','_ZIP','_DICT','_TUPLE','_GETATTR'])
def test_helper_alias_changed_during_import_cannot_authenticate_itself(alias):
    run('''
import sys,builtins
original_import=builtins.__import__;seen=[];events=[]
name='symox.polys.matrices._lll_native'
def importing(*args,**kwargs):
 result=original_import(*args,**kwargs)
 helper=sys.modules.get(name)
 if helper is not None and ALIAS in helper.__dict__ and not seen:
  seen.append(True);original=helper.__dict__[ALIAS]
  def observer(*values,**options):events.append(ALIAS);return original(*values,**options)
  helper.__dict__[ALIAS]=observer
  builtins.__import__=original_import
 return result
builtins.__import__=importing
try:
 import symox as s
 from symox.polys.matrices.ddm import DDM
 from symox.polys.matrices import _lll_native as p
 assert seen
 DDM([[s.ZZ(1)]],(1,1),s.ZZ).lll_transform()
 a=DDM([[s.ZZ(3),s.ZZ(2)],[s.ZZ(2),s.ZZ(0)]],(2,2),s.ZZ)
 events.clear();actual=a.lll_transform();trace=events[:];events.clear()
 p.prefix=lambda *args:False
 expected=a.lll_transform()
 assert actual==expected and trace==events==[],(ALIAS,trace,events)
finally:builtins.__import__=original_import
'''.replace('ALIAS',repr(alias)))


@pytest.mark.parametrize('tool',['_get_dump','get_cache_token'])
@pytest.mark.parametrize('replacement',['python','foreign_c'])
def test_abc_inspector_replaced_before_import_is_never_called_by_selection(tool,replacement):
    run('''
import _abc,math
original=getattr(_abc,TOOL);active=False;events=[]
def observer(*args):
 if active:events.append(TOOL);raise AssertionError('replacement ABC inspector called')
 return original(*args)
setattr(_abc,TOOL,observer if KIND=='python' else math.floor)
try:
 import symox as s
 from symox.polys.matrices.ddm import DDM
 from symox.polys.matrices import _lll_native as p
 # This real conversion targets warm admission without a prior LLL call.
 s.QQ.convert_from(s.ZZ(1),s.ZZ)
 a=DDM([[s.ZZ(3),s.ZZ(2)],[s.ZZ(2),s.ZZ(0)]],(2,2),s.ZZ)
 active=True;actual=a.lll_transform();trace=events[:];events.clear()
 p.prefix=lambda *args:False;expected=a.lll_transform()
 assert actual==expected and trace==events==[],(TOOL,KIND,trace,events)
finally:setattr(_abc,TOOL,original)
'''.replace('TOOL',repr(tool)).replace('KIND',repr(replacement)))
