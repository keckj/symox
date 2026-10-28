"""Probe metadata names are not coefficients and must not run encoding/equality."""
import subprocess
import sys
import types

import pytest
import symox as s


@pytest.mark.parametrize('key',['\ud800','unused_\udfff','mixed_\ud800_\udfff'])
def test_unread_surrogate_class_key_does_not_change_a_readonly_probe(key):
    import _abc
    from symox._basic import Basic
    # A private holder shares only the real opaque ABC data. No classifier runs.
    shadow=type('Shadow',(),{'_abc_impl':Basic.__dict__['_abc_impl'],key:None})
    fp=lambda state:tuple(tuple(sorted(map(id,g))) for g in state[:3])+(state[3],)
    before=fp(_abc._get_dump(types.SimpleNamespace(_abc_impl=Basic.__dict__['_abc_impl'])))
    result=s._core._lll_cache_is_cold(shadow,s._core.DomainElement,_abc._get_dump,_abc)
    after=fp(_abc._get_dump(types.SimpleNamespace(_abc_impl=Basic.__dict__['_abc_impl'])))
    assert type(result) is bool and before==after


@pytest.mark.parametrize('key',['\ud800','unused_\udfff','mixed_\ud800_\udfff'])
def test_true_public_cold_call_ignores_unread_surrogate_basic_attribute(key):
    script='''
import sys
class NoReference:
 def find_spec(self,name,*args):
  if name=='sympy' or name.startswith('sympy.'):raise AssertionError(name)
sys.meta_path.insert(0,NoReference())
import symox as s
from symox._basic import Basic
from symox.polys.matrices.ddm import DDM
from symox.polys.matrices import _lll_native as p
type.__setattr__(Basic,KEY,None)
a=DDM([[s.ZZ(3),s.ZZ(2)],[s.ZZ(2),s.ZZ(0)]],(2,2),s.ZZ)
actual=a.lll_transform()
p.prefix=lambda *args:False
expected=a.lll_transform()
assert actual==expected
assert [[int(c) for c in row] for row in actual[0]]==[[2,0],[-1,2]]
assert not any(n=='sympy' or n.startswith('sympy.') for n in sys.modules)
'''.replace('KEY',repr(key))
    result=subprocess.run([sys.executable,'-c',script],capture_output=True,text=True,timeout=40)
    assert result.returncode==0,result.stdout+result.stderr


@pytest.mark.parametrize('foreign_first',[False,True])
def test_all_keys_are_validated_before_the_abc_data_lookup(foreign_first):
    import _abc
    from symox._basic import Basic
    events=[];armed=False
    class Key(str):
        def __hash__(self):return hash('_abc_impl')
        def __eq__(self,other):
            if armed:events.append('equality');raise AssertionError('untrusted key equality')
            return False
    foreign=Key('foreign')
    pairs=[(foreign,None),('_abc_impl',Basic.__dict__['_abc_impl'])]
    if not foreign_first:pairs.reverse()
    shadow=type('Shadow',(),dict(pairs))
    assert any(type(k) is Key for k in type.__dict__['__dict__'].__get__(shadow))
    armed=True
    result=s._core._lll_cache_is_cold(shadow,s._core.DomainElement,_abc._get_dump,_abc)
    assert result is None and events==[]


def test_missing_abc_data_is_an_unsupported_probe_without_attribute_hooks():
    import _abc
    events=[]
    class Meta(type):
        def __getattribute__(cls,name):events.append(name);raise AssertionError(name)
    class Shadow(metaclass=Meta):pass
    assert s._core._lll_cache_is_cold(Shadow,s._core.DomainElement,_abc._get_dump,_abc) is None
    assert events==[]


@pytest.mark.parametrize('key',['\ud800','unused_\udfff','mixed_\ud800_\udfff'])
def test_unread_surrogate_module_key_does_not_change_a_readonly_probe(monkeypatch,key):
    import _abc
    from symox._basic import Basic
    fp=lambda state:tuple(tuple(sorted(map(id,g))) for g in state[:3])+(state[3],)
    before=fp(_abc._get_dump(Basic))
    monkeypatch.setitem(_abc.__dict__,key,None)
    result=s._core._lll_cache_is_cold(Basic,s._core.DomainElement,_abc._get_dump,_abc)
    assert type(result) is bool and fp(_abc._get_dump(Basic))==before


@pytest.mark.parametrize('key',['\ud800','unused_\udfff','mixed_\ud800_\udfff'])
def test_true_public_cold_call_ignores_unread_surrogate_abc_module_attribute(key):
    script='''
import symox as s,_abc
from symox.polys.matrices.ddm import DDM
from symox.polys.matrices import _lll_native as p
_abc.__dict__[KEY]=None
a=DDM([[s.ZZ(3),s.ZZ(2)],[s.ZZ(2),s.ZZ(0)]],(2,2),s.ZZ)
actual=a.lll_transform();p.prefix=lambda *args:False;expected=a.lll_transform()
assert actual==expected
assert [[int(c) for c in row] for row in actual[0]]==[[2,0],[-1,2]]
'''.replace('KEY',repr(key))
    result=subprocess.run([sys.executable,'-c',script],capture_output=True,text=True,timeout=40)
    assert result.returncode==0,result.stdout+result.stderr


@pytest.mark.parametrize('kind',['surrogate','object','str_subclass'])
def test_uninterpretable_module_name_declines_without_encoding_or_hooks(monkeypatch,kind):
    import _abc
    from symox._basic import Basic
    events=[]
    if kind=='surrogate':name='\ud800'
    elif kind=='object':
        class Name:
            def __str__(self):events.append('str');raise AssertionError('name coerced')
        name=Name()
    else:
        class Name(str):
            def __eq__(self,other):events.append('eq');raise AssertionError('subclass equality')
            __hash__=str.__hash__
        name=Name('_abc')
    monkeypatch.setitem(_abc.__dict__,'__name__',name)
    assert s._core._lll_cache_is_cold(Basic,s._core.DomainElement,_abc._get_dump,_abc) is None
    assert events==[]


def test_module_keys_are_all_validated_before_name_lookup():
    import _abc
    from symox._basic import Basic
    events=[];armed=False
    class Key:
        def __hash__(self):return hash('__name__')
        def __eq__(self,other):
            if armed:events.append('eq');raise AssertionError('foreign module key compared')
            return False
    key=Key();name=_abc.__dict__.pop('__name__')
    try:
        _abc.__dict__[key]=None;_abc.__dict__['__name__']=name;armed=True
        assert s._core._lll_cache_is_cold(Basic,s._core.DomainElement,_abc._get_dump,_abc) is None
        assert events==[]
    finally:
        armed=False;del _abc.__dict__[key];_abc.__dict__['__name__']=name


def test_readonly_probe_emits_no_audit_event_and_does_not_warm_cache():
    script='''
import symox as s,sys,_abc
from symox._basic import Basic
fp=lambda state:tuple(tuple(sorted(map(id,g))) for g in state[:3])+(state[3],)
before=fp(_abc._get_dump(Basic));events=[]
sys.addaudithook(lambda name,args:events.append(name))
result=s._core._lll_cache_is_cold(Basic,s._core.DomainElement,_abc._get_dump,_abc)
probe_events=events[:];events.clear()
after=fp(_abc._get_dump(Basic))
assert type(result) is bool and before==after
assert probe_events==[],probe_events
'''
    result=subprocess.run([sys.executable,'-c',script],capture_output=True,text=True,timeout=40)
    assert result.returncode==0,result.stdout+result.stderr
