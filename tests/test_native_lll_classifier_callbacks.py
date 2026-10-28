"""First-conversion classifiers may expose the entire scratch workspace."""
import importlib
import inspect
import subprocess
import sys
import pytest
import symox as s
from symox.polys.matrices.ddm import DDM

p=importlib.import_module('symox.polys.matrices._lll_native')
basic=importlib.import_module('symox._basic')


def raw(value):return [[str(c) for c in row] for row in value]


def test_temporary_classifier_preserves_observed_scratch_state(monkeypatch):
    old=basic._BasicMeta.__instancecheck__
    def run(enable):
        a=DDM([[s.ZZ(2),s.ZZ(0)],[s.ZZ(3),s.ZZ(2)]],(2,2),s.ZZ);scratch=[];events=[]
        def classify(cls,value):
            basic._BasicMeta.__instancecheck__=old
            frame=inspect.currentframe()
            while frame and frame.f_code.co_name!='_ddm_lll':frame=frame.f_back
            scratch.extend([frame.f_locals['mu'],frame.f_locals['y_star']]);events.append('classify')
            return old(cls,value)
        with monkeypatch.context() as patch:
            patch.setattr(basic._BasicMeta,'__instancecheck__',classify)
            if not enable:patch.setattr(p,'prefix',lambda *args:False)
            b,t=a.lll_transform()
        return raw(b),raw(t),[raw(v) for v in scratch],events
    actual=run(True);expected=run(False)
    assert actual==expected
    assert actual[2]==[[['0','0'],['-1/2','0']],[['2','0'],['0','2']]]
    assert actual[3]==['classify']


def test_untrusted_abc_alias_is_not_probed_by_selection(monkeypatch):
    import abc
    a=DDM([[s.ZZ(2),s.ZZ(0)],[s.ZZ(3),s.ZZ(2)]],(2,2),s.ZZ);events=[]
    class Foreign:
        def __getattribute__(self,name):events.append(name);raise AssertionError('selector probed an untrusted alias')
    monkeypatch.setattr(abc,'ABCMeta',Foreign())
    assert not p.prefix(a,s.QQ(3,4),True) and events==[]
    b,t=a.lll_transform()
    assert raw(b)==[['2','0'],['-1','2']] and events==[]


def test_an_abc_registry_callback_keeps_full_old_math():
    script='''
import symox as s,inspect,abc,importlib
from symox.polys.matrices.ddm import DDM
from symox._basic import Basic
p=importlib.import_module('symox.polys.matrices._lll_native')
scratch=[]
class Meta(type):
 def __subclasscheck__(cls,sub):
  if sub is s._core.DomainElement:
   frame=inspect.currentframe()
   while frame and frame.f_code.co_name!='_ddm_lll':frame=frame.f_back
   if frame:
    scratch.extend([frame.f_locals['mu'],frame.f_locals['y_star']])
  return False
class Foreign(metaclass=Meta):pass
Basic.register(Foreign)
a=DDM([[s.ZZ(2),s.ZZ(0)],[s.ZZ(3),s.ZZ(2)]],(2,2),s.ZZ)
b,t=a.lll_transform()
assert scratch
assert [[str(c) for c in row] for row in scratch[0]]==[['0','0'],['-1/2','0']]
assert [[str(c) for c in row] for row in scratch[1]]==[['2','0'],['0','2']]
'''
    subprocess.run([sys.executable,'-c',script],text=True,capture_output=True,check=True,timeout=30)


def test_removed_helper_reimport_keeps_the_complete_old_path():
    script='''
import symox as s,sys,importlib
from symox.polys.matrices.ddm import DDM
lll=importlib.import_module('symox.polys.matrices.lll')
name='symox.polys.matrices._lll_native'
sys.modules.pop(name)
helper=importlib.import_module(name)
lll._lll_native=helper
a=DDM([[s.ZZ(2),s.ZZ(0)],[s.ZZ(3),s.ZZ(2)]],(2,2),s.ZZ)
assert not helper.prefix(a,s.QQ(3,4),True)
b,t=a.lll_transform()
assert [[int(c) for c in row] for row in b]==[[2,0],[-1,2]]
'''
    subprocess.run([sys.executable,'-c',script],text=True,capture_output=True,check=True,timeout=30)


@pytest.mark.parametrize('where',['alias','class_init','base_init','class_new'])
def test_rank_exception_constructor_observes_original_reduction(monkeypatch,where):
    lll=importlib.import_module('symox.polys.matrices.lll');original=lll.DMRankError
    # This test targets the warm native route: a real tiny public calculation
    # certifies the classification without hiding a first-call speed claim.
    DDM([[s.ZZ(1)]],(1,1),s.ZZ).lll_transform()
    a=DDM([[s.ZZ(1),s.ZZ(2)],[s.ZZ(2),s.ZZ(4)]],(2,2),s.ZZ);events=[]
    assert p.prefix(a,s.QQ(3,4),True)
    def capture():
        frame=inspect.currentframe()
        while frame and frame.f_code.co_name!='_ddm_lll':frame=frame.f_back
        events.append([[int(c) for c in row] for row in frame.f_locals['y']])
    if where=='alias':
        class Custom(original):
            def __init__(self,message):capture();super().__init__(message)
        monkeypatch.setattr(lll,'DMRankError',Custom)
    elif where in ('class_init','base_init'):
        cls=original if where=='class_init' else original.__bases__[0];old=cls.__init__
        def initialize(self,*args):capture();old(self,*args)
        monkeypatch.setattr(cls,'__init__',initialize)
    else:
        old=original.__new__
        def create(cls,*args):capture();return old(cls,*args)
        monkeypatch.setattr(original,'__new__',staticmethod(create))
    def run():
        try:a.lll_transform()
        except original as error:return type(error),str(error)
        raise AssertionError('rank failure disappeared')
    actual=run();trace=events[:];events.clear()
    with monkeypatch.context() as patch:
        patch.setattr(p,'prefix',lambda *args:False);expected=run()
    assert actual==expected and trace==events==[[[1,2],[0,0]]]
    assert [[int(c) for c in row] for row in a]==[[1,2],[2,4]]
