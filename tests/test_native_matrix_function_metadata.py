"""Public function metadata and complete, non-auditing matrix admission."""
import importlib
import json
import os
from pathlib import Path
import subprocess
import sys
import types

import pytest
import symox as sx

AUDIT_PROBE='"""Fresh-process metadata audits; construct/import/AST preparation precede hook."""\nimport ast,importlib,json,sys,types\nfrom pathlib import Path\nimport symox as s\nnf=importlib.import_module(\'symox.polys.matrices.normalforms\')\ndense=importlib.import_module(\'symox.polys.matrices.dense\')\nfrom symox.polys.matrices.domainmatrix import DomainMatrix\nroute,arm=sys.argv[1:3]\nname={\'hnf\':\'_hermite_normal_form\',\'lu\':\'ddm_ilu\',\'solve\':\'ddm_ilu_solve\'}[route];module=nf if route==\'hnf\' else dense\nnode=next(n for n in ast.parse(Path(module.__file__).read_text()).body if isinstance(n,ast.FunctionDef) and n.name==name)\ndef native_node(n):\n if isinstance(n,ast.Assign):return any(isinstance(t,ast.Name) and t.id==\'native\' for t in n.targets)\n if isinstance(n,ast.If):\n  return (isinstance(n.test,ast.Compare) and isinstance(n.test.left,ast.Name) and n.test.left.id==\'native\') or (isinstance(n.test,ast.Call) and isinstance(n.test.func,ast.Name) and n.test.func.id==\'_native_ilu_solve\')\n return False\nnode.body=[n for n in node.body if not native_node(n)]\ncode=compile(ast.fix_missing_locations(ast.Module(body=[node],type_ignores=[])),module.__file__,\'exec\')\nold=types.FunctionType(next(c for c in code.co_consts if type(c)is types.CodeType),module.__dict__)\nfunction=old if arm==\'old\' else getattr(module,name)\nq=lambda rows:[[s.QQ(v) for v in row] for row in rows]\nif route==\'hnf\':\n a=DomainMatrix.from_list([[1,2],[3,4]],s.ZZ);invoke=lambda:function(a);values=lambda out:[[int(v) for v in row] for row in out.to_list()]\nelif route==\'lu\':\n a=q([[2,1],[4,3]]);invoke=lambda:function(a);values=lambda out:{\'matrix\':[[str(v) for v in row] for row in a],\'result\':out}\nelse:\n x=[[None],[None]];l=q([[1,0],[2,1]]);u=q([[2,1],[0,1]]);b=q([[4],[10]])\n invoke=lambda:function(x,l,u,[],b);values=lambda out:{\'matrix\':[[str(v) for v in row] for row in x],\'result\':out}\ncalls=[]\nkey={\'hnf\':\'_ddm_hnf\',\'lu\':\'_ddm_ilu\',\'solve\':\'_ddm_ilu_solve\'}[route]\nnative=getattr(s._core,key)\ndef trace(*args,**kwargs):\n result=native(*args,**kwargs);calls.append(result is not None if route!=\'solve\' else result);return result\nsetattr(s._core,key,trace)\nif len(sys.argv)>3 and sys.argv[3]==\'unavailable\':s._core._ddm_function_metadata=None\nevents=[]\ndef audit(event,args):\n if event==\'object.__getattr__\' and len(args)>1 and args[1] in (\'__code__\',\'__defaults__\',\'__kwdefaults__\',\'__closure__\',\'__class__\'):\n  events.append(args[1]);raise RuntimeError(\'metadata audit refusal: \'+args[1])\nsys.addaudithook(audit)\ntry:out=invoke()\nexcept Exception as error:result={\'status\':\'raised\',\'type\':type(error).__name__,\'message\':str(error)}\nelse:result={\'status\':\'returned\',\'value\':values(out)}\nassert not any(n==\'sympy\' or n.startswith(\'sympy.\') for n in sys.modules)\nprint(json.dumps({\'route\':route,\'arm\':arm,\'events\':events,\'result\':result,\'reference_imported\':False,\'native_calls\':calls}))\n'


def fresh(route,arm,unavailable=False):
    command=[sys.executable,'-c',AUDIT_PROBE,route,arm]
    if unavailable:command.append('unavailable')
    result=subprocess.run(command,capture_output=True,text=True,env=dict(os.environ,PYTHONPATH=str(Path(sx.__file__).resolve().parents[1])),timeout=45)
    assert result.returncode==0,result.stderr
    return json.loads(result.stdout)


@pytest.mark.parametrize('route',['hnf','lu','solve'])
def test_audit_refusal_has_no_metadata_event_and_matches_original(route):
    native=fresh(route,'native');old=fresh(route,'old')
    assert native['result']==old['result'] and native['result']['status']=='returned'
    assert native['events']==old['events']==[]
    assert native['native_calls']==[True] and old['native_calls']==[]
    assert native['reference_imported'] is old['reference_imported'] is False


@pytest.mark.parametrize('route',['hnf','lu','solve'])
def test_missing_metadata_uses_complete_original_body_without_audit(route):
    result=fresh(route,'native',True);old=fresh(route,'old')
    assert result['result']==old['result'] and result['events']==[]
    assert result['native_calls']==[]


def test_public_metadata_matches_five_real_fields_and_refreshes_without_cache():
    getter=sx._core._ddm_function_metadata
    def factory(captured):
        def function(x=2,*,keyword=3):return captured+x+keyword
        return function
    function=factory(7)
    first=getter(function)
    assert first is not None
    assert all(a is b for a,b in zip(first,(function.__code__,function.__globals__,function.__defaults__,function.__kwdefaults__,function.__closure__)))
    function.__defaults__=(11,);function.__kwdefaults__={'keyword':13}
    def alternate_factory(factor):
        def alternate(x=2,*,keyword=3):return factor*x+keyword
        return alternate
    replacement=alternate_factory(17);function.__code__=replacement.__code__
    again=getter(function)
    assert again is not first and again[0] is not first[0]
    assert all(a is b for a,b in zip(again,(function.__code__,function.__globals__,function.__defaults__,function.__kwdefaults__,function.__closure__)))
    function.__closure__[0].cell_contents=23
    assert getter(function)[4][0].cell_contents==23


def test_foreign_object_refused_without_getattribute():
    class Foreign:
        def __getattribute__(self,name):raise AssertionError(name)
    assert sx._core._ddm_function_metadata(Foreign()) is None
    assert sx._core._ddm_function_metadata(len) is None
    with pytest.raises(TypeError):sx._core._ddm_function_metadata()
    with pytest.raises(TypeError):sx._core._ddm_function_metadata(lambda:None,unexpected=True)


@pytest.mark.parametrize('attribute',['__code__','__defaults__','__kwdefaults__','__closure__','__class__'])
def test_getters_do_not_emit_attribute_audits(attribute):
    script="""import sys,json,symox as s
value=7
def factory(value):
 def function(x=2,*,kw=3):return value+x+kw
 return function
function=factory(value);getter=s._core._ddm_function_metadata;events=[]
def audit(event,args):
 if event=='object.__getattr__' and len(args)>1 and args[1]==sys.argv[1]:
  events.append(args[1]);raise RuntimeError(args[1])
sys.addaudithook(audit)
result=getter(function)
assert result is not None and len(result)==5
print(json.dumps(events))
"""
    p=subprocess.run([sys.executable,'-c',script,attribute],env=dict(os.environ,PYTHONPATH=str(Path(sx.__file__).resolve().parents[1])),capture_output=True,text=True,timeout=30)
    assert p.returncode==0,p.stderr
    assert json.loads(p.stdout)==[]


@pytest.mark.parametrize('helper',['_hnf_native','_lu_native','_lu_solve_native'])
def test_initial_helper_capture_does_not_read_audited_metadata(helper):
    script="""import sys,importlib,json,symox as s
from symox.polys.matrices import normalforms,dense
name='symox.polys.matrices.'+sys.argv[1]
sys.modules.pop(name,None);events=[]
def audit(event,args):
 if event=='object.__getattr__' and len(args)>1 and args[1] in ('__code__','__defaults__','__kwdefaults__','__closure__'):
  events.append(args[1]);raise RuntimeError(args[1])
sys.addaudithook(audit)
module=importlib.import_module(name)
assert module._METADATA_OK
print(json.dumps(events))
"""
    p=subprocess.run([sys.executable,'-c',script,helper],env=dict(os.environ,PYTHONPATH=str(Path(sx.__file__).resolve().parents[1])),capture_output=True,text=True,timeout=30)
    assert p.returncode==0,p.stderr
    assert json.loads(p.stdout)==[]


@pytest.mark.parametrize('route,helper,dispatch',[('hnf','_hnf_native','_native_hnf'),('lu','_lu_native','_native_ilu'),('solve','_lu_solve_native','_native_ilu_solve')])
def test_foreign_metadata_before_helper_import_has_no_speculative_callback(route,helper,dispatch):
    injection="""events_foreign=[]
class Foreign:
 def __getattribute__(self,name):
  events_foreign.append(name);raise AssertionError('speculative metadata access '+name)
s._core._ddm_function_metadata=Foreign()
module_name='symox.polys.matrices.'+%r
sys.modules.pop(module_name,None)
helper=importlib.import_module(module_name)
assert helper._METADATA_OK is False and events_foreign==[]
setattr(nf if %r=='hnf' else dense,%r,helper.try_hnf if %r=='hnf' else helper.try_ilu if %r=='lu' else helper.try_solve)
"""%(helper,route,dispatch,route,route)
    script=AUDIT_PROBE.replace('events=[]\ndef audit',injection+'events=[]\ndef audit').replace("'native_calls':calls}","'native_calls':calls,'foreign_events':events_foreign}")
    p=subprocess.run([sys.executable,'-c',script,route,'native'],env=dict(os.environ,PYTHONPATH=str(Path(sx.__file__).resolve().parents[1])),capture_output=True,text=True,timeout=30)
    assert p.returncode==0,p.stderr
    data=json.loads(p.stdout)
    assert data['result']==fresh(route,'old')['result']
    assert data['native_calls']==data['events']==data['foreign_events']==[]


@pytest.mark.parametrize('route,helper,dispatch',[('hnf','_hnf_native','_native_hnf'),('lu','_lu_native','_native_ilu'),('solve','_lu_solve_native','_native_ilu_solve')])
def test_foreign_core_module_before_helper_import_declines_without_getattr(route,helper,dispatch):
    injection="""from types import ModuleType
core=s._core;original_type=type(core);events_foreign=[]
class ForeignModule(ModuleType):
 def __getattribute__(self,name):
  events_foreign.append(name);raise AssertionError('speculative core access '+name)
core.__class__=ForeignModule
module_name='symox.polys.matrices.'+%r
sys.modules.pop(module_name,None)
helper=importlib.import_module(module_name)
assert helper._METADATA_OK is False and events_foreign==[]
setattr(nf if %r=='hnf' else dense,%r,helper.try_hnf if %r=='hnf' else helper.try_ilu if %r=='lu' else helper.try_solve)
"""%(helper,route,dispatch,route,route)
    script=AUDIT_PROBE.replace('events=[]\ndef audit',injection+'events=[]\ndef audit').replace("assert not any(n=='sympy'", "core.__class__=original_type\nassert not any(n=='sympy'").replace("'native_calls':calls}","'native_calls':calls,'foreign_events':events_foreign}")
    p=subprocess.run([sys.executable,'-c',script,route,'native'],env=dict(os.environ,PYTHONPATH=str(Path(sx.__file__).resolve().parents[1])),capture_output=True,text=True,timeout=30)
    assert p.returncode==0,p.stderr
    data=json.loads(p.stdout)
    assert data['result']==fresh(route,'old')['result']
    assert data['native_calls']==data['events']==data['foreign_events']==[]


@pytest.mark.parametrize('route,helper,dispatch,owner',[('hnf','_hnf_native','_native_hnf','_algebraic_methods'),('hnf','_hnf_native','_native_hnf','_element_print'),('lu','_lu_native','_native_ilu','_algebraic_methods'),('solve','_lu_solve_native','_native_ilu_solve','_algebraic_methods')])
def test_foreign_anchor_module_before_helper_import_has_no_attribute_callbacks(route,helper,dispatch,owner):
    injection="""from types import ModuleType
owner=importlib.import_module('symox.polys.domains.'+%r);original_type=type(owner);events_foreign=[]
class ForeignModule(ModuleType):
 def __getattribute__(self,name):
  events_foreign.append(name);raise AssertionError('speculative anchor access '+name)
owner.__class__=ForeignModule
module_name='symox.polys.matrices.'+%r
sys.modules.pop(module_name,None)
helper=importlib.import_module(module_name)
assert helper._METADATA_OK is False and events_foreign==[]
setattr(nf if %r=='hnf' else dense,%r,helper.try_hnf if %r=='hnf' else helper.try_ilu if %r=='lu' else helper.try_solve)
"""%(owner,helper,route,dispatch,route,route)
    script=AUDIT_PROBE.replace('events=[]\ndef audit',injection+'events=[]\ndef audit').replace("assert not any(n=='sympy'", "owner.__class__=original_type\nassert not any(n=='sympy'").replace("'native_calls':calls}","'native_calls':calls,'foreign_events':events_foreign}")
    p=subprocess.run([sys.executable,'-c',script,route,'native'],env=dict(os.environ,PYTHONPATH=str(Path(sx.__file__).resolve().parents[1])),capture_output=True,text=True,timeout=30)
    assert p.returncode==0,p.stderr
    data=json.loads(p.stdout)
    assert data['result']==fresh(route,'old')['result']
    assert data['native_calls']==data['events']==data['foreign_events']==[]


@pytest.mark.parametrize('route,helper,dispatch',[('hnf','_hnf_native','_native_hnf'),('lu','_lu_native','_native_ilu'),('solve','_lu_solve_native','_native_ilu_solve')])
@pytest.mark.parametrize('anchor',['Domain','DomainElement'])
def test_foreign_domain_anchor_before_helper_import_has_no_attribute_callbacks(route,helper,dispatch,anchor):
    injection="""events_foreign=[]
class Foreign:
 def __getattribute__(self,name):
  events_foreign.append(name);raise AssertionError('speculative domain anchor access '+name)
setattr(s._core,%r,Foreign())
module_name='symox.polys.matrices.'+%r
sys.modules.pop(module_name,None)
helper=importlib.import_module(module_name)
assert events_foreign==[]
setattr(nf if %r=='hnf' else dense,%r,helper.try_hnf if %r=='hnf' else helper.try_ilu if %r=='lu' else helper.try_solve)
"""%(anchor,helper,route,dispatch,route,route)
    script=AUDIT_PROBE.replace('events=[]\ndef audit',injection+'events=[]\ndef audit').replace("'native_calls':calls}","'native_calls':calls,'foreign_events':events_foreign}")
    p=subprocess.run([sys.executable,'-c',script,route,'native'],env=dict(os.environ,PYTHONPATH=str(Path(sx.__file__).resolve().parents[1])),capture_output=True,text=True,timeout=30)
    assert p.returncode==0,p.stderr
    data=json.loads(p.stdout)
    assert data['result']==fresh(route,'old')['result']
    assert data['native_calls']==data['events']==data['foreign_events']==[]


@pytest.mark.parametrize('route,helper,dispatch',[('hnf','_hnf_native','_native_hnf'),('lu','_lu_native','_native_ilu'),('solve','_lu_solve_native','_native_ilu_solve')])
@pytest.mark.parametrize('phase',['preimport','call'])
def test_foreign_core_key_never_invokes_equality_in_metadata_lookup(route,helper,dispatch,phase):
    injection="""events_foreign=[]
namespace=s._core.__dict__;key_hash=hash('_ddm_function_metadata')
namespace.pop('_ddm_function_metadata')
class ForeignKey:
 def __hash__(self):return key_hash
 def __eq__(self,other):
  events_foreign.append('equality');raise RuntimeError('foreign core key equality')
namespace[ForeignKey()]=17
if %r=='preimport':
 module_name='symox.polys.matrices.'+%r
 sys.modules.pop(module_name,None)
 helper=importlib.import_module(module_name)
 assert helper._METADATA_OK is False and events_foreign==[]
 setattr(nf if %r=='hnf' else dense,%r,helper.try_hnf if %r=='hnf' else helper.try_ilu if %r=='lu' else helper.try_solve)
"""%(phase,helper,route,dispatch,route,route)
    script=AUDIT_PROBE.replace('events=[]\ndef audit',injection+'events=[]\ndef audit').replace("'native_calls':calls}","'native_calls':calls,'foreign_events':events_foreign}")
    p=subprocess.run([sys.executable,'-c',script,route,'native'],env=dict(os.environ,PYTHONPATH=str(Path(sx.__file__).resolve().parents[1])),capture_output=True,text=True,timeout=30)
    assert p.returncode==0,p.stderr
    data=json.loads(p.stdout)
    assert data['result']==fresh(route,'old')['result']
    assert data['native_calls']==data['events']==data['foreign_events']==[]
