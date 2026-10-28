"""Identity proofs must not hide changed protocol dependencies or callbacks."""
import importlib
import sys
import types

import pytest

rt = importlib.import_module('symox.integrals.rationaltools')


def legacy(states):
    from types import FunctionType
    if rt._ATAN_TYPE(states) is not tuple:return False
    try:
        for state in states:
            if rt._ATAN_TYPE(state) is not tuple or len(state)!=6:return False
            function,code,closure,contents,defaults,kwdefaults=state
            if (rt._ATAN_TYPE(function) is not FunctionType or rt._ATAN_TYPE(contents) is not tuple
                    or function.__code__ is not code or function.__closure__ is not closure
                    or function.__defaults__ is not defaults or function.__kwdefaults__ is not kwdefaults
                    or len(closure or ())!=len(contents)
                    or not rt._ATAN_ALL(cell.cell_contents is value for cell,value in zip(closure or (),contents))):return False
    except (AttributeError,ValueError,TypeError):return False
    return True


def function():
    value=object()
    def f(x=object(), *, y=object()):return value,x,y
    return f


def snapshot(f):
    closure=f.__closure__
    return ((f,f.__code__,closure,tuple(c.cell_contents for c in closure or ()),f.__defaults__,f.__kwdefaults__),)


@pytest.mark.parametrize('change', ['none','code','defaults','kwdefaults','cell','empty'])
def test_function_identity_mutations(change):
    f=function(); states=snapshot(f)
    if change=='code':
        def factory():
            value=object()
            def other(x=0, *, y=0):return x,y,value
            return other
        f.__code__=factory().__code__
    elif change=='defaults':f.__defaults__=(object(),)
    elif change=='kwdefaults':f.__kwdefaults__={'y':object()}
    elif change=='cell':f.__closure__[0].cell_contents=object()
    elif change=='empty':del f.__closure__[0].cell_contents
    assert rt._atan_states_match(states) is legacy(states) is (change=='none')


@pytest.mark.parametrize('states', [None,[],(),((),),((object(),None,None,(),None,None),)])
def test_malformed_snapshots(states):
    assert rt._atan_states_match(states) is legacy(states)


def test_empty_closure_and_shared_mutable_contents():
    def f(x=[]):return x
    states=snapshot(f)
    f.__defaults__[0].append(1)
    assert rt._atan_states_match(states) is legacy(states) is True
    f.__defaults__=([],)
    assert rt._atan_states_match(states) is legacy(states) is False


def test_exact_tuple_and_no_foreign_getattr():
    class T(tuple):pass
    class Host:
        def __getattribute__(self,name):raise AssertionError(name)
    f=function(); states=snapshot(f)
    for value in [T(states),((Host(),None,None,(),None,None),), (tuple([*states[0][:3],T(states[0][3]),*states[0][4:]]),)]:
        assert rt._atan_states_match(value) is legacy(value) is False


@pytest.mark.parametrize('name',['_ATAN_TYPE','_ATAN_ALL'])
def test_changed_dependency_runs_original_callbacks(monkeypatch,name):
    calls=[]; original=getattr(rt,name)
    def callback(value):
        calls.append(value)
        return original(value)
    monkeypatch.setattr(rt,name,callback)
    states=snapshot(function())
    assert rt._atan_states_match(states)
    native_count=len(calls);calls.clear()
    assert legacy(states)
    assert len(calls)==native_count


@pytest.mark.parametrize('error',[KeyboardInterrupt,SystemExit,MemoryError])
def test_changed_all_interrupt_is_not_swallowed(monkeypatch,error):
    def interrupted(values):raise error('observable all')
    monkeypatch.setattr(rt,'_ATAN_ALL',interrupted)
    states=snapshot(function())
    with pytest.raises(error,match='observable all'):rt._atan_states_match(states)
    with pytest.raises(error,match='observable all'):legacy(states)


def test_function_type_alias_is_observed(monkeypatch):
    states=snapshot(function())
    monkeypatch.setattr(types,'FunctionType',object)
    assert rt._atan_states_match(states) is legacy(states) is False


def test_replaced_extension_helper_cannot_bless_stale_snapshot(monkeypatch):
    f=function();states=snapshot(f);f.__defaults__=(object(),)
    monkeypatch.setattr(rt._rational_core,'_atan_states_match',lambda *args:True)
    assert rt._atan_states_match(states) is legacy(states) is False


@pytest.mark.parametrize('name',['len','zip','tuple'])
def test_changed_module_builtins_use_legacy_loop(monkeypatch,name):
    states=snapshot(function());calls=[]
    import builtins
    original=getattr(builtins,name)
    def callback(*args):
        calls.append(args)
        return original(*args)
    monkeypatch.setattr(rt,name,callback,raising=False)
    if name=='tuple':
        # The original exact-type check observes the replaced alias immediately.
        assert rt._atan_states_match(states) is False
        assert calls==[]
    else:
        assert rt._atan_states_match(states)
        assert calls


def test_canonical_snapshot_selects_extension_without_python_generator():
    states=snapshot(function());seen=[]
    def profile(frame,event,arg):
        if event=='c_call' and arg is rt._ATAN_NATIVE_STATES_MATCH:seen.append('native')
        if event=='call' and frame.f_code.co_name=='<genexpr>' and frame.f_code.co_filename==rt.__file__:seen.append('generator')
    previous=sys.getprofile()
    try:
        sys.setprofile(profile)
        assert rt._atan_states_match(states)
    finally:sys.setprofile(previous)
    assert seen==['native']


def test_empty_cell_followed_by_valid_snapshot_recovers():
    f=function(); states=snapshot(f)
    del f.__closure__[0].cell_contents
    assert rt._atan_states_match(states) is False
    f.__closure__[0].cell_contents=states[0][3][0]
    assert rt._atan_states_match(states) is True


def test_native_identity_loop_propagates_pending_signal():
    import signal
    if not hasattr(signal,'setitimer'):pytest.skip('requires interval timer')
    states=snapshot(function())*1_000_000
    dependencies=(type,all,len,zip,tuple)
    helper=rt._ATAN_NATIVE_STATES_MATCH
    assert helper(states[:1],types.FunctionType,dependencies) is True
    def interrupted(signum,frame):raise KeyboardInterrupt('native proof interrupted')
    previous=signal.signal(signal.SIGALRM,interrupted)
    try:
        signal.setitimer(signal.ITIMER_REAL,0.005)
        with pytest.raises(KeyboardInterrupt,match='native proof interrupted'):
            helper(states,types.FunctionType,dependencies)
    finally:
        signal.setitimer(signal.ITIMER_REAL,0)
        signal.signal(signal.SIGALRM,previous)


@pytest.mark.parametrize('name',['all','len'])
@pytest.mark.parametrize('owner_kind',['none','foreign-module','hostile-module'])
def test_foreign_named_c_functions_at_bootstrap_keep_legacy_protocol(name,owner_kind):
    import subprocess
    script=r'''
import builtins,ctypes,importlib.abc,importlib.util,sys,types
class NoReference(importlib.abc.MetaPathFinder):
 def find_spec(self,fullname,path=None,target=None):
  if fullname=='sympy' or fullname.startswith('sympy.'):raise AssertionError('runtime SymPy')
sys.meta_path.insert(0,NoReference())
import numpy
name,owner_kind,path=sys.argv[1:]
observed=[]
class HostileModule(types.ModuleType):
 def __getattribute__(self,name):
  observed.append(name)
  raise KeyboardInterrupt('module metadata callback')
owner=(None if owner_kind=='none' else
       types.ModuleType('builtins') if owner_kind=='foreign-module' else HostileModule('builtins'))
# Model an ordinary foreign extension C function with the same public name.
# Copy a real CPython C implementation; no code patch or callback trampoline.
class MethodDef(ctypes.Structure):
 _fields_=[('name',ctypes.c_char_p),('method',ctypes.c_void_p),('flags',ctypes.c_int),('doc',ctypes.c_char_p)]
get_function=ctypes.pythonapi.PyCFunction_GetFunction
get_function.argtypes=[ctypes.py_object];get_function.restype=ctypes.c_void_p
get_flags=ctypes.pythonapi.PyCFunction_GetFlags
get_flags.argtypes=[ctypes.py_object];get_flags.restype=ctypes.c_int
create=ctypes.pythonapi.PyCFunction_NewEx
create.argtypes=[ctypes.POINTER(MethodDef),ctypes.c_void_p,ctypes.py_object];create.restype=ctypes.py_object
source=builtins.any if name=='all' else builtins.len
# The definition and its name storage remain alive for the C function's lifetime.
definition=MethodDef(name.encode(),get_function(source),get_flags(source),None)
# Borrow the retained owner's PyObject pointer without ctypes attribute lookup.
foreign=create(ctypes.byref(definition),None if owner is None else id(owner),'builtins')
assert type(foreign) is types.BuiltinFunctionType and foreign.__name__==name
original=getattr(builtins,name)
try:
 setattr(builtins,name,foreign)
 spec=importlib.util.spec_from_file_location('_core',path)
 core=importlib.util.module_from_spec(spec)
 spec.loader.exec_module(core)
finally:setattr(builtins,name,original)
f=lambda:None
states=((f,f.__code__,f.__closure__,(),f.__defaults__,f.__kwdefaults__),)
all_=foreign if name=='all' else all
len_=foreign if name=='len' else len
# Before this fix the native helper returns True, bypassing the foreign all.
assert core._atan_states_match(states,types.FunctionType,(type,all_,len_,zip,tuple)) is None
# The complete old loop continues to observe the captured function.
assert (bool(all_(iter(()))) if name=='all' else len_(states[0])==6) is (name=='len')
assert observed==[]
assert not any(k=='sympy' or k.startswith('sympy.') for k in sys.modules)
'''
    subprocess.run([sys.executable,'-c',script,name,owner_kind,rt._rational_core.__file__],check=True)
