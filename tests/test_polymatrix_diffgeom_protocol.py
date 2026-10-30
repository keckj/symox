"""Official differential-geometry forwarding on exact ground matrices."""
import importlib
import itertools
import os
import subprocess
import sys

import pytest
import symox as s

pm=importlib.import_module('symox.polys.polymatrix')


def matrix(kind):
    x,y,z=s.symbols('dg_x dg_y dg_z')
    field=s.QQ if kind=='QQ'else s.QQ.frac_field(y,z)
    values=[[1,2,3],[3,7,8]]if kind=='QQ'else[[y/(z+1),1,2],[2*y/(z+1),3,5]]
    return pm.PolyMatrix(values,x,ring=field)


def geometry():return importlib.import_module('symox.diffgeom._adapters')


def outcome(call):
    try:return ('value',call())
    except Exception as error:return ('error',type(error),str(error))


def compare_old(monkeypatch,a,operation,events):
    actual=outcome(lambda:getattr(a,operation)());observed=events[:]
    with monkeypatch.context()as patch:
        patch.setattr(pm,'_native_ground',lambda *args:None)
        events.clear();expected=outcome(lambda:getattr(a,operation)());legacy=events[:]
    assert actual==expected and observed==legacy
    return actual,observed


@pytest.mark.parametrize('kind',['QQ','RF'])
@pytest.mark.parametrize('operation',['rref','nullspace'])
def test_official_geometry_keeps_native_ground_values_domains_aliases(monkeypatch,kind,operation):
    geometry();a=matrix(kind);old_items=tuple(a._items);calls=[];native=s._core._polymatrix_ground
    def spy(*args):
        result=native(*args);calls.append(result is not None);return result
    monkeypatch.setattr(s._core,'_polymatrix_ground',spy)
    result=getattr(a,operation)();assert calls and all(calls)
    with monkeypatch.context()as patch:
        patch.setattr(pm,'_native_ground',lambda *args:None)
        assert getattr(a,operation)()==result
    results=result if operation=='nullspace'else[result[0]]
    assert all(v.domain is a.domain and v.gens==a.gens for v in results)
    assert all(v is w for v,w in zip(a._items,old_items))


@pytest.mark.parametrize('order',list(itertools.permutations(('diffgeom','vector','holonomic.recurrence'))))
def test_all_six_optional_import_orders_and_audit_reader(order):
    source='''
import importlib,sys
import symox as s
from symox.polys import polymatrix as p
for name in sys.argv[1:]:importlib.import_module('symox.'+name)
x,y=s.symbols('dg_x dg_y')
matrices=[p.PolyMatrix([[1,2,3],[3,7,8]],x,ring=s.QQ),p.PolyMatrix([[y,1,2],[2*y,3,5]],x,ring=s.QQ.frac_field(y))]
events=[]
def hook(event,args):
 if event=='object.__getattr__':events.append(args[1]);raise RuntimeError('selection audit callback')
sys.addaudithook(hook)
for a in matrices:
 assert p._native_ground(a) is not None
 for operation in ('rref','nullspace'):
  actual=getattr(a,operation)();saved=p._native_ground;p._native_ground=lambda*args:None
  try:expected=getattr(a,operation)()
  finally:p._native_ground=saved
  assert actual==expected
assert events==[],events
'''
    subprocess.run([sys.executable,'-c',source,*order],check=True,env=os.environ.copy())


@pytest.mark.parametrize('name',['__add__','__sub__','__mul__','__truediv__','atoms','rcall'])
@pytest.mark.parametrize('kind',['QQ','RF'])
def test_foreign_public_scalar_wrappers_keep_full_callback_log(monkeypatch,name,kind):
    geometry();a=matrix(kind);events=[];original=getattr(s._core.Expr,name)
    def callback(*args,**kwargs):events.append(name);return original(*args,**kwargs)
    monkeypatch.setattr(s._core.Expr,name,callback)
    assert pm._native_ground(a)is None and events==[]
    _,observed=compare_old(monkeypatch,a,'rref',events)
    assert bool(observed)==(name in ('__sub__','__mul__','__truediv__'))


@pytest.mark.parametrize('mode',['code','defaults','kwdefaults','delegate','tensor','name','empty_cell'])
def test_changed_actual_geometry_function_state_declines_without_callbacks(monkeypatch,mode):
    adapter=geometry();a=matrix('RF');events=[]
    # The geometry forwarder itself, not Expr.__mul__: a vector import after
    # the geometry one puts its own certified forwarder on top of it.
    operation=next(p[0] for p in adapter._SCALAR_DISPATCH_PROOFS if p[-1]=='__mul__')
    def foreign(*args):events.append('foreign');raise RuntimeError('foreign callback')
    if mode=='defaults':monkeypatch.setattr(operation,'__defaults__',(foreign,))
    elif mode=='kwdefaults':monkeypatch.setattr(operation,'__kwdefaults__',{'foreign':foreign})
    elif mode=='code':
        def factory(_Tensor,name,original):
            def changed(self,other):return original(self,other)if name else _Tensor
            return changed
        monkeypatch.setattr(operation,'__code__',factory(None,None,None).__code__)
    else:
        index=operation.__code__.co_freevars.index({'delegate':'original','tensor':'_Tensor','name':'name','empty_cell':'original'}[mode])
        cell=operation.__closure__[index]
        if mode=='empty_cell':
            value=cell.cell_contents;del cell.cell_contents
            try:assert pm._native_ground(a)is None and events==[]
            finally:cell.cell_contents=value
            return
        monkeypatch.setattr(cell,'cell_contents',foreign if mode!='name'else object())
    assert pm._native_ground(a)is None and events==[]


@pytest.mark.parametrize('mode',['holder_alias','holder_len','holder_contains','field_arithmetic','ordinary_field_key','host_reader','module_class','proof_name','rt_len','rt_zip','rt_all'])
def test_geometry_consumed_protocols_retain_errors_and_callback_counts(monkeypatch,mode):
    from types import ModuleType
    adapter=geometry();coordinates=importlib.import_module('symox.diffgeom._coordinates')
    rt=importlib.import_module('symox.integrals.rationaltools')
    a=matrix('RF');events=[]
    def foreign(*args,**kwargs):events.append(mode);raise RuntimeError('foreign geometry protocol')
    if mode=='holder_alias':
        class Holder:
            def __len__(self):return foreign()
        monkeypatch.setattr(s,'_diffgeom_field_atoms',Holder())
    elif mode=='holder_len':monkeypatch.setattr(coordinates._FieldAtomSet,'__len__',foreign,raising=False)
    elif mode=='holder_contains':monkeypatch.setattr(coordinates._FieldAtomSet,'__contains__',foreign)
    elif mode=='field_arithmetic':monkeypatch.setattr(s._diffgeom_callable_type,'__mul__',foreign)
    elif mode=='ordinary_field_key':
        value=s.Symbol('ordinary_geometry_key');coordinates._field_atoms.add(value)
        try:
            assert pm._native_ground(a)is None
            compare_old(monkeypatch,a,'rref',events)
        finally:coordinates._field_atoms.remove(value)
        return
    elif mode=='host_reader':monkeypatch.setattr(s._core,'_host_contains',foreign)
    elif mode=='module_class':
        class Hooks(ModuleType):
            def __getattribute__(self,name):events.append(name);raise RuntimeError('foreign module lookup')
        monkeypatch.setattr(adapter,'__class__',Hooks)
    elif mode=='proof_name':
        class Equal:
            def __eq__(self,other):return foreign()
            def __ne__(self,other):return foreign()
        rows=adapter._SCALAR_DISPATCH_PROOFS
        changed=rows[:1]+(rows[1][:-1]+(Equal(),),)+rows[2:]
        monkeypatch.setattr(adapter,'_SCALAR_DISPATCH_PROOFS',changed)
        monkeypatch.setattr(adapter.install,'_stock_proofs',changed)
    else:monkeypatch.setattr(rt,mode[3:],foreign,raising=False)
    assert pm._native_ground(a)is None and events==[]
    compare_old(monkeypatch,a,'rref',events)


@pytest.mark.parametrize('operation',['rref','nullspace'])
def test_items_subclass_is_not_speculatively_iterated(monkeypatch,operation):
    geometry();a=matrix('RF');events=[]
    key=None
    class Items(list):
        def __iter__(self):events.append('iter');raise RuntimeError('items iter')
        def __len__(self):events.append('len');raise RuntimeError('items len')
    a._items=Items(a._items)
    assert pm._native_ground(a)is None and events==[]
    compare_old(monkeypatch,a,operation,events)


def test_bare_native_field_symbol_keeps_host_payload_and_is_not_a_formal_rf_parameter(monkeypatch):
    from symox.diffgeom.rn import R2_r
    field=R2_r.base_scalar(0);raw=s._core.Expr.__new__(s._core.Expr,field._atom)
    assert field._atom==raw and hash(field._atom)==hash(raw)
    assert s._core._host_contains(field._atom)==s._core._host_contains(raw)==(True,True)
    x=s.Symbol('dg_x');a=pm.PolyMatrix([[raw,1,2],[2*raw,3,5]],x,ring=s.QQ.frac_field(raw))
    assert pm._native_ground(a)is None
    compare_old(monkeypatch,a,'rref',[])


@pytest.mark.parametrize('name',['_ddm_function_metadata','_host_contains'])
@pytest.mark.parametrize('mode',['alias','removed','collision'])
def test_core_reader_origin_and_namespace_collisions_never_call_foreign_protocol(monkeypatch,name,mode):
    geometry();a=matrix('RF');events=[];key=None
    def foreign(*args):events.append('callback');raise RuntimeError('foreign reader')
    if mode=='alias':monkeypatch.setattr(s._core,name,foreign)
    else:
        monkeypatch.delattr(s._core,name)
        if mode=='collision':
            hashed=hash(name)
            class Key:
                def __hash__(self):return hashed
                def __eq__(self,other):events.append('eq');raise RuntimeError('foreign key equality')
            key=Key();s._core.__dict__[key]=None
    try:
        events.clear()
        assert pm._native_ground(a)is None and events==[]
        compare_old(monkeypatch,a,'rref',events)
    finally:
        if key is not None:del s._core.__dict__[key]


@pytest.mark.parametrize('mode',['dict','closure'])
def test_actual_function_metadata_containers_are_checked_before_iteration(monkeypatch,mode):
    from types import FunctionType
    adapter=geometry();a=matrix('RF');events=[]
    row=next(p for p in adapter._SCALAR_DISPATCH_PROOFS if p[-1]=='__mul__')
    operation=row[0]
    if mode=='dict':
        class Attributes(dict):
            def __iter__(self):events.append('iter');raise RuntimeError('foreign function dictionary')
        monkeypatch.setattr(operation,'__dict__',Attributes(operation.__dict__))
    else:
        class Closure(tuple):
            def __iter__(self):events.append('iter');raise RuntimeError('foreign closure iteration')
            def __len__(self):events.append('len');raise RuntimeError('foreign closure length')
            def __getitem__(self,index):events.append('getitem');raise RuntimeError('foreign closure indexing')
        cells=Closure(operation.__closure__)
        replacement=FunctionType(operation.__code__,operation.__globals__,closure=cells)
        replacement.__dict__.update(operation.__dict__)
        states=((replacement,operation.__code__,cells,tuple(c.cell_contents for c in operation.__closure__),None,None),)
        changed=(replacement,row[1],states,*row[3:])
        replacement._symox_diffgeom_proof=changed
        rows=tuple(changed if p is row else p for p in adapter._SCALAR_DISPATCH_PROOFS)
        monkeypatch.setattr(adapter,'_SCALAR_DISPATCH_PROOFS',rows)
        monkeypatch.setattr(adapter.install,'_stock_proofs',rows)
        monkeypatch.setattr(s._core.Expr,'__mul__',replacement)
    assert pm._native_ground(a)is None and events==[]
    compare_old(monkeypatch,a,'rref',events)


def test_public_facade_replacement_cannot_turn_bare_numeric_operands_into_tensors(monkeypatch):
    geometry();tensors=importlib.import_module('symox.diffgeom._tensors')
    a=matrix('QQ');events=[]
    a._items=[s._core.Expr.__new__(s._core.Expr,value)for value in a._items]
    class TensorInteger(tensors._Tensor):
        def __rsub__(self,*args):events.append('tensor');raise RuntimeError('tensor facade')
        def __rmul__(self,*args):events.append('tensor');raise RuntimeError('tensor facade')
    monkeypatch.setattr(s,'Integer',TensorInteger)
    assert pm._native_ground(a)is None and events==[]
    compare_old(monkeypatch,a,'rref',events)


def test_preimport_metadata_wrapper_is_not_blessed_as_a_new_reader():
    source='''
import builtins
original=builtins.__import__;events=[]
def foreign(*args):events.append('reader');raise RuntimeError('foreign metadata')
def importing(name,globals=None,locals=None,fromlist=(),level=0):
 result=original(name,globals,locals,fromlist,level)
 if name=='_core' and globals and globals.get('__package__')=='symox.diffgeom':
  import sys
  sys.modules['symox._core']._ddm_function_metadata=foreign
 return result
builtins.__import__=importing
import symox as s
import symox.diffgeom
from symox.polys import polymatrix as p
x=s.Symbol('dg_x');a=p.PolyMatrix([[1,2,3],[3,7,8]],x,ring=s.QQ)
events.clear()
assert p._native_ground(a)is None
assert events==[]
actual=a.rref();p._native_ground=lambda *args:None
assert a.rref()==actual and events==[]
'''
    subprocess.run([sys.executable,'-c',source],check=True,env=os.environ.copy())


@pytest.mark.parametrize('operation',['rref','nullspace'])
@pytest.mark.parametrize('raises',[False,True])
def test_wild_property_hash_protocol_retains_complete_old_counts_and_errors(monkeypatch,operation,raises):
    from symox.diffgeom.rn import R2_r
    R2_r.base_scalar(0)
    events=[];active=False
    class Predicate:
        def __call__(self,value):return True
        def __hash__(self):
            if active:
                events.append('hash')
                if raises:raise RuntimeError('wild property hash')
            return 173
    w=s.Wild('geometry_wild',properties=(Predicate(),));x=s.Symbol('dg_x')
    a=pm.PolyMatrix([[w+1,1,2],[2*w+2,3,5]],x,ring=s.QQ.frac_field(w))
    active=True
    assert pm._native_ground(a)is None and events==[]
    actual,observed=compare_old(monkeypatch,a,operation,events)
    active=False
    assert observed
    if raises:assert actual[0]=='error'and actual[1]is RuntimeError and observed==['hash']


@pytest.mark.parametrize('kind',['QQ','RF'])
@pytest.mark.parametrize('operation',['rref','nullspace'])
def test_legitimate_rationaltools_reload_keeps_native_forwarding_and_exact_domains(monkeypatch,kind,operation):
    geometry();rt=importlib.import_module('symox.integrals.rationaltools')
    importlib.reload(rt)
    a=matrix(kind);calls=[];native=s._core._polymatrix_ground
    def spy(*args):
        result=native(*args);calls.append(result is not None);return result
    monkeypatch.setattr(s._core,'_polymatrix_ground',spy)
    actual=getattr(a,operation)()
    assert calls and all(calls)
    with monkeypatch.context()as patch:
        patch.setattr(pm,'_native_ground',lambda *args:None)
        expected=getattr(a,operation)()
    assert actual==expected
    values=actual if operation=='nullspace'else[actual[0]]
    assert all(v.domain is a.domain and v.gens==a.gens for v in values)


def test_geometry_state_checks_hold_when_a_later_vector_import_wraps_the_forwarder():
    # Regression: test_diffgeom_completion.py then test_gaps_misc_parity.py
    # import diffgeom and then vector, in one process. The certified vector
    # forwarder then wraps the geometry one, so Expr.__mul__ is no longer the
    # geometry forwarder that the state checks above mutate.
    source='''
import importlib,sys
import pytest
import symox as s
from symox.polys import polymatrix as p
importlib.import_module('symox.diffgeom')
importlib.import_module('symox.vector')
adapter=importlib.import_module('symox.diffgeom._adapters')
operation=next(row[0] for row in adapter._SCALAR_DISPATCH_PROOFS if row[-1]=='__mul__')
top=s._core.Expr.__mul__
assert top is not operation and '_Tensor' not in top.__code__.co_freevars
y,z=s.symbols('dg_y dg_z')
a=p.PolyMatrix([[y/(z+1),1,2],[2*y/(z+1),3,5]],s.Symbol('dg_x'),ring=s.QQ.frac_field(y,z))
assert p._native_ground(a) is not None
sys.exit(pytest.main(['-q','-p','no:cacheprovider',
    sys.argv[1]+'::test_changed_actual_geometry_function_state_declines_without_callbacks']))
'''
    subprocess.run([sys.executable,'-c',source,os.path.abspath(__file__)],check=True,
                   env=os.environ.copy(),timeout=300)


def _insertion_ordered_frozensets():
    # Two equal frozensets whose iteration orders differ. Stock CPython 3.10
    # marshal writes them differently (some distribution builds sort them).
    names=['k%d'%i for i in range(64)]
    for a,b in itertools.combinations(names,2):
        if list(frozenset((a,b)))!=list(frozenset((b,a))):return frozenset((a,b)),frozenset((b,a))
    raise AssertionError('no colliding pair of names')


def test_code_key_ignores_frozenset_iteration_order_only():
    from marshal import dumps
    from symox._code_key import code_key as key
    code_key=lambda value:key(value,dumps)
    first,second=_insertion_ordered_frozensets()
    code=compile('def f(x):\n    return x in {0}\n','m','exec').co_consts[0]
    index=code.co_consts.index(frozenset({0}))
    def with_constant(value):
        return code.replace(co_consts=code.co_consts[:index]+(value,)+code.co_consts[index+1:])
    a,b=with_constant(first),with_constant(second)
    assert list(a.co_consts[index])!=list(b.co_consts[index])
    assert code_key(a)==code_key(b)
    assert code_key(a)!=code_key(with_constant(first|{'k_other'}))
    outer=compile('def g():\n    def f(x):\n        return x in {0}\n    return f\n','m','exec').co_consts[0]
    position=next(i for i,c in enumerate(outer.co_consts) if type(c) is type(code))
    nested=lambda inner:outer.replace(co_consts=outer.co_consts[:position]+(inner,)+outer.co_consts[position+1:])
    assert code_key(nested(a))==code_key(nested(b))!=code_key(outer)
    events=[]
    class Foreign:
        def __eq__(self,other):events.append('eq');return True
        def __hash__(self):events.append('hash');return 0
        def __reduce_ex__(self,protocol):events.append('reduce');return (Foreign,())
    with pytest.raises(ValueError):code_key(code.replace(co_consts=code.co_consts+((Foreign(),),)))
    assert events==[]


def test_reload_reading_a_pyc_of_another_hash_seed_keeps_native_forwarding():
    # Regression: under the parallel suite, rationaltools was compiled from
    # source at import and reloaded from a .pyc written by another process.
    # The frozenset constant of _atan_vector_proof_with_metadata then
    # iterates, and on stock CPython 3.10 marshals, in another order, and the
    # reload authentication declined.
    names=('original','name','_BasisDependent')
    probe=('import itertools,sys\n'
           'n=sys.argv[1:];print(len({tuple(frozenset(p)) for p in itertools.permutations(n)}))')
    env=os.environ.copy()
    for seed in range(200):
        env['PYTHONHASHSEED']=str(seed)
        orders=subprocess.run([sys.executable,'-c',probe,*names],check=True,env=env,
                              capture_output=True,text=True,timeout=60).stdout
        if int(orders)>1:break
    else:pytest.skip('no hash seed reorders the frozenset constant')
    source='''
import importlib,itertools,types
import symox as s
from symox.polys import polymatrix as p
rt=importlib.import_module('symox.integrals.rationaltools')
importlib.reload(rt)
name='_atan_vector_proof_with_metadata';f=rt.__dict__[name];code=f.__code__
index=next(i for i,c in enumerate(code.co_consts) if type(c) is frozenset)
assert set(code.co_consts[index])==set(%r)
reordered=None
for order in itertools.permutations(code.co_consts[index]):
    candidate=code.replace(co_consts=code.co_consts[:index]+(frozenset(order),)+code.co_consts[index+1:])
    if tuple(candidate.co_consts[index])!=tuple(code.co_consts[index]):reordered=candidate;break
assert reordered is not None
# The same reloaded function, as read from a .pyc of another hash seed.
g=types.FunctionType(reordered,f.__globals__,f.__name__,f.__defaults__,f.__closure__)
g.__kwdefaults__=f.__kwdefaults__;g.__qualname__=f.__qualname__
rt.__dict__[name]=g
x,y,z=s.symbols('dg_x dg_y dg_z')
native=s._core._polymatrix_ground;calls=[]
def spy(*args):
    result=native(*args);calls.append(result is not None);return result
s._core._polymatrix_ground=spy
for a in (p.PolyMatrix([[1,2,3],[3,7,8]],x,ring=s.QQ),
          p.PolyMatrix([[y/(z+1),1,2],[2*y/(z+1),3,5]],x,ring=s.QQ.frac_field(y,z))):
    for operation in ('rref','nullspace'):
        calls.clear();actual=getattr(a,operation)()
        assert calls and all(calls),(operation,calls)
        saved=p._native_ground;p._native_ground=lambda *args:None
        try:expected=getattr(a,operation)()
        finally:p._native_ground=saved
        assert actual==expected
'''%(names,)
    subprocess.run([sys.executable,'-c',source],check=True,env=env,timeout=300)
