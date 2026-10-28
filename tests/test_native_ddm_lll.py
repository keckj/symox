"""Exact public LLL: complete reference fixtures and independent certificates."""
from fractions import Fraction
import importlib
from pathlib import Path
import subprocess
import sys
import json
import pytest
import symox as s
from symox.polys.matrices.ddm import DDM
from symox.polys.matrices.domainmatrix import DomainMatrix
from symox.polys.matrices.exceptions import DMRankError,DMShapeError,DMValueError

lll=importlib.import_module('symox.polys.matrices.lll')
protocol=importlib.import_module('symox.polys.matrices._lll_native')
methods=importlib.import_module('symox.polys.domains._methods')
ddm=importlib.import_module('symox.polys.matrices.ddm')
fixtures=[line.split('\t') for line in (Path(__file__).parents[1]/'crates/symox-core/tests/fixtures/lll_sympy_1_14.tsv').read_text().splitlines()]


def matrix(rows,n=None):
    n=len(rows[0]) if rows and n is None else (n or 0)
    return DDM([[s.ZZ(v) for v in row] for row in rows],(len(rows),n),s.ZZ)


def ints(value):return [[int(c) for c in row] for row in value]


def certificate(a,b,t,delta):
    m=len(a);n=len(a[0]) if a else 0
    assert [[sum(t[i][k]*a[k][j] for k in range(m)) for j in range(n)] for i in range(m)]==b
    u=[[Fraction(c) for c in row] for row in t];det=Fraction(1)
    for i in range(m):
        pivot=next(j for j in range(i,m) if u[j][i])
        if pivot!=i:u[i],u[pivot]=u[pivot],u[i];det=-det
        v=u[i][i];det*=v
        for j in range(i+1,m):
            c=u[j][i]/v
            for k in range(i+1,m):u[j][k]-=c*u[i][k]
    assert abs(det)==1
    stars=[];norms=[];mus=[]
    for row in b:
        origin=list(map(Fraction,row));star=origin[:];mu=[]
        for prev,norm in zip(stars,norms):
            c=sum(a*b for a,b in zip(origin,prev))/norm;mu.append(c)
            star=[a-c*b for a,b in zip(star,prev)]
        stars.append(star);norms.append(sum(v*v for v in star));mus.append(mu)
        assert all(abs(c)<=Fraction(1,2) for c in mu)
    assert all(norms[i]>=(delta-mus[i][i-1]**2)*norms[i-1] for i in range(1,m))


@pytest.mark.parametrize('case',range(len(fixtures)))
def test_all_recorded_exact_bases_transforms_and_certificates(case,monkeypatch):
    f=fixtures[case];m,n=int(f[0]),int(f[1]);dn,dd=int(f[2]),int(f[3])
    flat=list(map(int,f[5].split(','))) if f[5] else []
    a=matrix([flat[i*n:(i+1)*n] for i in range(m)],n)
    original=ints(a);rows=[id(row) for row in a];cells=[[id(c) for c in row] for row in a]
    delta=s.QQ(dn,dd);calls=[];native=s._core._ddm_lll_native
    def spy(*args):r=native(*args);calls.append(r is not None);return r
    monkeypatch.setattr(s._core,'_ddm_lll_native',spy)
    error={'rank':DMRankError,'shape':DMShapeError,'delta':DMValueError}.get(f[4])
    if error:
        with pytest.raises(error):a.lll_transform(delta)
    else:
        b,t=a.lll_transform(delta)
        expected_b=list(map(int,f[6].split(','))) if f[6] else []
        expected_t=list(map(int,f[7].split(','))) if f[7] else []
        assert [c for row in ints(b) for c in row]==expected_b
        assert [c for row in ints(t) for c in row]==expected_t
        assert ints(a.lll(delta))==ints(b)
        certificate(original,ints(b),ints(t),Fraction(dn,dd))
        if m:assert calls and all(calls)
        assert type(b) is DDM and type(t) is DDM and b.domain is s.ZZ and t.domain is s.ZZ
    assert ints(a)==original and [id(row) for row in a]==rows
    assert [[id(c) for c in row] for row in a]==cells


@pytest.mark.parametrize('entry',['ddm','sdm','dm_dense','dm_sparse'])
def test_actual_public_domainmatrix_and_sparse_consumers(entry,monkeypatch):
    from symox.polys.matrices.sdm import SDM
    a=matrix([[2,0,1],[3,2,0]])
    obj={'ddm':a,'sdm':a.to_sdm(),'dm_dense':DomainMatrix.from_rep(a),'dm_sparse':DomainMatrix.from_rep(a.to_sdm())}[entry]
    calls=[];native=s._core._ddm_lll_native
    def spy(*args):r=native(*args);calls.append(r is not None);return r
    monkeypatch.setattr(s._core,'_ddm_lll_native',spy)
    b,t=obj.lll_transform()
    # Preserve the first actual public result; also check its warmed route.
    if not calls:assert obj.lll_transform()==(b,t)
    if isinstance(b,DomainMatrix):b,t=b.to_list(),t.to_list()
    elif isinstance(b,SDM):b,t=b.to_list(),t.to_list()
    assert ints(b)==[[2,0,1],[1,2,-1]] and ints(t)==[[1,0],[-1,1]] and calls and all(calls)


@pytest.mark.parametrize('hook',['copy','zeros','eye','constructor','elements','convert_from','convert','call',
                                  'floor','add','sub','mul','div','abs','float','parent','ge','attribute'])
def test_custom_observable_protocols_keep_the_entire_original_trace(monkeypatch,hook):
    a=matrix([[2,0,1],[3,2,0]]);events=[]
    if hook in ('copy','constructor'):
        name='copy' if hook=='copy' else '__init__';old=getattr(DDM,name)
        def callback(*args,**kw):events.append(hook);return old(*args,**kw)
        monkeypatch.setattr(DDM,name,callback)
    elif hook in ('zeros','eye'):
        old=getattr(DDM,hook)
        def callback(cls,*args,**kw):events.append(hook);return old(*args,**kw)
        monkeypatch.setattr(DDM,hook,classmethod(callback))
    elif hook=='elements':
        old=ddm._native_elements
        def callback(*args):events.append(hook);return old(*args)
        monkeypatch.setattr(ddm,'_native_elements',callback)
    elif hook in ('convert_from','convert','call','attribute'):
        name={'call':'__call__','attribute':'__getattribute__'}.get(hook,hook);old=getattr(type(s.QQ),name)
        def callback(*args,**kw):events.append(hook);return old(*args,**kw)
        monkeypatch.setattr(type(s.QQ),name,callback)
    elif hook=='floor':
        old=lll.mfloor
        def callback(*args):events.append(hook);return old(*args)
        monkeypatch.setattr(lll,'mfloor',callback)
    else:
        name={'add':'__add__','sub':'__sub__','mul':'__mul__','div':'__truediv__','abs':'__abs__',
              'float':'__float__','parent':'parent','ge':'__ge__'}[hook];old=getattr(s._core.DomainElement,name)
        def callback(*args,**kw):events.append(hook);return old(*args,**kw)
        monkeypatch.setattr(s._core.DomainElement,name,callback)
    b,t=a.lll_transform();trace=events[:];events.clear()
    with monkeypatch.context() as patch:
        patch.setattr(protocol,'prefix',lambda *args:False)
        expected=a.lll_transform()
    assert (ints(b),ints(t))==tuple(ints(v) for v in expected) and trace==events and trace


@pytest.mark.parametrize('what',['source','value'])
def test_temporary_conversion_hook_cannot_lose_its_actual_first_result(monkeypatch,what):
    original=type(s.QQ).convert_from
    def run(enable):
        a=matrix([[2,0,1],[3,2,0]]);events=[]
        def convert(self,value,base):
            events.append('convert');type(s.QQ).convert_from=original
            result=original(self,value,base)
            return result+s.QQ.one if what=='value' else result
        with monkeypatch.context() as patch:
            patch.setattr(type(s.QQ),'convert_from',convert)
            if not enable:patch.setattr(protocol,'prefix',lambda *args:False)
            result=a.lll_transform()
        return tuple(ints(v) for v in result),events
    assert run(True)==run(False)


def test_irrelevant_domain_sum_hook_does_not_disable_the_native_route(monkeypatch):
    def forbidden(*args):raise AssertionError('LLL does not call Domain.sum')
    monkeypatch.setattr(type(s.QQ),'sum',forbidden)
    a=matrix([[2,0],[3,2]]);calls=[];native=s._core._ddm_lll_native
    def spy(*args):r=native(*args);calls.append(r is not None);return r
    monkeypatch.setattr(s._core,'_ddm_lll_native',spy)
    a.lll_transform()
    if not calls:a.lll_transform()
    assert calls and all(calls)


@pytest.mark.parametrize('data',[([[9,0],[0,1]]),([[1,0],[2,1]]),([[1,0],[0,1]])])
def test_native_publication_preserves_untouched_rows_and_cell_objects(data):
    a=matrix(data);y=a.copy();t=a.eye(len(a),s.ZZ)
    y_rows=list(y);t_rows=list(t);y_cells=[list(row) for row in y];t_cells=[list(row) for row in t]
    result=s._core._ddm_lll_native(y,len(a[0]),s.QQ(3,4),t,[s.QQ.convert_from(v,s.ZZ) for v in y[0]],DDM,DMRankError)
    assert result is not None
    basis,transform=result
    if data[0][0]==9:sources=[1,0]
    elif data[1][0]==2:sources=[0,None]
    else:sources=[0,1]
    for i,source in enumerate(sources):
        if source is None:
            assert all(basis[i] is not row for row in y_rows)
            assert all(transform[i] is not row for row in t_rows)
        else:
            assert basis[i] is y_rows[source] and transform[i] is t_rows[source]
            assert all(v is old for v,old in zip(basis[i],y_cells[source]))
            assert all(v is old for v,old in zip(transform[i],t_cells[source]))
    assert list(y)==y_rows and list(t)==t_rows


def test_actual_changed_conversion_value_declines_without_replay():
    a=matrix([[1,0],[2,1]]);before=ints(a)
    assert s._core._ddm_lll_native(a,2,s.QQ(3,4),None,[s.QQ(2),s.QQ(0)],DDM,DMRankError) is None
    assert ints(a)==before


@pytest.mark.parametrize('where',['mfloor','convert_from','copy'])
@pytest.mark.parametrize('error',[ValueError,KeyboardInterrupt])
def test_callback_exceptions_keep_the_original_error_and_trace(monkeypatch,where,error):
    a=matrix([[2,0],[3,2]]);events=[]
    def callback(*args,**kw):events.append(where);raise error('LLL host failure')
    if where=='mfloor':monkeypatch.setattr(lll,'mfloor',callback)
    elif where=='convert_from':monkeypatch.setattr(type(s.QQ),'convert_from',callback)
    else:monkeypatch.setattr(DDM,'copy',callback)
    def run():
        try:a.lll_transform()
        except BaseException as exc:return type(exc),str(exc)
        raise AssertionError('host error lost')
    actual=run();trace=events[:];events.clear()
    with monkeypatch.context() as patch:
        patch.setattr(protocol,'prefix',lambda *args:False);expected=run()
    assert actual==expected==(error,'LLL host failure') and events==trace==[where]


@pytest.mark.parametrize('target',['native_call','div_delegate','convert_code'])
def test_forged_stock_code_closure_and_delegate_preserve_the_old_path(monkeypatch,target):
    a=matrix([[2,0],[3,2]]);events=[]
    if target=='native_call':
        old=methods._native_call
        def callback(*args):events.append(target);return old(*args)
        monkeypatch.setattr(methods,'_native_call',callback)
    elif target=='div_delegate':
        f=s._core.DomainElement.__truediv__;cell=f.__closure__[f.__code__.co_freevars.index('native')];old=cell.cell_contents
        def callback(*args):events.append(target);return old(*args)
        cell.cell_contents=callback
    else:
        f=methods.convert_from;old=f.__code__
        def callback(self,value,base):raise ValueError('forged conversion code')
        monkeypatch.setattr(f,'__code__',callback.__code__)
    try:
        assert not protocol.prefix(a,s.QQ(3,4),True)
        events.clear()  # Exclude this test's explicit QQ construction above.
        def run():
            try:return tuple(ints(v) for v in a.lll_transform())
            except Exception as e:return type(e),str(e)
        actual=run();trace=events[:];events.clear()
        with monkeypatch.context() as patch:
            patch.setattr(protocol,'prefix',lambda *args:False);expected=run()
        assert actual==expected and events==trace
        if target!='convert_code':assert trace
    finally:
        if target=='div_delegate':cell.cell_contents=old


@pytest.mark.parametrize('delta',[s.QQ(1,4),s.QQ(1),s.QQ(-1),False,True,0,1])
@pytest.mark.parametrize('bad',['shape','domain','both'])
def test_validation_order_is_unchanged(delta,bad,monkeypatch):
    a=matrix([[1],[2]]) if bad in ('shape','both') else matrix([[1,2]])
    if bad in ('domain','both'):a.domain=s.QQ
    def run():
        try:a.lll_transform(delta)
        except Exception as e:return type(e),str(e)
        raise AssertionError('invalid case accepted')
    actual=run()
    with monkeypatch.context() as patch:
        patch.setattr(protocol,'prefix',lambda *args:False);expected=run()
    assert actual==expected


@pytest.mark.parametrize('module',['symox.stats','symox.vector','symox.holonomic.recurrence','reload_lll','reload_protocol','reload_ddm'])
def test_official_imports_and_reload_keep_the_stock_public_route_native(module):
    script='''
import symox as s,importlib
from symox.polys.matrices.ddm import DDM
module=MODULE
if module=='reload_lll':importlib.reload(importlib.import_module('symox.polys.matrices.lll'))
elif module=='reload_protocol':importlib.reload(importlib.import_module('symox.polys.matrices._lll_native'))
elif module=='reload_ddm':DDM=importlib.reload(importlib.import_module('symox.polys.matrices.ddm')).DDM
else:importlib.import_module(module)
a=DDM([[s.ZZ(2),s.ZZ(0)],[s.ZZ(3),s.ZZ(2)]],(2,2),s.ZZ)
calls=[];native=s._core._ddm_lll_native
def spy(*args):
 r=native(*args);calls.append(r is not None);return r
s._core._ddm_lll_native=spy
b,t=a.lll_transform()
assert [[int(c) for c in row] for row in b]==[[2,0],[-1,2]]
# A cold first call retains the old complete calculation. The next actual
# public call is eligible; this assertion describes only that warmed route.
if not calls:assert a.lll_transform()==(b,t)
assert calls and all(calls),calls
'''.replace('MODULE',repr(module))
    subprocess.run([sys.executable,'-c',script],text=True,capture_output=True,check=True,timeout=30)


def test_runtime_no_reference_imports_and_honest_cold_then_warm_routes():
    script='''
import sys
class Block:
 def find_spec(self,name,*args):
  if name=='sympy' or name.startswith('sympy.'):raise AssertionError('reference imported')
sys.meta_path.insert(0,Block())
import symox as s
from symox.polys.matrices.ddm import DDM
from symox.polys.matrices.domainmatrix import DomainMatrix
calls=[];native=s._core._ddm_lll_native
def spy(*args):
 r=native(*args);calls.append(r is not None);return r
s._core._ddm_lll_native=spy
a=DDM([[s.ZZ(2),s.ZZ(0)],[s.ZZ(3),s.ZZ(2)]],(2,2),s.ZZ)
p=__import__('importlib').import_module('symox.polys.matrices._lll_native')
initial=p.prefix(a,s.QQ(3,4),True)
dm=DomainMatrix.from_rep(a)
b,t=dm.lll_transform()
assert [[int(c) for c in row] for row in b.to_list()]==[[2,0],[-1,2]]
assert bool(calls)==initial
assert dm.lll_transform()==(b,t)
assert calls and all(calls)
'''
    subprocess.run([sys.executable,'-c',script],text=True,capture_output=True,check=True,timeout=30)


def test_native_interruption_is_not_a_declined_optimization_or_partial_result():
    import signal,time
    script='''
import symox as s,random,signal
from symox.polys.matrices.ddm import DDM
class Stop(BaseException):pass
def stop(*args):raise Stop()
signal.signal(signal.SIGINT,stop)
rng=random.Random(170319);n=22
rows=[[s.ZZ(rng.getrandbits(512)-(1<<511)) for _ in range(n+2)] for _ in range(n)]
a=DDM(rows,(n,n+2),s.ZZ);before=[[int(c) for c in row] for row in a]
native=s._core._ddm_lll_native
announced=False
def spy(*args):
 global announced
 if not announced:announced=True;print('native',flush=True)
 return native(*args)
# This interruption test targets the warm native kernel, not cold routing.
DDM([[s.ZZ(1)]],(1,1),s.ZZ).lll_transform()
s._core._ddm_lll_native=spy
try:a.lll_transform()
except Stop:
 assert announced and [[int(c) for c in row] for row in a]==before
 print('interrupted',flush=True)
else:raise AssertionError('interruption lost')
'''
    process=subprocess.Popen([sys.executable,'-c',script],stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    try:
        assert process.stdout.readline().strip()=='native'
        time.sleep(.05);process.send_signal(signal.SIGINT)
        out,error=process.communicate(timeout=30)
        assert process.returncode==0 and out.strip()=='interrupted',error
    finally:
        if process.poll() is None:process.kill();process.wait()


@pytest.mark.parametrize('name',['shape','domain','__new__','of_type','elements_alias'])
def test_constructor_and_metadata_descriptors_keep_old_callbacks(monkeypatch,name):
    a=matrix([[2,0],[3,2]]);events=[]
    if name in ('shape','domain'):
        def get(self):events.append(name);return self.__dict__[name]
        def put(self,value):self.__dict__[name]=value
        monkeypatch.setattr(DDM,name,property(get,put),raising=False)
    elif name=='__new__':
        def new(cls,*args,**kwargs):events.append(name);return list.__new__(cls)
        monkeypatch.setattr(DDM,name,staticmethod(new))
    elif name=='of_type':
        old=type(s.ZZ).of_type
        def of_type(self,value):events.append(name);return old(self,value)
        monkeypatch.setattr(type(s.ZZ),'of_type',of_type)
    else:
        old=ddm.DDM
        def build(*args,**kwargs):events.append(name);return old(*args,**kwargs)
        monkeypatch.setattr(ddm,'DDM',build)
    actual=tuple(ints(v) for v in a.lll_transform());trace=events[:];events.clear()
    with monkeypatch.context() as patch:
        patch.setattr(protocol,'prefix',lambda *args:False)
        expected=tuple(ints(v) for v in a.lll_transform())
    assert actual==expected and trace==events and trace


def test_temporary_constructor_hook_does_not_hide_exposed_scratch_mutations(monkeypatch):
    descriptor=DDM.__dict__['zeros'];old=DDM.zeros
    def run(enable):
        a=matrix([[2,0],[3,2]]);scratch=[];calls=[];native=s._core._ddm_lll_native
        def zeros(cls,*args):
            DDM.zeros=descriptor
            result=old(*args);scratch.append(result);return result
        def spy(*args):calls.append(True);return native(*args)
        with monkeypatch.context() as patch:
            patch.setattr(DDM,'zeros',classmethod(zeros));patch.setattr(s._core,'_ddm_lll_native',spy)
            if not enable:patch.setattr(protocol,'prefix',lambda *args:False)
            result=tuple(ints(v) for v in a.lll_transform())
        return result,[ints(v) for v in scratch],calls
    assert run(True)==run(False)


def test_reload_code_authentication_never_probes_arbitrary_constants():
    events=[]
    class Foreign:
        def __eq__(self,other):events.append('eq');raise AssertionError('metadata comparison callback')
    code=DDM.copy.__code__
    forged=code.replace(co_consts=(Foreign(),)+code.co_consts[1:])
    assert not protocol._same_code(forged,code) and events==[]


def test_constructor_function_with_forged_builtins_keeps_the_old_callbacks():
    script='''
import symox as s,importlib,types,builtins
m=importlib.import_module('symox.polys.matrices.ddm')
p=importlib.import_module('symox.polys.matrices._lll_native')
l=importlib.import_module('symox.polys.matrices.lll')
register=p.register;p.register=lambda *args:None
m=importlib.reload(m);p.register=register
f=m.DDM.__init__;events=[]
def fake_all(values):events.append('all');return builtins.all(values)
old=m.__dict__['__builtins__'];fake=dict(builtins.__dict__);fake['all']=fake_all
m.__dict__['__builtins__']=fake
try:g=types.FunctionType(f.__code__,m.__dict__,f.__name__,f.__defaults__,f.__closure__)
finally:m.__dict__['__builtins__']=old
m.DDM.__init__=g
p.register(m.DDM,l)
a=m.DDM([[s.ZZ(2),s.ZZ(0)],[s.ZZ(3),s.ZZ(2)]],(2,2),s.ZZ)
assert not p.prefix(a,s.QQ(3,4),True)
events.clear();actual=a.lll_transform();trace=events[:];events.clear()
p.prefix=lambda *args:False
expected=a.lll_transform()
assert actual==expected and trace==events and trace
'''
    subprocess.run([sys.executable,'-c',script],text=True,capture_output=True,check=True,timeout=30)


def test_scratch_finalizers_observe_the_original_complete_matrices(monkeypatch):
    a=matrix([[2,0],[3,2]]);events=[]
    def finalize(self):events.append([[str(c) for c in row] for row in self])
    monkeypatch.setattr(DDM,'__del__',finalize,raising=False)
    actual=tuple(ints(v) for v in a.lll_transform());trace=events[:];events.clear()
    with monkeypatch.context() as patch:
        patch.setattr(protocol,'prefix',lambda *args:False)
        expected=tuple(ints(v) for v in a.lll_transform())
    assert actual==expected and events==trace and trace
