"""Exact FF dense/sparse domains, denominator identity and public presentation."""
import importlib
import random

import pytest
import sympy as sp
import symox as sx
from symox.polys.domains import ZZ,QQ
from symox.polys.matrices.ddm import DDM
from symox.polys.matrices.sdm import SDM
from symox.polys.matrices.domainmatrix import DomainMatrix
from sympy.polys.matrices.dense import ddm_irref_den as dense_reference
from sympy.polys.matrices.sdm import sdm_rref_den as sparse_reference

dense=importlib.import_module('symox.polys.matrices.dense')
sparse=importlib.import_module('symox.polys.matrices.sdm')


def ground(raw,K):
    return [[K(int(sp.numer(v)),int(sp.denom(v))) if K is QQ else K(int(v)) for v in row] for row in raw]


def matrix(rows):return sp.Matrix([[sp.Rational(str(v)) for v in row] for row in rows])
def dict_rows(rows):return {i:{j:v for j,v in enumerate(row) if v} for i,row in enumerate(rows) if any(row)}
def sparse_matrix(rows,m,n):return sp.Matrix(m,n,lambda i,j:sp.Rational(str(rows.get(i,{}).get(j,0))))


def aliases(rows,den):
    classes={}
    return [classes.setdefault(id(v),len(classes)) for row in rows for v in row]+[classes.setdefault(id(den),len(classes))]


@pytest.mark.parametrize('K',[ZZ,QQ])
@pytest.mark.parametrize('seed',range(22))
def test_dense_denominator_cells_types_aliases_and_rank_match_reference(K,seed):
    rng=random.Random(77641+seed);m=1+seed%5;n=1+(seed*3)%8
    raw=[[sp.Rational(rng.randrange(-7,8),rng.randrange(1,7) if K is QQ else 1) for j in range(n)] for i in range(m)]
    if seed%4==0:raw[0]=[0]*n
    if seed%3==0 and m>1:raw[-1]=raw[0][:]
    a=ground(raw,K);old=[row[:] for row in a]
    ids={id(row):i for i,row in enumerate(a)};old_ids={id(row):i for i,row in enumerate(old)}
    expected_den,expected_pivots=dense_reference(old,K)
    actual=sx._core._ddm_irref_den(a,K)
    assert actual is not None
    den,pivots=actual
    assert (den,pivots)==(expected_den,expected_pivots) and a==old
    assert [ids[id(row)] for row in a]==[old_ids[id(row)] for row in old]
    assert aliases(a,den)==aliases(old,expected_den)
    reference,p=sp.Matrix(raw).rref()
    assert matrix(a)/sp.Rational(str(den))==reference and pivots==list(p)
    assert den.parent()==K and all(v.parent()==K for row in a for v in row)
    if pivots:assert den is a[0][pivots[0]]


@pytest.mark.parametrize('K',[ZZ,QQ])
@pytest.mark.parametrize('seed',range(24))
def test_sparse_denominator_values_key_order_and_no_input_mutation(K,seed):
    rng=random.Random(188223+seed);m=2+seed%5;n=2+(seed*5)%9
    raw=[[sp.Rational(rng.randrange(-5,6),rng.randrange(1,8) if K is QQ else 1) for j in range(n)] for i in range(m)]
    if seed%4==0:raw[-1]=raw[0][:]
    if seed%5==0:raw[0]=[0]*n
    a=dict_rows(ground(raw,K))
    # Exercise noncanonical insertion order rather than only dense conversions.
    order=list(a);rng.shuffle(order)
    a={i:dict(rng.sample(list(a[i].items()),len(a[i]))) for i in order}
    before={i:row.copy() for i,row in a.items()};ids={i:{j:id(v) for j,v in row.items()} for i,row in a.items()}
    expected,expected_den,expected_pivots=sparse_reference(a,K)
    native=sx._core._sdm_rref_den(a,K)
    if len(a)==1:assert native is None;actual=sparse.sdm_rref_den(a,K)
    else:assert native is not None;actual=native
    out,den,pivots=actual
    assert (out,den,pivots)==(expected,expected_den,expected_pivots)
    assert repr(out)==repr(expected) and list(out)==list(expected)
    assert [list(row) for row in out.values()]==[list(row) for row in expected.values()]
    assert a==before and {i:{j:id(v) for j,v in row.items()} for i,row in a.items()}==ids
    reduced,p=sp.Matrix(raw).rref()
    assert sparse_matrix(out,m,n)/sp.Rational(str(den))==reduced and pivots==list(p)
    assert den.parent()==K
    assert all(v.parent()==K for row in out.values() for v in row.values())
    assert all(out[i][pivot] is den for i,pivot in enumerate(pivots))


@pytest.mark.parametrize('K',[ZZ,QQ])
@pytest.mark.parametrize('kind',['dense','sparse'])
def test_augmented_transformation_and_denominator_have_exact_certificate(K,kind):
    raw=[[0,2,4,6],[1,3,5,7],[2,6,10,14],[0,0,0,0]];m=len(raw);n=len(raw[0])
    augmented=[row+[int(i==j) for j in range(m)] for i,row in enumerate(raw)]
    if kind=='dense':
        a=ground(raw,K);den,pivots=dense.ddm_irref_den(a,K)
        aug=ground(augmented,K);ad,ap=dense.ddm_irref_den(aug,K)
        reduced=matrix(a)/sp.Rational(str(den));full=matrix(aug)/sp.Rational(str(ad))
    else:
        a,den,pivots=sparse.sdm_rref_den(dict_rows(ground(raw,K)),K)
        aug,ad,ap=sparse.sdm_rref_den(dict_rows(ground(augmented,K)),K)
        reduced=sparse_matrix(a,m,n)/sp.Rational(str(den));full=sparse_matrix(aug,m,n+m)/sp.Rational(str(ad))
    transform=full[:,n:]
    assert transform.det()!=0 and transform*sp.Matrix(raw)==reduced
    assert [p for p in ap if p<n]==pivots==[0,1]
    assert reduced.rank()==len(pivots)


@pytest.mark.parametrize('K',[ZZ,QQ])
def test_sparse_large_column_indices_stay_sparse_and_match_presentation(K):
    c=2**50
    a={9:{c:K(3),1:K(2)},3:{1:K(4),c:K(5)}}
    expected=sparse_reference(a,K);actual=sx._core._sdm_rref_den(a,K)
    assert actual==expected and repr(actual[0])==repr(expected[0])
    assert actual[2]==[1,c]
    assert all(len(row)==1 for row in actual[0].values())


@pytest.mark.parametrize('K',[ZZ,QQ])
def test_sparse_one_row_preserves_the_shallow_copy_and_existing_pivot_object(K):
    a={7:{9:K(2),1:K(3),4:K(5)}}
    assert sx._core._sdm_rref_den(a,K) is None
    out,den,pivots=sparse.sdm_rref_den(a,K)
    assert pivots==[1] and den is a[7][1]
    assert out[0] is not a[7] and list(out[0])==[9,1,4]
    assert all(out[0][j] is value for j,value in a[7].items())


@pytest.mark.parametrize('K',[ZZ,QQ])
def test_dense_shared_previous_pivots_and_returned_denominator(K):
    a=ground([[2,0,0],[0,3,0],[0,0,4]],K)
    den,pivots=dense.ddm_irref_den(a,K)
    assert pivots==[0,1,2] and den==K(24)
    assert den is a[0][0] is a[1][1]
    assert den is not a[2][2]


@pytest.mark.parametrize('K',[ZZ,QQ])
def test_auto_ff_cd_uses_native_sparse_and_public_presentation_is_identical(K,monkeypatch):
    raw=[[2,1,3],[4,3,7],[1,2,4]]
    rows=ground(raw,K)
    a=DomainMatrix(rows,(3,3),K)
    calls=[];original=sx._core._sdm_rref_den
    def native(*args):
        value=original(*args);calls.append(value);return value
    monkeypatch.setattr(sx._core,'_sdm_rref_den',native)
    actual,den,pivots=a.rref_den()
    assert calls and all(value is not None for value in calls)
    monkeypatch.setattr(sx._core,'_sdm_rref_den',lambda *args:None)
    expected,ed,ep=a.rref_den()
    assert den==ed and pivots==ep and actual.to_dok()==expected.to_dok()
    assert repr(actual)==repr(expected) and str(actual)==str(expected)
    assert actual.to_list()==expected.to_list()


@pytest.mark.parametrize('K',[ZZ,QQ])
def test_explicit_dense_route_and_ddm_preserve_input(K,monkeypatch):
    a=DDM(ground([[2,1,3],[4,3,7],[1,2,4]],K),(3,3),K);before=[row[:] for row in a]
    calls=[];original=sx._core._ddm_irref_den
    def native(*args):
        value=original(*args);calls.append(value);return value
    monkeypatch.setattr(sx._core,'_ddm_irref_den',native)
    result,den,pivots=a.rref_den()
    assert calls and all(value is not None for value in calls)
    assert list(a)==before and result.domain==K
    dm=DomainMatrix.from_rep(a)
    value,d,p=dm.rref_den(method='FF_dense')
    assert value.to_list()==list(result) and d==den and p==tuple(pivots)


@pytest.mark.parametrize('K',[ZZ,QQ])
@pytest.mark.parametrize('raw',[[],[[]],[[],[]],[[0,0],[0,0]]])
def test_empty_and_zero_dense_matrices_return_typed_unit(K,raw):
    a=ground(raw,K);expected=ground(raw,K)
    den,p=dense.ddm_irref_den(a,K)
    ed,ep=dense_reference(expected,K)
    assert (a,den,p)==(expected,ed,ep) and den.parent()==K
    out,den,p=sparse.sdm_rref_den({},K)
    assert out=={} and p==[] and den==K.one and den.parent()==K


@pytest.mark.parametrize('K',[ZZ,QQ])
def test_raw_ints_and_aliases_decline_without_changing_the_old_contract(K):
    a=[[2,0],[0,3]];old=[r[:] for r in a]
    assert sx._core._ddm_irref_den(a,K) is None and a==old
    actual=dense.ddm_irref_den(a,K);expected=dense_reference(old,K)
    assert (a,actual)==(old,expected)
    assert [[type(v) for v in row] for row in a]==[[type(v) for v in row] for row in old]
    assert type(actual[0]) is type(expected[0])
    row=[K(2),K(3)];a=[row,row];oldrow=row[:];old=[oldrow,oldrow]
    assert sx._core._ddm_irref_den(a,K) is None
    assert dense.ddm_irref_den(a,K)==dense_reference(old,K) and a==old
    assert a[0] is a[1]
    row={0:K(2),1:K(3)};a={0:row,1:row}
    assert sx._core._sdm_rref_den(a,K) is None
    assert sparse.sdm_rref_den(a,K)==sparse_reference(a,K) and a[0] is a[1]


@pytest.mark.parametrize('K',[ZZ,QQ])
def test_invalid_sparse_storage_declines_and_preserves_error_behavior(K):
    for a in ({0:{0:K(0)},1:{1:K(1)}},{0:{},1:{1:K(1)}},
              {-1:{0:K(1)},2:{1:K(1)}},{2**80:{0:K(1)},1:{1:K(1)}}):
        assert sx._core._sdm_rref_den(a,K) is None
        try:expected=sparse_reference(a,K)
        except Exception as error:
            with pytest.raises(type(error),match='^'+__import__('re').escape(str(error))+'$'):
                sparse.sdm_rref_den(a,K)
        else:assert sparse.sdm_rref_den(a,K)==expected


@pytest.mark.parametrize('kind',['dense','sparse'])
def test_custom_domain_callbacks_and_interruptions_keep_the_original_prefix(kind):
    class Domain:
        def __init__(self):self.calls=[]
        def __getattr__(self,name):
            self.calls.append(name)
            if name=='exquo':
                def interrupted(a,b):raise KeyboardInterrupt('custom exquo')
                return interrupted
            return getattr(ZZ,name)
    K=Domain();a=ground([[2,0],[0,3]],ZZ)
    if kind=='sparse':a=dict_rows(a)
    helper=sx._core._ddm_irref_den if kind=='dense' else sx._core._sdm_rref_den
    assert helper(a,K) is None and K.calls==[]
    fn=dense.ddm_irref_den if kind=='dense' else sparse.sdm_rref_den
    ref=dense_reference if kind=='dense' else sparse_reference
    before=ground([[2,0],[0,3]],ZZ)
    if kind=='sparse':before=dict_rows(before)
    oldK=Domain()
    with pytest.raises(KeyboardInterrupt,match='custom exquo'):ref(before,oldK)
    with pytest.raises(KeyboardInterrupt,match='custom exquo'):fn(a,K)
    assert K.calls==oldK.calls and a==before


@pytest.mark.parametrize('wrapped',[False,True])
def test_real_exquo_override_is_observed_even_with_copied_function_metadata(monkeypatch,wrapped):
    import functools
    T=type(ZZ);original=T.exquo;calls=[]
    def exquo(self,a,b):calls.append((a,b));raise KeyboardInterrupt('override')
    if wrapped:exquo=functools.wraps(original)(exquo)
    monkeypatch.setattr(T,'exquo',exquo)
    a=ground([[2,0],[0,3]],ZZ)
    with pytest.raises(KeyboardInterrupt,match='override'):dense.ddm_irref_den(a,ZZ)
    assert calls and a==ground([[2,0],[0,6]],ZZ)
    calls.clear();a={0:{0:ZZ(2)},1:{1:ZZ(3)}}
    with pytest.raises(KeyboardInterrupt,match='override'):sparse.sdm_rref_den(a,ZZ)
    assert calls and a=={0:{0:ZZ(2)},1:{1:ZZ(3)}}


def test_stock_exquo_delegate_override_remains_observable(monkeypatch):
    module=importlib.import_module('symox.polys.domains._algebraic_methods');calls=[]
    def callback(*args):calls.append(args);raise KeyboardInterrupt('delegate')
    monkeypatch.setattr(module,'_native_exquo',callback)
    with pytest.raises(KeyboardInterrupt,match='delegate'):
        dense.ddm_irref_den(ground([[2,0],[0,3]],ZZ),ZZ)
    assert calls


@pytest.mark.parametrize('flag',['is_Exact','is_Composite'])
def test_exactness_properties_are_not_evaluated_by_eligibility(monkeypatch,flag):
    T=type(ZZ);calls=[]
    def flag_value(self):calls.append(self);raise KeyboardInterrupt(flag)
    monkeypatch.setattr(T,flag,property(flag_value))
    guard=importlib.import_module('symox.polys.matrices._fraction_free')
    assert guard.stock_domain(ZZ) is False and not calls
    if flag=='is_Exact':
        with pytest.raises(KeyboardInterrupt,match=flag):
            sparse.sdm_rref_den({0:{0:ZZ(2)},1:{1:ZZ(3)}},ZZ)
        assert len(calls)==1
    # Empty sparse storage never needs the exactness/Composite property.
    assert sparse.sdm_rref_den({},ZZ)==({},ZZ.one,[])


@pytest.mark.parametrize('method',['items','values'])
def test_sdm_instance_methods_and_class_hooks_are_not_bypassed(monkeypatch,method):
    a=SDM({0:{0:ZZ(2)},1:{1:ZZ(3)}},(2,2),ZZ);calls=[]
    def callback(*args):calls.append(method);raise KeyboardInterrupt(method)
    a.__dict__[method]=callback
    assert sparse._ff_stock_sdm(a) is None and calls==[]
    # The general branch uses items, while the one-row branch uses values.
    if method=='values':a.pop(1)
    with pytest.raises(KeyboardInterrupt,match=method):sparse.sdm_rref_den(a,ZZ)
    assert calls==[method]
    a=SDM({0:{0:ZZ(2)},1:{1:ZZ(3)}},(2,2),ZZ);calls.clear()
    monkeypatch.setattr(SDM,'items',property(lambda self:callback()),raising=False)
    assert sparse._ff_stock_sdm(a) is None and not calls
    with pytest.raises(KeyboardInterrupt,match=method):sparse.sdm_rref_den(a,ZZ)
    assert calls==[method]


def test_sdm_instance_metadata_keys_are_not_compared_by_eligibility():
    calls=[]
    class Key:
        def __hash__(self):return hash('items')
        def __eq__(self,other):calls.append(other);raise AssertionError('metadata equality')
    a=SDM({0:{0:ZZ(2)},1:{1:ZZ(3)}},(2,2),ZZ)
    a.__dict__[Key()]=1
    assert sparse._ff_stock_sdm(a) is None and calls==[]


def test_row_and_key_subclasses_are_refused_without_protocol_calls():
    calls=[]
    class Row(dict):
        def items(self):calls.append('items');raise AssertionError
    class Key(int):
        def __int__(self):calls.append('int');raise AssertionError
    a={0:Row({0:ZZ(1)}),1:{1:ZZ(2)}}
    assert sx._core._sdm_rref_den(a,ZZ) is None
    a={Key(0):{0:ZZ(1)},1:{1:ZZ(2)}}
    assert sx._core._sdm_rref_den(a,ZZ) is None and calls==[]


@pytest.mark.parametrize('kind,name',[('dense','len'),('dense','range'),('dense','enumerate'),
    ('sparse','len'),('sparse','min'),('sparse','sorted'),('sparse','zip'),('sparse','set'),
    ('sparse','list'),('sparse','dict'),('sparse','enumerate')])
def test_builtin_overrides_keep_the_complete_fallback(monkeypatch,kind,name):
    import builtins
    module=dense if kind=='dense' else sparse;calls=[];original=getattr(builtins,name)
    def callback(*args,**kwargs):calls.append(name);return original(*args,**kwargs)
    monkeypatch.setattr(module,name,callback,raising=False)
    # Reach each builtin naturally: set/list/dict appear after cancellation.
    a=ground([[2,1,3],[4,3,7],[1,2,4]],ZZ)
    if kind=='sparse':a=dict_rows(a)
    helper='_ddm_irref_den' if kind=='dense' else '_sdm_rref_den'
    def forbidden(*args):raise AssertionError('native route bypassed an override')
    monkeypatch.setattr(sx._core,helper,forbidden)
    fn=dense.ddm_irref_den if kind=='dense' else sparse.sdm_rref_den
    fn(a,ZZ)
    # Some builtins belong to only a branch of the algorithm; refusal itself
    # is the contract, and observed callbacks remain the old algorithm's choice.


def test_native_fraction_free_has_no_runtime_sympy_dependency():
    import os,subprocess,sys
    script="""import sys
sys.modules['sympy']=None
from symox.polys.domains import ZZ,QQ
from symox.polys.matrices import DomainMatrix
for K in (ZZ,QQ):
 a=DomainMatrix([[K(2),K(1)],[K(4),K(3)]],(2,2),K)
 for method in ('auto','FF_dense'):
  r,d,p=a.rref_den(method=method)
  assert p==(0,1) and r.to_list()==[[d,K.zero],[K.zero,d]]
"""
    subprocess.run([sys.executable,'-c',script],env=os.environ.copy(),check=True)


@pytest.mark.parametrize('name',['RR','GF5','GF6'])
def test_other_domains_keep_the_original_dense_and_sparse_routes(name):
    from symox.polys.domains import RR,GF
    K={'RR':RR,'GF5':GF(5),'GF6':GF(6)}[name]
    def outcome(fn,a):
        try:return ('result',fn(a,K),a)
        except Exception as error:return ('error',type(error),str(error),a)
    a=[[K(2),K(1)],[K(3),K(4)]];old=[row[:] for row in a]
    assert sx._core._ddm_irref_den(a,K) is None
    assert outcome(dense.ddm_irref_den,a)==outcome(dense_reference,old)
    a=dict_rows([[K(2),K(1)],[K(3),K(4)]]);old={i:row.copy() for i,row in a.items()}
    assert sx._core._sdm_rref_den(a,K) is None
    assert outcome(sparse.sdm_rref_den,a)==outcome(sparse_reference,old)


def test_added_inplace_multiplication_hook_is_not_bypassed(monkeypatch):
    calls=[]
    def imul(self,other):calls.append((self,other));raise KeyboardInterrupt('imul')
    monkeypatch.setattr(type(ZZ(1)),'__imul__',imul,raising=False)
    a=ground([[1,2,4],[0,0,3]],ZZ);old=[r[:] for r in a]
    with pytest.raises(KeyboardInterrupt,match='imul'):dense_reference(old,ZZ)
    calls.clear()
    with pytest.raises(KeyboardInterrupt,match='imul'):dense.ddm_irref_den(a,ZZ)
    assert calls and a==old


def test_native_signal_interruptions_preserve_dense_row_identity_and_sparse_input():
    import os,subprocess,sys
    script="""import signal
import symox as sx
from symox.polys.domains import ZZ,QQ
from symox.polys.matrices.dense import ddm_irref_den
from symox.polys.matrices.sdm import sdm_rref_den
signal.signal(signal.SIGALRM,lambda *args:(_ for _ in ()).throw(KeyboardInterrupt()))
for K in (ZZ,QQ):
 for kind in ('dense','sparse'):
  n=350
  if kind=='dense':
   a=[[K(int(i==j)) for j in range(n)] for i in range(n)]
   rows=[id(row) for row in a]
  else:
   a={i:{j:K((i+j)%5+1) for j in range(n)} for i in range(n)}
   before={i:{j:id(v) for j,v in row.items()} for i,row in a.items()}
  helper='_ddm_irref_den' if kind=='dense' else '_sdm_rref_den'
  original=getattr(sx._core,helper);entered=[]
  def native(*args):entered.append(True);return original(*args)
  setattr(sx._core,helper,native)
  signal.setitimer(signal.ITIMER_REAL,0.03)
  try:
   (ddm_irref_den if kind=='dense' else sdm_rref_den)(a,K)
  except KeyboardInterrupt:pass
  else:raise AssertionError('signal was swallowed')
  finally:signal.setitimer(signal.ITIMER_REAL,0)
  assert entered
  setattr(sx._core,helper,original)
  if kind=='dense':
   assert [id(row) for row in a]==rows
   assert all(v==K(int(i==j)) and v.parent()==K for i,row in enumerate(a) for j,v in enumerate(row))
  else:
   assert {i:{j:id(v) for j,v in row.items()} for i,row in a.items()}==before
"""
    subprocess.run([sys.executable,'-c',script],env=os.environ.copy(),check=True,timeout=30)
