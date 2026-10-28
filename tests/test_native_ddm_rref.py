"""Exact QQ Gauss-Jordan reduction, row objects and rank certificates."""
import copy
import importlib
import random

import pytest
import sympy as sp
import symox as sx
from symox.polys.domains import QQ,ZZ
from symox.polys.matrices.ddm import DDM
from symox.polys.matrices.domainmatrix import DomainMatrix
from sympy.polys.matrices.dense import ddm_irref as reference

dense=importlib.import_module('symox.polys.matrices.dense')


def matrix(rows):return sp.Matrix([[sp.Rational(str(v)) for v in row] for row in rows])
def qq(rows):return [[QQ(int(sp.numer(v)),int(sp.denom(v))) for v in row] for row in rows]


@pytest.mark.parametrize('partial',[False,True])
@pytest.mark.parametrize('seed',range(14))
def test_exact_rref_row_permutations_and_pivots_match_both_oracles(seed,partial):
    rng=random.Random(821901+seed);m=1+seed%5;n=1+(seed*3)%7
    raw=[[sp.Rational(rng.randrange(-5,6),rng.randrange(1,8)) for _ in range(n)] for _ in range(m)]
    if seed%3==0 and m>1:raw[-1]=raw[0][:]
    if seed%4==0:raw[0]=[sp.S.Zero]*n
    a=qq(raw);old=[row[:] for row in a]
    identities={id(row):i for i,row in enumerate(a)}
    old_ids={id(row):i for i,row in enumerate(old)}
    expected_pivots=reference(old,_partial_pivot=partial)
    actual=sx._core._ddm_irref(a,partial)
    assert actual==expected_pivots
    assert a==old
    assert [identities[id(row)] for row in a]==[old_ids[id(row)] for row in old]
    reduced,pivots=sp.Matrix(raw).rref()
    assert matrix(a)==reduced and actual==list(pivots)
    assert all(type(v) is sx._core.DomainElement and v.parent()==QQ for row in a for v in row)


@pytest.mark.parametrize('partial',[False,True])
def test_augmented_transformation_certifies_exact_row_space_and_rank(partial):
    raw=[[0,2,4,6],[1,3,5,7],[2,6,10,14],[0,0,0,0]]
    a=qq(raw);pivots=dense.ddm_irref(a,_partial_pivot=partial)
    m=len(raw);n=len(raw[0]);augmented=qq([row+[int(i==j) for j in range(m)] for i,row in enumerate(raw)])
    augmented_pivots=dense.ddm_irref(augmented,_partial_pivot=partial)
    result=matrix(augmented);transform=result[:,n:]
    assert transform*sp.Matrix(raw)==matrix(a)
    assert transform.det()!=0
    assert [p for p in augmented_pivots if p<n]==pivots==[0,1]
    assert matrix(a).rank()==len(pivots)


def test_partial_pivot_compares_sub_float_rational_differences_and_first_ties():
    tiny=QQ(1,10**350)
    a=[[QQ(1),QQ(0)],[QQ(1)+tiny,QQ(1)]];rows=list(a)
    assert dense.ddm_irref(a,_partial_pivot=True)==[0,1]
    assert a[0] is rows[1] and a[1] is rows[0]
    b=[[QQ(1),QQ(0)],[QQ(-1),QQ(1)]];rows=list(b)
    assert dense.ddm_irref(b,_partial_pivot=True)==[0,1]
    assert b[0] is rows[0] and b[1] is rows[1]


@pytest.mark.parametrize('partial',[False,True])
def test_nearly_dependent_exact_rows_never_lose_rank(partial):
    tiny=QQ(1,10**400)
    a=[[QQ(1),QQ(1)],[QQ(1),QQ(1)+tiny]]
    assert dense.ddm_irref(a,_partial_pivot=partial)==[0,1]
    assert a==[[QQ(1),QQ(0)],[QQ(0),QQ(1)]]


def test_stock_ddm_and_domainmatrix_dense_routes_use_native_reduction(monkeypatch):
    calls=[];original=sx._core._ddm_irref
    def native(*args):
        value=original(*args);calls.append(value);return value
    monkeypatch.setattr(sx._core,'_ddm_irref',native)
    a=DDM(qq([[0,2,4],[1,3,5]]),(2,3),QQ)
    saved=[row[:] for row in a];rows=list(a)
    reduced,pivots=a.rref()
    assert pivots==[0,1] and type(reduced) is DDM and reduced.domain==QQ
    assert list(a)==saved and all(row is old for row,old in zip(a,rows))
    dm=DomainMatrix.from_rep(a)
    result,pivots=dm.rref(method='GJ_dense')
    assert pivots==(0,1) and matrix(result.to_list())==matrix(reduced)
    assert len(calls)==2 and all(value is not None for value in calls)


def test_domainmatrix_inverse_uses_reduction_and_has_exact_two_sided_certificate(monkeypatch):
    calls=[];original=sx._core._ddm_irref
    def native(*args):
        value=original(*args);calls.append(value);return value
    monkeypatch.setattr(sx._core,'_ddm_irref',native)
    raw=[[0,2,1],[1,1,0],[2,0,1]]
    a=DomainMatrix(qq(raw),(3,3),QQ)
    inverse=a.inv();actual=matrix(inverse.to_list());source=sp.Matrix(raw)
    assert source*actual==sp.eye(3) and actual*source==sp.eye(3)
    assert calls==[[0,1,2]]


@pytest.mark.parametrize('a',[[],[[]],[[],[]]])
@pytest.mark.parametrize('partial',[False,True])
def test_empty_shapes_have_no_cell_or_row_mutation(a,partial):
    a=copy.deepcopy(a);rows=list(a)
    assert sx._core._ddm_irref(a,partial)==[]
    assert all(row is old for row,old in zip(a,rows))


@pytest.mark.parametrize('partial',[False,True])
def test_aliased_rows_decline_untouched_then_keep_the_general_feedback(partial):
    row=[QQ(2),QQ(4)];a=[row,row];old=copy.deepcopy(a)
    assert sx._core._ddm_irref(a,partial) is None
    assert a[0] is a[1] is row and a==[[QQ(2),QQ(4)]]*2
    expected=reference(old,_partial_pivot=partial)
    assert dense.ddm_irref(a,_partial_pivot=partial)==expected
    assert a==old and a[0] is a[1] is row


@pytest.mark.parametrize('a',[
    [[QQ(2),QQ(4)],[QQ(1)]],
    [[],[QQ(1)]],
])
def test_nonrectangular_shapes_preserve_values_before_refusal_and_error_order(a):
    old=copy.deepcopy(a);before=copy.deepcopy(a)
    assert sx._core._ddm_irref(a,False) is None and a==before
    try:expected=reference(old)
    except Exception as error:
        with pytest.raises(type(error)) as actual:dense.ddm_irref(a)
        assert str(actual.value)==str(error)
    else:assert dense.ddm_irref(a)==expected
    assert a==old


@pytest.mark.parametrize('a',[
    [[1,2],[3,4]],[[QQ(1),0],[QQ(0),QQ(1)]],
    [[ZZ(2),ZZ(4)],[ZZ(1),ZZ(3)]],[[sx.Float('0.1',60)]],
])
def test_other_scalar_representations_keep_the_old_divisions_and_types(a):
    old=copy.deepcopy(a);before=copy.deepcopy(a)
    assert sx._core._ddm_irref(a,False) is None and a==before
    try:expected=reference(old)
    except Exception as error:
        with pytest.raises(type(error)) as actual:dense.ddm_irref(a)
        assert str(actual.value)==str(error)
    else:assert dense.ddm_irref(a)==expected
    assert a==old
    assert [[type(value) for value in row] for row in a]==[[type(value) for value in row] for row in old]


def test_huge_python_int_keeps_its_original_overflow_instead_of_promoting_to_qq():
    a=[[2**10000]]
    assert sx._core._ddm_irref(a,False) is None
    with pytest.raises(OverflowError):dense.ddm_irref(a)
    assert a==[[2**10000]]


def test_custom_truth_of_partial_pivot_is_never_probed_for_native_eligibility():
    class Pivot:
        calls=0
        def __bool__(self):self.calls+=1;raise KeyboardInterrupt
    partial=Pivot();a=[[QQ(1)]]
    assert sx._core._ddm_irref(a,partial) is None and partial.calls==0
    with pytest.raises(KeyboardInterrupt):dense.ddm_irref(a,_partial_pivot=partial)
    assert partial.calls==1
    partial.calls=0
    assert dense.ddm_irref([],_partial_pivot=partial)==[] and partial.calls==0


def test_custom_rows_decline_without_any_getter_or_iterator_calls():
    class Row(list):
        calls=0
        def __iter__(self):self.calls+=1;raise AssertionError('row iterator')
        def __len__(self):self.calls+=1;raise AssertionError('row length')
    row=Row([QQ(1)])
    assert sx._core._ddm_irref([row],False) is None and row.calls==0
    with pytest.raises(AssertionError,match='row length'):dense.ddm_irref([row])
    assert row.calls==1


@pytest.mark.parametrize('name',['max','range','abs','len','enumerate'])
def test_replaced_builtin_preserves_its_hook(name,monkeypatch):
    def interrupted(*args,**kwargs):raise KeyboardInterrupt
    monkeypatch.setattr(dense,name,interrupted,raising=False)
    with pytest.raises(KeyboardInterrupt):dense.ddm_irref([[QQ(1)]],_partial_pivot=True)


def test_added_inplace_coefficient_hook_is_observable(monkeypatch):
    def interrupted(*args):raise KeyboardInterrupt
    monkeypatch.setattr(sx._core.DomainElement,'__imul__',interrupted,raising=False)
    with pytest.raises(KeyboardInterrupt):dense.ddm_irref([[QQ(2)]])


def test_changed_stock_ddm_row_assignment_is_not_bypassed(monkeypatch):
    a=DDM(qq([[0,2],[1,3]]),(2,2),QQ);rows=list(a)
    def interrupted(*args):raise KeyboardInterrupt
    monkeypatch.setattr(DDM,'__setitem__',interrupted)
    with pytest.raises(KeyboardInterrupt):dense.ddm_irref(a)
    assert all(row is old for row,old in zip(a,rows))


@pytest.mark.parametrize('hook',['__bool__','__pow__'])
def test_coefficient_hook_before_reload_stays_observable(hook):
    import subprocess,sys
    source=r'''
import importlib
import symox as S
from symox.polys.domains import QQ
module=importlib.import_module('symox.polys.matrices.dense')
a=[[QQ(2)]]
def interrupted(*args):raise KeyboardInterrupt
setattr(S._core.DomainElement,HOOK,interrupted)
importlib.reload(module)
try:module.ddm_irref(a)
except KeyboardInterrupt:pass
else:raise AssertionError('coefficient callback bypassed after reload')
assert a[0][0]==QQ(2)
'''.replace('HOOK',repr(hook))
    subprocess.run([sys.executable,'-c',source],check=True,timeout=30)


@pytest.mark.parametrize('partial',[False,True])
def test_native_interruption_preserves_completed_normalization_prefix(partial):
    import signal,subprocess,sys
    if not hasattr(signal,'setitimer'):pytest.skip('POSIX interval timer unavailable')
    source=r'''
import signal
import symox as S
from symox.polys.domains import QQ
n=300;zero=QQ(0);one=QQ(1);two=QQ(2)
a=[[zero]*n for _ in range(n)]
for i in range(n):a[i][i]=two
rows=list(a)
def interrupted(*args):raise KeyboardInterrupt
signal.signal(signal.SIGALRM,interrupted)
try:
    signal.setitimer(signal.ITIMER_REAL,0.001)
    try:S._core._ddm_irref(a,PARTIAL)
    except KeyboardInterrupt:pass
    else:raise AssertionError('native reduction lost interruption')
finally:signal.setitimer(signal.ITIMER_REAL,0)
assert all(row is old for row,old in zip(a,rows))
seen_unmodified=False
for i,row in enumerate(a):
    for j,value in enumerate(row):
        assert type(value) is S._core.DomainElement and value.parent()==QQ
        if i!=j:assert value==zero
        elif value==two:seen_unmodified=True
        else:assert value==one and not seen_unmodified
'''.replace('PARTIAL',repr(partial))
    subprocess.run([sys.executable,'-c',source],check=True,timeout=30)


@pytest.mark.parametrize('kernel,name',[
    ('rref','max'),('rref','range'),('rref','abs'),('rref','len'),('rref','enumerate'),
    ('matmul','sum'),('matmul','zip'),('matmul','list'),('matmul','map'),
])
def test_builtin_hook_before_reload_never_becomes_a_trusted_identity(kernel,name):
    import subprocess,sys
    source=r'''
import importlib
from symox.polys.domains import QQ
module=importlib.import_module('symox.polys.matrices.dense')
a=[[QQ(2)]]
def interrupted(*args,**kwargs):raise KeyboardInterrupt
setattr(module,NAME,interrupted)
importlib.reload(module)
try:
    if KERNEL=='rref':module.ddm_irref(a,_partial_pivot=True)
    else:module.ddm_imatmul(a,[[QQ(2)]],[[QQ(3)]])
except KeyboardInterrupt:pass
else:raise AssertionError('a builtin callback became trusted after reload')
assert a[0][0]==QQ(1 if NAME=='enumerate' else 2)
'''.replace('KERNEL',repr(kernel)).replace('NAME',repr(name))
    subprocess.run([sys.executable,'-c',source],check=True,timeout=30)
