"""In-place exact dense products, scalar representations and alias protocols."""
import copy
from fractions import Fraction
import importlib
import random

import pytest
import sympy as sp
import symox as sx
from symox.polys.domains import ZZ,QQ
from symox.polys.matrices.ddm import DDM
from sympy.polys.matrices.dense import ddm_imatmul as reference

dense=importlib.import_module('symox.polys.matrices.dense')


def rational_matrix(rows):return sp.Matrix([[sp.Rational(str(value)) for value in row] for row in rows])


def converted(rows,domain):
    if domain is int:return [[int(value) for value in row] for row in rows]
    return [[domain(value.numerator,value.denominator) if isinstance(value,Fraction) else domain(value) for value in row] for row in rows]


@pytest.mark.parametrize('domain',[int,ZZ,QQ])
@pytest.mark.parametrize('shape',[(1,1,1),(2,3,4),(5,4,2)])
def test_exact_accumulation_matches_independent_product_and_keeps_row_identity(domain,shape):
    m,k,n=shape;rng=random.Random(1957+m*100+k*10+n)
    def values(rows,cols):
        return [[Fraction(rng.randrange(-9,10),rng.randrange(1,8)) if domain is QQ else rng.randrange(-9,10) for j in range(cols)] for i in range(rows)]
    a=converted(values(m,n),domain);b=converted(values(m,k),domain);c=converted(values(k,n),domain)
    expected=rational_matrix(a)+rational_matrix(b)*rational_matrix(c)
    rows=list(a);left=[row[:] for row in b];right=[row[:] for row in c]
    input_ids=[[id(value) for value in row] for row in b+c]
    assert sx._core._ddm_imatmul(a,b,c) is True
    assert all(row is saved for row,saved in zip(a,rows))
    assert rational_matrix(a)==expected
    assert b==left and c==right
    assert [[id(value) for value in row] for row in b+c]==input_ids
    if domain is int:assert all(type(value) is int for row in a for value in row)
    else:assert all(type(value) is sx._core.DomainElement and value.parent()==domain for row in a for value in row)


@pytest.mark.parametrize('domain',[ZZ,QQ])
def test_stock_ddm_matmul_uses_native_accumulation(domain,monkeypatch):
    a=DDM([[1,2,3],[4,5,6]],(2,3),domain)
    b=DDM([[1,2],[3,4],[5,6]],(3,2),domain)
    observed=[];original=sx._core._ddm_imatmul
    def native(*args):
        result=original(*args);observed.append((args[3],result));return result
    monkeypatch.setattr(sx._core,'_ddm_imatmul',native)
    product=a.matmul(b)
    assert product==DDM([[22,28],[49,64]],(2,2),domain)
    assert type(product) is DDM and product.domain==domain
    assert observed==[(DDM,True)]


def test_arbitrary_size_integers_and_rational_denominators_are_exact():
    n=2**600
    a=[[n,-n]];b=[[n,3]];c=[[n,5],[-7,n]]
    expected=rational_matrix(a)+rational_matrix(b)*rational_matrix(c)
    dense.ddm_imatmul(a,b,c)
    assert rational_matrix(a)==expected and all(type(v) is int for v in a[0])
    a=[[QQ(1,2),QQ(0)]];b=[[QQ(n,37)]];c=[[QQ(-29,31),QQ(701,41)]]
    expected=rational_matrix(a)+rational_matrix(b)*rational_matrix(c)
    dense.ddm_imatmul(a,b,c)
    assert rational_matrix(a)==expected
    assert all(v.parent()==QQ for v in a[0])


@pytest.mark.parametrize('case',['same-left','same-right','shared-left-row','repeated-output-row'])
def test_output_aliases_decline_untouched_and_keep_reference_feedback(case):
    a=[[1,2],[3,4]];b=[[2,1],[0,3]];c=[[1,2],[3,1]]
    if case=='same-left':b=a
    elif case=='same-right':c=a
    elif case=='shared-left-row':b[1]=a[0]
    else:a[1]=a[0]
    expected=copy.deepcopy((a,b,c));rows=list(a);before=[row[:] for row in a]
    assert sx._core._ddm_imatmul(a,b,c) is False
    assert a==before and all(row is old for row,old in zip(a,rows))
    reference(*expected)
    assert dense.ddm_imatmul(a,b,c) is None
    assert a==expected[0] and all(row is old for row,old in zip(a,rows))


def test_shared_read_only_rows_and_immutable_entries_are_safe():
    entry=QQ(3,7);left_row=[entry,entry];right_row=[entry,entry]
    b=[left_row,left_row];c=[right_row,right_row];a=[[QQ(1),QQ(2)],[QQ(3),QQ(4)]]
    expected=rational_matrix(a)+rational_matrix(b)*rational_matrix(c)
    assert sx._core._ddm_imatmul(a,b,c) is True
    assert rational_matrix(a)==expected
    assert all(value is entry for row in b+c for value in row)


@pytest.mark.parametrize('a,b,c',[
    ([],[],[]),([[]],[[]],[]),([[5]],[[]],[]),([[]],[[1]],[[]]),
    ([[2,3]],[[5]],[[7,11],[13,17]]),
    ([[2,3]],[[5,7]],[[7,11],[13]]),
    ([[2,3],[4,5]],[[5]],[[7,11]]),
])
def test_empty_and_nonstandard_shapes_keep_existing_zip_semantics(a,b,c):
    a,b,c=copy.deepcopy((a,b,c));expected=copy.deepcopy((a,b,c))
    before=copy.deepcopy(a)
    assert sx._core._ddm_imatmul(a,b,c) is False
    assert a==before
    reference(*expected);dense.ddm_imatmul(a,b,c)
    assert a==expected[0]


def test_shape_exception_keeps_reference_partial_mutation_order():
    a=[[1],[]];rows=list(a);b=[[1],[1]];c=[[2]]
    assert sx._core._ddm_imatmul(a,b,c) is False and a==[[1],[]]
    with pytest.raises(IndexError):dense.ddm_imatmul(a,b,c)
    assert a==[[3],[]] and all(row is old for row,old in zip(a,rows))


class Untouched:
    def __init__(self):self.calls=[]
    def __iter__(self):self.calls.append('iter');raise AssertionError('iterator consumed')
    def __mul__(self,other):self.calls.append('mul');raise KeyboardInterrupt
    def __rmul__(self,other):self.calls.append('rmul');raise KeyboardInterrupt


@pytest.mark.parametrize('position',range(3))
def test_custom_outer_protocols_and_generators_decline_without_consumption(position):
    args=[[[0]],[[2]],[[3]]];probe=Untouched();args[position]=probe
    assert sx._core._ddm_imatmul(*args) is False
    assert probe.calls==[]
    generator=(row for row in [[2]])
    args[position]=generator
    assert sx._core._ddm_imatmul(*args) is False
    assert next(generator)==[2]


def test_custom_rows_and_coefficients_preserve_hooks_and_interruptions():
    class Row(list):
        def __iter__(self):raise AssertionError('row iterator')
        def __getitem__(self,key):raise AssertionError('row getter')
        def __setitem__(self,key,value):raise AssertionError('row setter')
    assert sx._core._ddm_imatmul([Row([0])],[[2]],[[3]]) is False
    assert sx._core._ddm_imatmul([[0]],[Row([2])],[[3]]) is False
    probe=Untouched();a=[[0]]
    assert sx._core._ddm_imatmul(a,[[probe]],[[3]]) is False
    assert probe.calls==[] and a==[[0]]
    with pytest.raises(KeyboardInterrupt):dense.ddm_imatmul(a,[[probe]],[[3]])
    assert probe.calls==['mul'] and a==[[0]]


@pytest.mark.parametrize('value',[sx.Float('0.1',60),Fraction(1,3),sx.Rational(1,3),sx.GF(7)(2),sx.RR(0.5)])
def test_other_coefficients_and_domains_keep_the_general_arithmetic(value):
    a=[[value]];b=[[value]];c=[[value]]
    assert sx._core._ddm_imatmul(a,b,c) is False
    expected=[[value]];reference(expected,b,c)
    dense.ddm_imatmul(a,b,c)
    assert a==expected and type(a[0][0]) is type(expected[0][0])
    assert getattr(a[0][0],'_prec',None)==getattr(expected[0][0],'_prec',None)


def test_mixed_coefficient_representations_are_not_normalized():
    a=[[0]];b=[[ZZ(2)]];c=[[QQ(3)]]
    assert sx._core._ddm_imatmul(a,b,c) is False
    expected=[[0]];reference(expected,b,c)
    dense.ddm_imatmul(a,b,c)
    assert a==expected and type(a[0][0]) is type(expected[0][0])
    assert a[0][0].parent()==expected[0][0].parent()


def test_replaced_operator_keeps_the_original_callback(monkeypatch):
    def interrupted(*args):raise KeyboardInterrupt
    monkeypatch.setattr(dense,'mul',interrupted)
    a=[[0]]
    with pytest.raises(KeyboardInterrupt):dense.ddm_imatmul(a,[[2]],[[3]])
    assert a==[[0]]


def test_changed_stock_ddm_iteration_uses_the_original_protocol(monkeypatch):
    calls=[];original=DDM.__iter__
    def iterator(value):calls.append(value);return original(value)
    monkeypatch.setattr(DDM,'__iter__',iterator)
    a=DDM([[0]],(1,1),ZZ);b=DDM([[2]],(1,1),ZZ);c=DDM([[3]],(1,1),ZZ)
    calls.clear()
    dense.ddm_imatmul(a,b,c)
    assert a[0][0]==ZZ(6)
    assert len(calls)==3


def test_list_subclass_outer_protocol_is_not_consumed_during_eligibility():
    class Rows(list):
        calls=0
        def __iter__(self):self.calls+=1;raise AssertionError('outer subclass iteration')
    rows=Rows([[2]])
    assert sx._core._ddm_imatmul([[0]],rows,[[3]]) is False
    assert rows.calls==0


@pytest.mark.parametrize('name,expected_calls',[('sum',1),('map',1),('zip',2),('list',1)])
def test_replaced_builtin_preserves_its_original_call_order(name,expected_calls,monkeypatch):
    import builtins
    calls=[];original=getattr(builtins,name)
    def callback(*args):calls.append(name);return original(*args)
    monkeypatch.setattr(dense,name,callback,raising=False)
    a=[[1]]
    dense.ddm_imatmul(a,[[2]],[[3]])
    assert a==[[7]] and calls==[name]*expected_calls


def test_rhs_generator_precedes_coefficient_interruption_exactly_once():
    events=[]
    class Coefficient:
        def __mul__(self,other):events.append('mul');raise KeyboardInterrupt
    def right():events.append('right');yield [3]
    a=[[1]]
    with pytest.raises(KeyboardInterrupt):dense.ddm_imatmul(a,[[Coefficient()]],right())
    assert events==['right','mul'] and a==[[1]]


def test_coefficient_callback_before_dense_reload_is_not_treated_as_stock():
    import subprocess,sys
    subprocess.run([sys.executable,'-c',r'''
import importlib
import symox as S
from symox.polys.domains import QQ
module=importlib.import_module('symox.polys.matrices.dense')
a=[[QQ(0)]];b=[[QQ(2)]];c=[[QQ(3)]]
def interrupted(self,other):raise KeyboardInterrupt
S._core.DomainElement.__mul__=interrupted
importlib.reload(module)
try:module.ddm_imatmul(a,b,c)
except KeyboardInterrupt:pass
else:raise AssertionError('a coefficient callback was bypassed after reload')
assert a[0][0]==QQ(0)
'''],check=True,timeout=30)


def test_native_interruption_keeps_completed_entries_and_row_identity():
    import signal,subprocess,sys
    if not hasattr(signal,'setitimer'):pytest.skip('POSIX interval timer unavailable')
    subprocess.run([sys.executable,'-c',r'''
import signal
import symox as S
n=180
a=[[0]*n for _ in range(n)]
b=[[1]*n for _ in range(n)]
c=[[1]*n for _ in range(n)]
rows=list(a)
def interrupted(*args):raise KeyboardInterrupt
signal.signal(signal.SIGALRM,interrupted)
try:
    signal.setitimer(signal.ITIMER_REAL,0.001)
    try:S._core._ddm_imatmul(a,b,c)
    except KeyboardInterrupt:pass
    else:raise AssertionError('native computation did not propagate the interruption')
finally:signal.setitimer(signal.ITIMER_REAL,0)
assert all(row is old for row,old in zip(a,rows))
seen_unmodified=False
for row in a:
    for value in row:
        assert value in (0,n)
        if value==0:seen_unmodified=True
        else:assert not seen_unmodified
'''],check=True,timeout=30)
