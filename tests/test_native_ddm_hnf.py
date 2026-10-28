"""Ordinary ZZ HNF, exact column lattices, publication and full fallbacks."""
import ast
import copy
import importlib
import random
import types
from fractions import Fraction
from pathlib import Path

import pytest
import sympy as sp
import symox as sx
from sympy.matrices.normalforms import hermite_normal_form as reference
from symox.polys.matrices.ddm import DDM
from symox.polys.matrices.domainmatrix import DomainMatrix

nf=importlib.import_module('symox.polys.matrices.normalforms')
public=importlib.import_module('symox.matrices.normalforms')
printer=importlib.import_module('symox.polys.domains._element_print')


def domain_matrix(rows,shape):
    return DomainMatrix.from_rep(DDM([[sx.ZZ(v) for v in row] for row in rows],shape,sx.ZZ))


def integers(rows):return [[int(v) for v in row] for row in rows]


def old_function(namespace=None):
    module=ast.parse(Path(nf.__file__).read_text())
    node=next(n for n in module.body if isinstance(n,ast.FunctionDef) and n.name=='_hermite_normal_form')
    node.body=[n for n in node.body if not (isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='native' for t in n.targets))
               and not (isinstance(n,ast.If) and isinstance(n.test,ast.Compare) and isinstance(n.test.left,ast.Name) and n.test.left.id=='native')]
    code=compile(ast.fix_missing_locations(ast.Module(body=[node],type_ignores=[])),nf.__file__,'exec')
    return types.FunctionType(next(c for c in code.co_consts if type(c)is types.CodeType),nf.__dict__ if namespace is None else namespace)


def multiply(a,b,n):
    return [[sum(a[i][k]*b[k][j] for k in range(n)) for j in range(n)] for i in range(len(a))]


def determinant(a):
    n=len(a);a=[list(map(Fraction,row)) for row in a];result=Fraction(1)
    for j in range(n):
        p=next((i for i in range(j,n) if a[i][j]),None)
        if p is None:return 0
        if p!=j:a[p],a[j]=a[j],a[p];result=-result
        pivot=a[j][j];result*=pivot
        for i in range(j+1,n):
            q=a[i][j]/pivot
            for k in range(j+1,n):a[i][k]-=q*a[j][k]
    return result


def lattice_certificate(source,full,k,operations,n):
    transform=[[int(i==j) for j in range(n)] for i in range(n)]
    for i,j,a,b,c,d in operations:
        element=[[int(p==q) for q in range(n)] for p in range(n)]
        if i==j:
            assert (a,b,c,d)==(-1,0,-1,0);element[i][i]=-1
        else:
            assert a*d-b*c==1
            element[i][i],element[j][i],element[i][j],element[j][j]=a,b,c,d
        transform=multiply(transform,element,n)
    assert abs(determinant(transform))==1
    assert multiply(source,transform,n)==full
    assert all(value==0 for row in full for value in row[:k])
    h=[row[k:] for row in full];pivots=[]
    for j in range(n-k):
        pivot=max(i for i,row in enumerate(h) if row[j])
        assert h[pivot][j]>0
        assert all(0<=h[pivot][right]<h[pivot][j] for right in range(j+1,n-k))
        pivots.append(pivot)
    assert pivots==sorted(set(pivots))


@pytest.mark.parametrize('shape',[(0,0),(0,5),(1,0),(4,0),(1,1),(1,7),(4,2),(2,4),(4,4),(4,6),(6,4),(6,6)])
@pytest.mark.parametrize('seed',range(12))
def test_exact_canonical_hnf_and_independent_unimodular_lattice_certificate(shape,seed,monkeypatch):
    m,n=shape;rng=random.Random(98317+seed*31+m*17+n)
    rows=[[rng.randrange(-15,16) for _ in range(n)] for _ in range(m)]
    if seed%5==0:rows=[[0]*n for _ in range(m)]
    if seed%4==0 and n>1:
        for row in rows:row[-1]=row[0]*3
    a=domain_matrix(rows,shape);before=a.to_list();calls=[];native=sx._core._ddm_hnf
    def record(*args,**kwargs):
        result=native(*args,**kwargs);calls.append(result);return result
    monkeypatch.setattr(sx._core,'_ddm_hnf',record)
    actual=nf.hermite_normal_form(a)
    expected=reference(sp.Matrix(m,n,[v for row in rows for v in row]))
    assert actual.shape==expected.shape and integers(actual.to_list())==expected.tolist()
    assert a.to_list()==before and calls and all(value is not None for value in calls)
    work=a.to_ddm().copy();k,operations=native(work,shape,DDM,True)
    full=integers(work);lattice_certificate(rows,full,k,operations,n)
    assert actual.shape==(m,n-k) and [row[k:] for row in full]==expected.tolist()


@pytest.mark.parametrize('entry',['matrix','immutable','dense','sparse'])
def test_real_public_routes_select_native_without_changing_return_kind(entry,monkeypatch):
    rows=[[12,6,4],[3,9,6],[2,16,14]];calls=[];native=sx._core._ddm_hnf
    def record(*args):
        result=native(*args);calls.append(result);return result
    monkeypatch.setattr(sx._core,'_ddm_hnf',record)
    if entry in ('matrix','immutable'):
        source=(sx.Matrix if entry=='matrix' else sx.ImmutableMatrix)(rows)
        result=public.hermite_normal_form(source)
        assert result==sx.Matrix([[10,0,2],[0,15,3],[0,0,2]]) and source==sx.Matrix(rows)
    else:
        source=domain_matrix(rows,(3,3))
        if entry=='sparse':source=source.to_sparse()
        result=nf.hermite_normal_form(source)
        assert integers(result.to_list())==[[10,0,2],[0,15,3],[0,0,2]] and result.domain==sx.ZZ
        assert integers(source.to_list())==rows
    assert calls==[0]


@pytest.mark.parametrize('digits',[100,500,1500])
def test_large_integer_and_negative_floor_remainders_preserve_exactness(digits):
    big=10**digits
    rows=[[big+1,-big-3,7],[big-1,big+9,-11],[3,-5,13]]
    a=domain_matrix(rows,(3,3));out=nf.hermite_normal_form(a)
    assert integers(out.to_list())==reference(sp.Matrix(rows)).tolist()


def test_scratch_row_and_cell_publication_matches_original_identity_masks():
    rows=[[12,6,4],[3,9,6],[2,16,14]]
    source=domain_matrix(rows,(3,3));a=source.to_ddm().copy();b=source.to_ddm().copy()
    original_a=[row[:] for row in a];original_b=[row[:] for row in b];a_rows=list(a);b_rows=list(b)
    # Run the preserved scalar body with copy returning its actual supplied
    # workspace, so publication and untouched cells are directly observable.
    original_copy=DDM.copy
    DDM.copy=lambda value:value
    try:
        expected=old_function()(DomainMatrix.from_rep(b))
        actual=nf._hermite_normal_form(DomainMatrix.from_rep(a))
    finally:DDM.copy=original_copy
    assert actual==expected and integers(a)==integers(b)
    assert all(left is right for left,right in zip(a,a_rows)) and all(left is right for left,right in zip(b,b_rows))
    assert [[value is original_a[i][j] for j,value in enumerate(row)] for i,row in enumerate(a)]==[[value is original_b[i][j] for j,value in enumerate(row)] for i,row in enumerate(b)]


@pytest.mark.parametrize('name',['_gcdex','add_columns','range','len'])
def test_actual_algorithm_helpers_and_builtins_keep_callbacks(name,monkeypatch):
    events=[];original=nf.__dict__.get(name,getattr(__import__('builtins'),name,None))
    def hook(*args):events.append(name);return original(*args)
    monkeypatch.setattr(nf,name,hook,raising=False)
    a=domain_matrix([[-12,6,4],[3,9,6],[2,16,14]],(3,3))
    expected=old_function()(a);before=events[:];events.clear();actual=nf._hermite_normal_form(a);after=events[:]
    assert actual==expected and after==before and after


@pytest.mark.parametrize('name',['__add__','__mul__','__rmul__','__neg__','__lt__','__floordiv__','__mod__','__eq__','__ne__','parent'])
def test_consumed_coefficient_protocols_are_not_bypassed(name,monkeypatch):
    events=[];original=getattr(sx._core.DomainElement,name)
    def hook(*args):events.append(name);return original(*args)
    monkeypatch.setattr(sx._core.DomainElement,name,hook)
    a=domain_matrix([[-12,6,4],[3,9,6],[2,16,14]],(3,3))
    expected=old_function()(a);before=events[:];events.clear();actual=nf._hermite_normal_form(a);after=events[:]
    assert integers(actual.to_list())==integers(expected.to_list()) and after==before and after


def test_zz_gcdex_callback_and_gcdex_function_global_are_preserved(monkeypatch):
    events=[];original=sx._core.Domain.gcdex
    def gcd(domain,*args):events.append('gcdex');return original(domain,*args)
    monkeypatch.setattr(sx._core.Domain,'gcdex',gcd)
    a=domain_matrix([[12,6,4],[3,9,6],[2,16,14]],(3,3))
    expected=old_function()(a);before=events[:];events.clear();actual=nf._hermite_normal_form(a);after=events[:]
    assert actual==expected and after==before and after


def test_clone_namespace_does_not_hide_helper_own_range_override(monkeypatch):
    namespace=dict(nf.__dict__);current=types.FunctionType(nf._hermite_normal_form.__code__,namespace)
    old=old_function(namespace);events=[];real=range
    def hook(*args):events.append(args);return real(*args)
    monkeypatch.setattr(nf,'range',hook,raising=False)
    a=domain_matrix([[2,3],[4,5]],(2,2))
    expected=old(a);before=events[:];events.clear();actual=current(a);after=events[:]
    assert actual==expected and after==before and after


@pytest.mark.parametrize('name',['_gcdex','add_columns'])
def test_missing_or_nonfunction_late_helper_is_not_touched_for_zero_lattice(name):
    class Hook:
        def __getattribute__(self,attribute):raise AssertionError('new helper attribute call')
        def __call__(self,*args):raise AssertionError('new helper call')
    namespace=dict(nf.__dict__);namespace[name]=Hook()
    current=types.FunctionType(nf._hermite_normal_form.__code__,namespace)
    a=domain_matrix([[0,0],[0,0]],(2,2))
    assert current(a).shape==(2,0)
    namespace.pop(name)
    assert current(a).shape==(2,0)


@pytest.mark.parametrize('option',['valid_D','invalid_D','check_rank_fallback','custom_check'])
def test_outer_modular_and_rank_option_behavior_stays_on_its_original_path(option,monkeypatch):
    calls=[];native=sx._core._ddm_hnf
    def record(*args):
        result=native(*args);calls.append(result);return result
    monkeypatch.setattr(sx._core,'_ddm_hnf',record)
    if option=='check_rank_fallback':
        a=domain_matrix([[1,2],[2,4]],(2,2));out=nf.hermite_normal_form(a,D=sx.ZZ(3),check_rank=True)
        assert integers(out.to_list())==reference(sp.Matrix([[1,2],[2,4]])).tolist() and calls
    elif option=='invalid_D':
        with pytest.raises(Exception) as error:nf.hermite_normal_form(domain_matrix([[1,2],[3,4]],(2,2)),D=sx.ZZ(0))
        assert type(error.value).__name__=='DMDomainError' and calls==[]
    elif option=='custom_check':
        class Check:
            def __bool__(self):raise AssertionError('unused check_rank')
        nf.hermite_normal_form(domain_matrix([[1,2],[3,4]],(2,2)),check_rank=Check())
        assert calls
    else:
        out=nf.hermite_normal_form(domain_matrix([[1,2],[3,4]],(2,2)),D=sx.ZZ(2))
        assert integers(out.to_list())==reference(sp.Matrix([[1,2],[3,4]]),D=2).tolist() and calls==[]


def test_other_domains_and_raw_python_integer_representations_keep_full_paths():
    a=DomainMatrix.from_list([[1,2],[3,4]],sx.QQ)
    with pytest.raises(Exception) as error:nf.hermite_normal_form(a)
    assert type(error.value).__name__=='DMDomainError'
    raw=DomainMatrix.from_rep(DDM([[1,2],[3,4]],(2,2),sx.ZZ))
    assert nf.hermite_normal_form(raw)==old_function()(raw)


def test_unconsumed_bool_and_presentation_hooks_do_not_disable_math(monkeypatch):
    calls=[];native=sx._core._ddm_hnf
    def record(*args):
        result=native(*args);calls.append(result);return result
    def forbidden(*args):raise AssertionError('unconsumed hook')
    monkeypatch.setattr(sx._core,'_ddm_hnf',record)
    for name in ('__bool__','__repr__','__str__'):
        monkeypatch.setattr(sx._core.DomainElement,name,forbidden)
    out=nf.hermite_normal_form(domain_matrix([[2,3],[4,5]],(2,2)))
    assert out.shape==(2,2) and calls==[0]


def test_print_equality_wrapper_global_type_callback_is_preserved(monkeypatch):
    events=[];real=type
    def hook(value):events.append('type');return real(value)
    monkeypatch.setattr(printer,'type',hook,raising=False)
    a=domain_matrix([[2,3],[4,5]],(2,2))
    expected=old_function()(a);before=events[:];events.clear();actual=nf._hermite_normal_form(a);after=events[:]
    assert integers(actual.to_list())==integers(expected.to_list()) and after==before and after


def test_actual_working_copy_prefix_and_custom_column_publication_are_not_replayed(monkeypatch):
    events=[];original=DDM.copy;old_columns=nf.add_columns
    def copied(self):events.append('copy');return original(self)
    def columns(*args):events.append('columns');return old_columns(*args)
    monkeypatch.setattr(DDM,'copy',copied);monkeypatch.setattr(nf,'add_columns',columns)
    a=domain_matrix([[2,3],[4,5]],(2,2))
    expected=old_function()(a);before=events[:];events.clear();actual=nf._hermite_normal_form(a);after=events[:]
    assert actual==expected and after==before and after.count('copy')==before.count('copy')


def test_floordiv_native_delegate_cell_change_keeps_actual_callback(monkeypatch):
    method=sx._core.DomainElement.__floordiv__;cell=method.__closure__[0];original=cell.cell_contents;events=[]
    def hook(*args):events.append('divide');return original(*args)
    cell.cell_contents=hook
    try:
        a=domain_matrix([[2,3],[4,5]],(2,2));expected=old_function()(a);before=events[:];events.clear()
        actual=nf._hermite_normal_form(a);after=events[:]
        assert actual==expected and after==before and after
    finally:cell.cell_contents=original


def test_aliases_in_a_custom_copy_fall_back_to_the_original_feedback(monkeypatch):
    monkeypatch.setattr(DDM,'copy',lambda value:value)
    def matrix():
        row=[sx.ZZ(2),sx.ZZ(3)];rep=DDM([row,row],(2,2),sx.ZZ);rep[1]=rep[0]
        return DomainMatrix.from_rep(rep)
    a=matrix();b=matrix();expected=old_function()(a);actual=nf._hermite_normal_form(b)
    assert actual==expected and b.rep[0] is b.rep[1]


def test_numberfields_and_modular_routes_are_not_relabelled_as_native(monkeypatch):
    numberfields=importlib.import_module('symox.polys.numberfields._normalforms')
    def forbidden(*args,**kwargs):raise AssertionError('unported route')
    monkeypatch.setattr(sx._core,'_ddm_hnf',forbidden)
    a=domain_matrix([[1,2],[3,4]],(2,2))
    assert integers(numberfields.hermite_normal_form(a).to_list())==reference(sp.Matrix([[1,2],[3,4]])).tolist()
    assert integers(nf.hermite_normal_form(a,D=sx.ZZ(2)).to_list())==reference(sp.Matrix([[1,2],[3,4]]),D=2).tolist()
