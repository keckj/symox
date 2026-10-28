"""Exact dense Berkowitz coefficients, Cayley-Hamilton and protocol fallbacks."""
import builtins
import importlib
import random

import pytest
import sympy as sp
import symox as sx
from sympy.polys.matrices.dense import ddm_berk as reference
from symox.polys.matrices.ddm import DDM
from symox.polys.matrices.sdm import SDM
from symox.polys.matrices.domainmatrix import DomainMatrix

dense = importlib.import_module('symox.polys.matrices.dense')


def elements(raw,k):
    return [[k(int(sp.numer(v)),int(sp.denom(v))) if k is sx.QQ else k(int(v)) for v in row] for row in raw]


def alias_partition(rows,inputs):
    classes = {}
    for value in inputs:
        classes.setdefault(id(value),len(classes))
    return [classes.setdefault(id(v),len(classes)) for row in rows for v in row]


@pytest.mark.parametrize('k',(sx.ZZ,sx.QQ))
@pytest.mark.parametrize('seed',range(22))
def test_complete_reference_coefficients_shapes_aliases_and_cayley_hamilton(k,seed):
    rng = random.Random(931721+seed); n = seed%8
    raw = [[sp.Rational(rng.randrange(-9,10),rng.randrange(1,13) if k is sx.QQ else 1)
            for _ in range(n)] for _ in range(n)]
    if n>1 and seed%3==0: raw[-1]=raw[0][:]
    if seed%5==0: raw = [[0]*n for _ in range(n)]
    a = elements(raw,k)
    if n>1 and seed%4==0: a[-1] = a[0]  # Read-only row aliases are valid.
    old = [row[:] for row in a]
    original_rows = [id(row) for row in a]
    original_cells = [[id(v) for v in row] for row in a]
    actual = sx._core._ddm_berk(a,k)
    expected = reference(old,k)
    assert actual is not None and actual == expected
    assert type(actual) is list and len(actual)==n+1
    assert all(type(row) is list and len(row)==1 for row in actual)
    assert len({id(row) for row in actual})==n+1
    source = [v for row in a for v in row]
    assert alias_partition(actual,source)==alias_partition(expected,source)
    assert all(type(row[0]) is type(k.one) and row[0].parent()==k for row in actual)
    assert original_rows == [id(row) for row in a]
    assert original_cells == [[id(v) for v in row] for row in a]
    rational_matrix = sp.Matrix([[sp.Rational(str(v)) for v in row] for row in a])
    coeffs = [sp.Rational(str(row[0])) for row in actual]
    assert coeffs == rational_matrix.charpoly().all_coeffs()
    polynomial = sp.zeros(n)
    for coefficient in coeffs:
        polynomial = polynomial*rational_matrix+coefficient*sp.eye(n)
    assert polynomial == sp.zeros(n)


@pytest.mark.parametrize('k',(sx.ZZ,sx.QQ))
def test_large_coefficients_and_coprime_denominators(k):
    raw = [[sp.Rational(2**211+3,17) if k is sx.QQ else 2**211+3, 2**187-9, 0],
           [sp.Rational(-3**81,5**47) if k is sx.QQ else -3**81, 1, 2],
           [5,sp.Rational(7**61,11**31) if k is sx.QQ else 7**61,3]]
    a=elements(raw,k)
    actual=dense.ddm_berk(a,k)
    assert actual==reference([row[:] for row in a],k)
    assert [sp.Rational(str(row[0])) for row in actual]==sp.Matrix(raw).charpoly().all_coeffs()


@pytest.mark.parametrize('entry',('DDM','charpoly_berk','charpoly'))
@pytest.mark.parametrize('k',(sx.ZZ,sx.QQ))
def test_stock_dense_public_paths_use_native_and_preserve_input(monkeypatch,entry,k):
    a=DomainMatrix(elements([[1,2,3],[4,5,6],[7,8,10]],k),(3,3),k)
    before=a.copy(); ids=[[id(v) for v in row] for row in a.rep]
    calls=[]; original=sx._core._ddm_berk
    def spy(*args):
        result=original(*args);calls.append(result is not None);return result
    monkeypatch.setattr(sx._core,'_ddm_berk',spy)
    actual=a.rep.charpoly() if entry=='DDM' else getattr(a,entry)()
    assert actual==[row[0] for row in reference([row[:] for row in a.rep],k)]
    assert calls==[True] and a==before and ids==[[id(v) for v in row] for row in a.rep]


@pytest.mark.parametrize('hook',('one','zero','neg','mul','add','sum','python_matmul','native_matmul','recursion','ddm_getitem'))
def test_skipped_protocols_keep_exact_old_callback_trace(monkeypatch,hook):
    a=elements([[1,2,3],[4,5,6],[7,8,10]],sx.QQ)
    original=dense.ddm_berk; events=[]
    with monkeypatch.context() as patch:
        if hook in ('one','zero'):
            cls=type(sx.QQ); descriptor=cls.__dict__[hook]
            def get(self):
                events.append(hook);return descriptor.__get__(self,cls)
            patch.setattr(cls,hook,property(get))
        elif hook in ('neg','mul','add'):
            cls=type(sx.QQ.one);name={'neg':'__neg__','mul':'__mul__','add':'__add__'}[hook];old=getattr(cls,name)
            def callback(self,*args):
                events.append(hook);return old(self,*args)
            patch.setattr(cls,name,callback)
        elif hook=='sum':
            def callback(*args):
                events.append('sum');return builtins.sum(*args)
            patch.setattr(dense,'sum',callback,raising=False)
        elif hook=='python_matmul':
            old=dense.ddm_imatmul
            def callback(*args):
                events.append('matmul');return old(*args)
            patch.setattr(dense,'ddm_imatmul',callback)
        elif hook=='native_matmul':
            old=sx._core._ddm_imatmul
            def callback(*args):
                events.append('native_matmul');return old(*args)
            patch.setattr(sx._core,'_ddm_imatmul',callback)
        elif hook=='recursion':
            def callback(rows,k):
                events.append(('recursive',len(rows)));return original(rows,k)
            patch.setattr(dense,'ddm_berk',callback)
        else:
            a=DDM(a,(3,3),sx.QQ);old=DDM.__getitem__
            def get(self,key):
                events.append(('get',key));return old(self,key)
            patch.setattr(DDM,'__getitem__',get)
        actual=original(a,sx.QQ);trace=events[:];events.clear()
        with monkeypatch.context() as disable:
            disable.setattr(sx._core,'_ddm_berk',lambda *args:None)
            expected=original(a,sx.QQ)
        assert actual==expected and events==trace and trace


@pytest.mark.parametrize('storage',('outer','row','tuple','integers','mixed','domain_mismatch','wide','tall','ragged'))
def test_noncanonical_storage_retains_complete_old_result_or_error(monkeypatch,storage):
    events=[]
    class Custom(list):
        def __getitem__(self,key):
            events.append(key);return super().__getitem__(key)
    def make():
        a=elements([[1,2],[3,4]],sx.ZZ)
        if storage=='outer':a=Custom(a)
        elif storage=='row':a[1]=Custom(a[1])
        elif storage=='tuple':a=tuple(tuple(row) for row in a)
        elif storage=='integers':a=[[1,2],[3,4]]
        elif storage=='mixed':a[1][0]=3
        elif storage=='domain_mismatch':a[1][0]=sx.QQ(1,2)
        elif storage=='wide':a[0].append(sx.ZZ.one);a[1].append(sx.ZZ.one)
        elif storage=='tall':a.append(a[0][:])
        elif storage=='ragged':a[1].pop()
        return a
    def run():
        try:
            result=dense.ddm_berk(make(),sx.ZZ)
            return ('value',[[(str(v),str(v.parent()) if hasattr(v,'parent') else type(v).__name__) for v in row] for row in result])
        except Exception as error:return (type(error),str(error))
    actual=run();trace=events[:];events.clear()
    with monkeypatch.context() as patch:
        patch.setattr(sx._core,'_ddm_berk',lambda *args:None)
        expected=run()
    assert actual==expected and events==trace


def test_custom_callback_exception_identity(monkeypatch):
    a=elements([[1,2],[3,4]],sx.QQ);error=RuntimeError('berkowitz callback')
    def fail(*args):raise error
    monkeypatch.setattr(dense,'ddm_imatmul',fail)
    with pytest.raises(RuntimeError) as caught:dense.ddm_berk(a,sx.QQ)
    assert caught.value is error


@pytest.mark.parametrize('domain',('RR','CC'))
@pytest.mark.parametrize('precision',(1,53,139,140,211))
def test_inexact_domains_keep_previous_grouping_bits(monkeypatch,domain,precision):
    k=(sx.RealField if domain=='RR' else sx.ComplexField)(precision)
    a=[[k(sx.Rational(i*5+j+1,7) + (sx.I*sx.Rational(i+j+1,11) if domain=='CC' else 0)) for j in range(3)] for i in range(3)]
    def payload(result):
        return [[(str(value.parent()),sx.srepr(k.to_sympy(value)),
                  tuple(sorted((repr(atom._mpf_),atom._prec) for atom in k.to_sympy(value).atoms(sx.Float))))
                 for value in row] for row in result]
    actual=payload(dense.ddm_berk(a,k))
    with monkeypatch.context() as patch:
        patch.setattr(sx._core,'_ddm_berk',lambda *args:None)
        assert actual==payload(dense.ddm_berk(a,k))


@pytest.mark.parametrize('k',(sx.GF(5),sx.GF(4),sx.QQ.poly_ring(sx.Symbol('t')),sx.EX))
def test_other_domains_retain_complete_legacy_algorithm(k):
    a=[[k(1),k(2)],[k(3),k(4)]]
    assert dense.ddm_berk(a,k)==reference(a,k)


def test_sparse_berkowitz_stays_sparse_without_new_dense_conversion(monkeypatch):
    a=SDM({0:{0:sx.QQ(2)},1:{1:sx.QQ(3)},2:{2:sx.QQ(5)}},(3,3),sx.QQ)
    def forbidden(*args):raise AssertionError('new dense conversion')
    monkeypatch.setattr(SDM,'to_ddm',forbidden)
    monkeypatch.setattr(sx._core,'_ddm_berk',forbidden)
    assert a.charpoly()==[sx.QQ(1),sx.QQ(-10),sx.QQ(31),sx.QQ(-30)]


@pytest.mark.parametrize('poison', ('shape','domain','rows','cols'))
def test_poisoned_ddm_metadata_retains_old_low_level_path(monkeypatch,poison):
    a=DDM(elements([[1,2],[3,4]],sx.QQ),(2,2),sx.QQ)
    setattr(a,poison, sx.ZZ if poison=='domain' else (5,5) if poison=='shape' else 5)
    calls=[];original=sx._core._ddm_berk
    def spy(*args):
        result=original(*args);calls.append(result is not None);return result
    monkeypatch.setattr(sx._core,'_ddm_berk',spy)
    actual=dense.ddm_berk(a,sx.QQ)
    assert calls[0] is False
    with monkeypatch.context() as patch:
        patch.setattr(sx._core,'_ddm_berk',lambda *args:None)
        expected=dense.ddm_berk(a,sx.QQ)
    assert actual==expected


def test_ddm_metadata_getter_is_not_observed_by_new_eligibility(monkeypatch):
    a=DDM(elements([[1,2],[3,4]],sx.QQ),(2,2),sx.QQ)
    events=[]
    def get_shape(self):
        events.append('shape');return (2,2)
    monkeypatch.setattr(DDM,'shape',property(get_shape),raising=False)
    actual=dense.ddm_berk(a,sx.QQ)
    assert events==[]
    with monkeypatch.context() as patch:
        patch.setattr(sx._core,'_ddm_berk',lambda *args:None)
        expected=dense.ddm_berk(a,sx.QQ)
    assert actual==expected and events==[]
