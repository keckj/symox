"""DDM Bareiss values, in-place state, rows, aliases and protocol fallbacks."""
import importlib
import random

import pytest
import sympy as sp
import symox as sx
from sympy.polys.matrices.dense import ddm_idet as reference
from symox.polys.matrices.ddm import DDM
from symox.polys.matrices.domainmatrix import DomainMatrix

dense = importlib.import_module('symox.polys.matrices.dense')


def cells(raw, k):
    return [[k(int(sp.numer(value)), int(sp.denom(value))) if k is sx.QQ else k(int(value))
             for value in row] for row in raw]


def aliases(rows, result):
    classes = {}
    return [classes.setdefault(id(value), len(classes)) for row in rows for value in row] + [classes.setdefault(id(result),len(classes))]


@pytest.mark.parametrize('k', (sx.ZZ, sx.QQ))
@pytest.mark.parametrize('seed', range(24))
def test_exact_values_full_mutation_and_aliases_match_reference(k, seed):
    rng = random.Random(191129+seed)
    n = seed % 7
    raw = [[sp.Rational(rng.randrange(-9,10),rng.randrange(1,8) if k is sx.QQ else 1)
            for _ in range(n)] for _ in range(n)]
    if n > 1 and seed % 3 == 0:
        raw[0][0] = 0
    if n > 1 and seed % 4 == 0:
        for row in raw:
            row[0] = 0
    if n > 2 and seed % 5 == 0:
        raw[-1] = raw[0][:]
    a = cells(raw,k); old = [row[:] for row in a]
    a_rows = {id(row): i for i,row in enumerate(a)}
    old_rows = {id(row): i for i,row in enumerate(old)}
    expected = reference(old,k)
    actual = sx._core._ddm_idet(a,k)
    assert actual is not None and actual == expected
    assert a == old
    assert [a_rows[id(row)] for row in a] == [old_rows[id(row)] for row in old]
    assert aliases(a,actual) == aliases(old,expected)
    assert actual.parent() == k
    assert sp.Rational(str(actual)) == sp.Matrix(raw).det()
    assert all(value.parent() == k for row in a for value in row)


@pytest.mark.parametrize('k', (sx.ZZ, sx.QQ))
@pytest.mark.parametrize('seed', range(5))
def test_big_coefficients_and_multiple_pivot_swaps(k, seed):
    big = 2**(151+seed*17)+3**57
    raw = [[0,0,sp.Rational(big,17) if k is sx.QQ else big,3],
           [0,5,7,11], [13,17,19,23], [29,31,37,41]]
    a = cells(raw,k); old = [row[:] for row in a]
    expected = reference(old,k)
    actual = dense.ddm_idet(a,k)
    assert actual == expected and a == old
    assert sp.Rational(str(actual)) == sp.Matrix(raw).det()


@pytest.mark.parametrize('fmt', ('dense','sparse'))
@pytest.mark.parametrize('k', (sx.ZZ,sx.QQ))
def test_domainmatrix_and_ddm_call_native_without_mutating_input(monkeypatch,fmt,k):
    a = DomainMatrix(cells([[0,2,3],[4,5,6],[7,8,10]],k),(3,3),k)
    if fmt == 'sparse': a = a.to_sparse()
    before = a.copy(); source_ids = {(i,j):id(value) for i,row in a.to_dod().items() for j,value in row.items()}
    calls = []
    old = sx._core._ddm_idet
    def spy(*args):
        result = old(*args); calls.append(result is not None); return result
    monkeypatch.setattr(sx._core,'_ddm_idet',spy)
    assert a.det() == k(-5)
    assert a == before and source_ids == {(i,j):id(value) for i,row in a.to_dod().items() for j,value in row.items()}
    assert calls == [True]


@pytest.mark.parametrize('kind', ('alias','row','outer','integers','mixed','ragged','wide','tall','zero_columns'))
def test_nonstandard_representations_preserve_complete_old_path(monkeypatch,kind):
    events = []
    class Row(list):
        def __getitem__(self,key):
            events.append(('get',key)); return super().__getitem__(key)
        def __setitem__(self,key,value):
            events.append(('set',key)); return super().__setitem__(key,value)
    class Outer(list):
        def __getitem__(self,key):
            events.append(('outer',key)); return super().__getitem__(key)
    def make():
        a = cells([[0,2,3],[4,5,6],[7,8,10]],sx.ZZ)
        if kind == 'alias': a[2] = a[1]
        if kind == 'row': a[1] = Row(a[1])
        if kind == 'outer': a = Outer(a)
        if kind == 'integers': a = [[0,2,3],[4,5,6],[7,8,10]]
        if kind == 'mixed': a[1][0] = 4
        if kind == 'ragged': a[1].pop()
        if kind == 'wide': a = [row+[sx.ZZ.one] for row in a]
        if kind == 'tall': a.append(a[0][:])
        if kind == 'zero_columns': a = [[],[]]
        return a
    def run():
        a = make()
        try:
            result = dense.ddm_idet(a,sx.ZZ)
            answer = ('value', str(result))
        except Exception as error:
            answer = (type(error),str(error))
        return answer, [[str(value) for value in row] for row in a]
    actual = run(); trace = events[:]; events.clear()
    with monkeypatch.context() as patch:
        patch.setattr(sx._core,'_ddm_idet',lambda *args:None)
        expected = run()
    assert actual == expected and events == trace


@pytest.mark.parametrize('hook', ('exquo','one','bool','multiply','subtract','list_setitem'))
def test_custom_operator_callbacks_keep_order_and_types(monkeypatch,hook):
    events = []
    k = sx.QQ
    a0 = cells([[0,2,3],[4,5,6],[7,8,10]],k)
    def run():
        a = [row[:] for row in a0]
        if hook == 'list_setitem': a = DDM(a,(3,3),k)
        result = dense.ddm_idet(a,k)
        return str(result), [[str(value) for value in row] for row in a]
    with monkeypatch.context() as patch:
        if hook in ('exquo','one'):
            cls = type(k); name = hook; old = cls.__dict__[name]
            if hook == 'one':
                def get(self):
                    events.append('one'); return old.__get__(self,cls)
                patch.setattr(cls,name,property(get))
            else:
                def callback(self,*args):
                    events.append('exquo'); return old(self,*args)
                patch.setattr(cls,name,callback)
        elif hook == 'list_setitem':
            old = DDM.__setitem__
            def setitem(self,key,value):
                events.append(('setitem',key)); return old(self,key,value)
            patch.setattr(DDM,'__setitem__',setitem)
        else:
            cls = type(k.one)
            name = {'bool':'__bool__','multiply':'__mul__','subtract':'__sub__'}[hook]
            old = getattr(cls,name)
            def callback(self,*args):
                events.append(hook); return old(self,*args)
            patch.setattr(cls,name,callback)
        actual = run(); trace = events[:]; events.clear()
        with monkeypatch.context() as disable:
            disable.setattr(sx._core,'_ddm_idet',lambda *args:None)
            expected = run()
        assert actual == expected and events == trace and trace


def test_exception_identity_and_published_prefix_match_reference(monkeypatch):
    error = RuntimeError('division callback')
    observed = []
    old = type(sx.ZZ).exquo
    def fail_second(self,a,b):
        observed.append((str(a),str(b)))
        if len(observed) == 2: raise error
        return old(self,a,b)
    monkeypatch.setattr(type(sx.ZZ),'exquo',fail_second)
    a = cells([[1,2,3],[4,5,6],[7,8,10]],sx.ZZ)
    with pytest.raises(RuntimeError) as caught:
        dense.ddm_idet(a,sx.ZZ)
    assert caught.value is error and a[1][1] == sx.ZZ(-3) and a[1][2] == sx.ZZ(6)


@pytest.mark.parametrize('k', (sx.RR,sx.QQ.poly_ring(sx.Symbol('t')),sx.GF(5)))
def test_domains_outside_route_keep_complete_original_algorithm(k):
    a = [[k(1),k(2)],[k(3),k(4)]]
    old = [row[:] for row in a]
    assert dense.ddm_idet(a,k) == reference(old,k)
    assert a == old


def test_nonsquare_public_error_is_unchanged():
    from symox.polys.matrices.exceptions import DMNonSquareMatrixError
    a = DomainMatrix(cells([[1,2,3],[4,5,6]],sx.QQ),(2,3),sx.QQ)
    with pytest.raises(DMNonSquareMatrixError): a.det()
    with pytest.raises(DMNonSquareMatrixError): a.rep.det()
