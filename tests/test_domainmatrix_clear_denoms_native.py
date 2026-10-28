"""Exact QQ preparation, representation contracts and complete protocol fallbacks."""
import pytest
import sympy as sp
import symox as sx
from symox.polys.matrices.domainmatrix import DomainMatrix
from symox.polys.matrices.ddm import DDM
from symox.polys.matrices.sdm import SDM
from symox.polys.matrices import _fraction_free as adapter
from sympy.polys.matrices.domainmatrix import DomainMatrix as ReferenceMatrix


METHODS = ('clear_denoms_rowwise', 'clear_denoms')
DATA = (
    ([], (0, 3)),
    ([[], [], []], (3, 0)),
    ([[0, 0], [0, 0]], (2, 2)),
    ([[1, -3], [0, 7]], (2, 2)),
    ([[(1, 6), (-5, 14), 0], [0, 0, 0], [(7, 15), (1, 35), (-2, 9)]], (3, 3)),
    ([[(2**251+7, 2**127-1), (-2**191, 3**81)], [(3**99, 5**61), 0]], (2, 2)),
)


def matrix(rows, shape, fmt, module=sx):
    k = module.QQ
    values = [[k(*v) if isinstance(v, tuple) else k(v) for v in row] for row in rows]
    if module is sx:
        result = DomainMatrix(values, shape, k)
    else:
        result = ReferenceMatrix(values, shape, k)
    return result.to_sparse() if fmt == 'sparse' else result.to_dense()


def snapshot(result, method):
    d, b = result
    numerator = tuple(tuple(str(v) for v in row) for row in b.to_list())
    den = tuple(tuple(str(v) for v in row) for row in d.to_list()) if method.endswith('rowwise') else str(d.element)
    return den, numerator, str(b.domain), b.shape, b.rep.fmt, type(b.rep).__name__


@pytest.mark.parametrize('rows,shape', DATA)
@pytest.mark.parametrize('fmt', ('dense', 'sparse'))
@pytest.mark.parametrize('convert', (False, True))
@pytest.mark.parametrize('method', METHODS)
def test_exact_denominators_reference_and_certificate(rows, shape, fmt, convert, method):
    a = matrix(rows, shape, fmt)
    reference = matrix(rows, shape, fmt, sp)
    actual = getattr(a, method)(convert=convert)
    expected = getattr(reference, method)(convert=convert)
    assert snapshot(actual, method) == snapshot(expected, method)
    d, b = actual
    assert b.domain == (sx.ZZ if convert else sx.QQ)
    assert d.domain == sx.ZZ
    if method.endswith('rowwise'):
        assert d.rep.fmt == 'sparse'
        assert d.convert_to(sx.QQ) * a == b.convert_to(sx.QQ)
    else:
        assert d.convert_to(sx.QQ) * a == b.convert_to(sx.QQ)
    for row in b.to_list():
        for value in row:
            assert type(value) is type(sx.QQ.one)
            assert value.parent() == b.domain


@pytest.mark.parametrize('fmt', ('dense', 'sparse'))
@pytest.mark.parametrize('method', METHODS)
def test_native_route_does_not_call_python_denominator_math(monkeypatch, fmt, method):
    a = matrix(DATA[4][0], DATA[4][1], fmt)
    calls = []
    original = sx._core._dm_clear_denoms
    def spy(*args):
        result = original(*args)
        calls.append(result is not None)
        return result
    monkeypatch.setattr(sx._core, '_dm_clear_denoms', spy)
    getattr(a, method)()
    assert calls == [True]


@pytest.mark.parametrize('method', METHODS)
def test_sparse_order_and_existing_element_aliases(monkeypatch, method):
    shared = sx.QQ(13)
    a = DomainMatrix.from_rep(SDM({3: {2: shared, 0: shared}, 1: {1: sx.QQ(-7)}}, (5, 4), sx.QQ))
    actual = getattr(a, method)()
    with monkeypatch.context() as patch:
        patch.setattr(adapter, 'native_clear_denoms', lambda *args: None)
        old = getattr(a, method)()
    assert snapshot(actual, method) == snapshot(old, method)
    assert list(actual[1].rep) == [3, 1]
    assert list(actual[1].rep[3]) == [2, 0]
    assert actual[1].rep[3][2] is shared
    assert actual[1].rep[3][0] is shared


@pytest.mark.parametrize('convert', (False, True))
@pytest.mark.parametrize('method', METHODS)
def test_noncanonical_sparse_zero_and_empty_rows_keep_old_path(monkeypatch, convert, method):
    for raw in ({1: {2: sx.QQ.zero, 0: sx.QQ(1, 2)}}, {1: {}}):
        a = DomainMatrix.from_rep(SDM(raw, (3, 3), sx.QQ))
        def run():
            try:
                return 'value', snapshot(getattr(a, method)(convert=convert), method)
            except Exception as error:
                return type(error), str(error)
        actual = run()
        with monkeypatch.context() as patch:
            patch.setattr(adapter, 'native_clear_denoms', lambda *args: None)
            assert run() == actual


@pytest.mark.parametrize('method', METHODS)
@pytest.mark.parametrize('protocol', ('matrix', 'representation', 'domain', 'convert', 'truth', 'subclass'))
def test_callbacks_keep_old_order_values_and_exceptions(monkeypatch, method, protocol):
    a = matrix(DATA[4][0], DATA[4][1], 'dense')
    events = []
    with monkeypatch.context() as patch:
        if protocol == 'matrix':
            name = 'to_dod' if method.endswith('rowwise') else 'to_flat_nz'
            old = getattr(DomainMatrix, name)
            def callback(self):
                events.append(('matrix', name)); return old(self)
            patch.setattr(DomainMatrix, name, callback)
        elif protocol == 'representation':
            name = 'to_dod' if method.endswith('rowwise') else 'to_flat_nz'
            old = getattr(DDM, name)
            def callback(self):
                events.append(('rep', name)); return old(self)
            patch.setattr(DDM, name, callback)
        elif protocol == 'domain':
            old = type(sx.QQ).lcm
            def callback(self, x, y):
                events.append(('lcm', str(x), str(y))); return old(self, x, y)
            patch.setattr(type(sx.QQ), 'lcm', callback)
        elif protocol == 'convert':
            convert = False
            old = type(sx.QQ).convert
            def callback(self, *args):
                events.append(('convert', str(self))); return old(self, *args)
            patch.setattr(type(sx.QQ), 'convert', callback)
        elif protocol == 'truth':
            class Truth:
                def __bool__(self):
                    events.append(('truth',)); return True
            convert = Truth()
        else:
            class Custom(DomainMatrix):
                def to_dod(self):
                    events.append(('subclass',)); return super().to_dod()
                def to_flat_nz(self):
                    events.append(('subclass',)); return super().to_flat_nz()
            a = Custom.from_rep(a.rep)
        convert = locals().get('convert', True)
        actual = getattr(a, method)(convert=convert)
        trace = events[:]; events.clear()
        with monkeypatch.context() as disable:
            disable.setattr(adapter, 'native_clear_denoms', lambda *args: None)
            expected = getattr(a, method)(convert=convert)
        assert events == trace and trace
        assert snapshot(actual, method) == snapshot(expected, method)


@pytest.mark.parametrize('method', METHODS)
def test_callback_exception_identity(monkeypatch, method):
    a = matrix(DATA[4][0], DATA[4][1], 'sparse')
    error = RuntimeError('denominator callback')
    def fail(*args):
        raise error
    monkeypatch.setattr(type(sx.QQ), 'denom', fail)
    with pytest.raises(RuntimeError) as caught:
        getattr(a, method)()
    assert caught.value is error


@pytest.mark.parametrize('method', METHODS)
def test_other_domains_retain_complete_fallback(method):
    for k in (sx.ZZ, sx.RR, sx.QQ.poly_ring(sx.Symbol('t'))):
        a = DomainMatrix([[k(1), k(2)]], (1, 2), k)
        d, b = getattr(a, method)(convert=True)
        assert b.to_Matrix() == a.to_Matrix()


@pytest.mark.parametrize('method', METHODS)
def test_sparse_builtin_key_objects_survive(method):
    row, column = int('1025'), int('1027')
    coefficient = sx.QQ(7)
    a = DomainMatrix.from_rep(SDM({row: {column: coefficient}}, (1030, 1030), sx.QQ))
    _, b = getattr(a, method)()
    assert next(iter(b.rep)) is row
    assert next(iter(b.rep[row])) is column
    assert b.rep[row][column] is coefficient


@pytest.mark.parametrize('method', METHODS)
def test_metadata_descriptor_is_not_probed_by_eligibility(monkeypatch, method):
    a = matrix(DATA[4][0], DATA[4][1], 'dense')
    events = []
    def get_rep(self):
        events.append('rep'); return self.__dict__['rep']
    def set_rep(self, value):
        self.__dict__['rep'] = value
    monkeypatch.setattr(DomainMatrix, 'rep', property(get_rep, set_rep), raising=False)
    actual = getattr(a, method)()
    trace = events[:]; events.clear()
    with monkeypatch.context() as patch:
        patch.setattr(adapter, 'native_clear_denoms', lambda *args: None)
        expected = getattr(a, method)()
    assert events == trace
    assert snapshot(actual, method) == snapshot(expected, method)


@pytest.mark.parametrize('method', METHODS)
def test_coefficient_truth_callback_and_poisoned_rep_keep_fallback(monkeypatch, method):
    a = matrix([[1, 0], [0, (1, 3)]], (2, 2), 'dense')
    events = []
    old = type(sx.QQ.one).__bool__
    def truth(value):
        events.append(str(value)); return old(value)
    monkeypatch.setattr(type(sx.QQ.one), '__bool__', truth)
    actual = getattr(a, method)()
    trace = events[:]; events.clear()
    with monkeypatch.context() as patch:
        patch.setattr(adapter, 'native_clear_denoms', lambda *args: None)
        expected = getattr(a, method)()
    assert events == trace and trace
    assert snapshot(actual, method) == snapshot(expected, method)


@pytest.mark.parametrize('fmt', ('dense', 'sparse'))
@pytest.mark.parametrize('convert', (False, True))
@pytest.mark.parametrize('method', METHODS)
def test_full_output_alias_partition_matches_baseline(monkeypatch, fmt, convert, method):
    shared = sx.QQ(5)
    a = DomainMatrix([[shared, shared, sx.QQ.zero], [sx.QQ.zero, sx.QQ.zero, sx.QQ.zero],
                      [sx.QQ(1, 6), sx.QQ(1, 3), sx.QQ.zero]], (3, 3), sx.QQ)
    if fmt == 'sparse':
        a = a.to_sparse()
    def aliases(result):
        d, b = result
        values = list(v for row in b.to_list() for v in row)
        values += list(v for row in d.to_list() for v in row) if method.endswith('rowwise') else [d.element]
        classes = {id(shared): 0}
        return [classes.setdefault(id(value), len(classes)) for value in values]
    actual = getattr(a, method)(convert=convert)
    with monkeypatch.context() as patch:
        patch.setattr(adapter, 'native_clear_denoms', lambda *args: None)
        expected = getattr(a, method)(convert=convert)
    assert aliases(actual) == aliases(expected)


def test_global_empty_shape_beyond_native_index_capacity_keeps_complete_fallback():
    shape = (2**90, 2**91)
    a = DomainMatrix.from_rep(SDM({}, shape, sx.QQ))
    d, b = a.clear_denoms(convert=True)
    assert d.element == sx.ZZ.one and d.domain == sx.ZZ
    assert b.shape == shape and b.domain == sx.ZZ and b.rep == {}


@pytest.mark.parametrize('fmt', ('dense', 'sparse'))
@pytest.mark.parametrize('entry', ('rref', 'rref_den'))
@pytest.mark.parametrize('algorithm', ('auto', 'CD'))
def test_auto_and_cd_fraction_free_use_native_preparation(monkeypatch, fmt, entry, algorithm):
    rows = [[(7 if i == j else 1+i+j, 6) for j in range(6)] for i in range(6)]
    a = matrix(rows, (6, 6), fmt)
    reference = matrix(rows, (6, 6), fmt, sp)
    calls = []
    old = sx._core._dm_clear_denoms
    def spy(*args):
        result = old(*args)
        calls.append(result is not None)
        return result
    monkeypatch.setattr(sx._core, '_dm_clear_denoms', spy)
    actual = getattr(a, entry)(method=algorithm)
    expected = getattr(reference, entry)(method=algorithm)
    if entry == 'rref':
        b, pivots = actual
        b_ref, pivots_ref = expected
    else:
        b, den, pivots = actual
        b_ref, den_ref, pivots_ref = expected
        b = b.to_field() / den
        b_ref = b_ref.to_field() / den_ref
    assert b.to_Matrix() == sx.Matrix([[sx.Rational(str(v)) for v in row] for row in b_ref.to_list()])
    assert pivots == pivots_ref and calls == [True]
