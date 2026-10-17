"""Structural transforms: SymPy conventions and independent inverse certificates."""
import itertools
import pickle
import random
import subprocess
import sys

import pytest
import sympy as sp
from sympy.combinatorics import Permutation
import symox as sx

CLASSES = [sx.Matrix, sx.ImmutableMatrix, sx.SparseMatrix, sx.ImmutableSparseMatrix]


def entries(matrix):
    return [[matrix[i, j] for j in range(matrix.cols)] for i in range(matrix.rows)]


def check(actual, expected):
    assert actual.shape == expected.shape
    assert all(sp.expand(actual[i, j].to_sympy() - expected[i, j]) == 0
               for i in range(actual.rows) for j in range(actual.cols))


@pytest.mark.parametrize('cls', CLASSES)
@pytest.mark.parametrize('orientation', ['rows', 'cols', 'columns'])
@pytest.mark.parametrize('direction', ['forward', 'backward', 'forwards', 'backwards'])
@pytest.mark.parametrize('perm', [[], [1, 0], [2, 0, 1], [[0, 1], [0, 2]], [[0, 1, 2]], [[0], []]])
def test_permutation_conventions_and_aliases(cls, orientation, direction, perm):
    original = cls(3, 3, list(range(9)))
    result = original.permute(perm, orientation, direction)
    check(result, sp.Matrix(3, 3, range(9)).permute(perm, orientation, direction))
    assert entries(original) == [[0, 1, 2], [3, 4, 5], [6, 7, 8]]
    assert result._immutable == original._immutable
    if orientation == 'rows':
        assert result == original.permute_rows(perm, direction)
        assert original.permuteFwd(perm) == original.permute(perm)
        assert original.permuteBkwd(perm) == original.permute(perm, direction='backward')
    else:
        assert result == original.permute_cols(perm, direction)


@pytest.mark.parametrize('cls', CLASSES)
def test_rectangular_inverse_certificates_and_foreign_permutation_protocol(cls):
    rng = random.Random(104)
    for rows, cols in [(1, 4), (4, 1), (3, 5), (5, 3)]:
        original = cls(rows, cols, list(range(rows*cols)))
        for orientation, size in [('rows', rows), ('cols', cols)]:
            array = list(range(size)); rng.shuffle(array)
            inverse = [array.index(i) for i in range(size)]
            transformed = original.permute(array, orientation)
            assert transformed.permute(inverse, orientation) == original
            assert original.permute(Permutation(array), orientation) == transformed
    m = cls([[1], [2], [3]])
    # Backward reverses a sequence of cycles, not each individual cycle or array.
    assert m.permuteBkwd([1, 2, 0]) == m.permuteFwd([1, 2, 0])
    assert m.permuteFwd([[0, 1], [0, 2]]).permuteBkwd([[0, 1], [0, 2]]) == m
    assert m.permute_rows([1, 0, 2, 3]) == m.permute_rows([1, 0])


@pytest.mark.parametrize('cls', CLASSES)
@pytest.mark.parametrize('axis', ['row', 'col'])
@pytest.mark.parametrize('op,kwargs', [
    ('n->kn', dict(index=0, k=sp.Rational(2, 3))),
    ('n->kn', dict(first=1, k=0)),
    ('n<->m', dict(first=0, second=2)),
    ('n<->m', dict(index=0, k=2)),
    ('n<->m', dict(first=0, second=2, k=sp.Symbol('irrelevant'))),
    ('n->n+km', dict(first=0, second=2, k=-3)),
    ('n->n+km', dict(index=0, first=2, k=sp.Symbol('k'))),
])
def test_elementary_operations_match_sympy_and_preserve_type(cls, axis, op, kwargs):
    mapping = dict(index=axis, first=axis+'1', second=axis+'2')
    kwargs = {mapping.get(k, k): v for k, v in kwargs.items()}
    native_kwargs = {k: sx.sympify(v) if k == 'k' else v for k, v in kwargs.items()}
    original = cls(3, 3, range(9))
    result = getattr(original, 'elementary_'+axis+'_op')(op, **native_kwargs)
    check(result, getattr(sp.Matrix(3, 3, range(9)), 'elementary_'+axis+'_op')(op, **kwargs))
    assert result._immutable == original._immutable
    assert pickle.loads(pickle.dumps(result)) == result


@pytest.mark.parametrize('cls', CLASSES)
@pytest.mark.parametrize('axis', ['row', 'col'])
def test_elementary_inverse_certificates_rectangles_and_scalar_scaling(cls, axis):
    a, k = sx.symbols('a k')
    original = cls([[a, 2, 0], [3, a, 1]])
    name = 'elementary_'+axis+'_op'
    kw = {axis: 0, axis+'2': 1, 'k': 3}
    transformed = getattr(original, name)('n->n+km', **kw)
    kw['k'] = -3
    recovered = getattr(transformed, name)('n->n+km', **kw)
    assert recovered.applyfunc(sx.expand) == original
    scaled = getattr(original, name)('n->kn', **{axis: 0, 'k': k})
    assert scaled[0, 0] == k*a
    added = getattr(original, name)('n->n+km', **{axis: 0, axis+'2': 1, 'k': k})
    i, j = (0, 1) if axis == 'row' else (1, 0)
    assert added[i, j] == original[i, j]+k*a


@pytest.mark.parametrize('cls', CLASSES)
@pytest.mark.parametrize('shape', [(0, 0), (0, 4), (4, 0)])
def test_empty_shapes_preserved(cls, shape):
    m = cls(*shape, [])
    for orientation in ['rows', 'cols']:
        assert m.permute([], orientation).shape == shape
        dimension = shape[0 if orientation == 'rows' else 1]
        if dimension:
            axis = 'row' if orientation == 'rows' else 'col'
            assert getattr(m, 'elementary_'+axis+'_op')(**{axis: 0, 'k': 2}).shape == shape


@pytest.mark.parametrize('cls', CLASSES)
@pytest.mark.parametrize('perm,error', [(1, ValueError), ([0, 0], ValueError), ([1], ValueError),
    ([-1, 0], IndexError), ([3, 1, 2, 0], IndexError), ([[0, 0]], ValueError),
    ([[0, 3]], IndexError), ([0.0, 1], TypeError)])
def test_invalid_permutations(cls, perm, error):
    with pytest.raises(error): cls.eye(3).permute(perm)


@pytest.mark.parametrize('cls', CLASSES)
def test_invalid_options_and_elementary_signatures(cls):
    m = cls.eye(3)
    for kwargs in [dict(orientation='diagonal'), dict(direction='sideways')]:
        with pytest.raises(TypeError): m.permute([], **kwargs)
    for name in ['elementary_row_op', 'elementary_col_op']:
        axis = 'row' if 'row' in name else 'col'
        for op, kwargs in [('bad', {}), ('n->kn', {}), ('n->kn', {axis: -1, 'k': 2}),
            ('n->kn', {axis: 3, 'k': 2}), ('n<->m', {axis+'1': 0, axis+'2': 0}),
            ('n->n+km', {axis+'1': 0, axis+'2': 0, 'k': 2})]:
            with pytest.raises(ValueError): getattr(m, name)(op, **kwargs)
        with pytest.raises(TypeError): getattr(m, name)(**{axis: 0.5, 'k': 2})


@pytest.mark.parametrize('cls', [sx.Matrix, sx.SparseMatrix])
def test_callbacks_are_sequential_allow_self_reads_and_support_negative_indices(cls):
    original = [[1, 2, 0], [3, 4, 5], [6, 7, 8]]
    m, oracle = cls(original), sp.Matrix(original)
    assert m.row_op(-1, lambda v, j: v+2*m[0, j]) is None
    oracle.row_op(-1, lambda v, j: v+2*oracle[0, j])
    check(m, oracle)
    assert m.col_op(-1, lambda v, i: v-m[i, 0]) is None
    oracle.col_op(-1, lambda v, i: v-oracle[i, 0])
    assert m.zip_row_op(1, 0, lambda a, b: a-3*b) is None
    oracle.zip_row_op(1, 0, lambda a, b: a-3*b)
    check(m, oracle)
    m.row_op(0, lambda v, j: v+sum(m[0, k] for k in range(j)))
    oracle.row_op(0, lambda v, j: v+sum(oracle[0, k] for k in range(j)))
    check(m, oracle)


@pytest.mark.parametrize('cls', [sx.Matrix, sx.SparseMatrix])
def test_callback_failure_keeps_earlier_updates_and_empty_callbacks_are_not_invoked(cls):
    m = cls([[1, 2, 3]])
    def update(value, j):
        if j == 1: raise RuntimeError('expected')
        return value+10
    with pytest.raises(RuntimeError, match='expected'): m.row_op(0, update)
    assert entries(m) == [[11, 2, 3]]
    assert cls(3, 0, []).row_op(999, None) is None
    assert cls(0, 3, []).col_op(999, None) is None
    assert cls(0, 0, []).zip_row_op(999, 888, None) is None
    with pytest.raises(IndexError): cls(0, 3, []).row_op(0, update)


@pytest.mark.parametrize('cls', [sx.ImmutableMatrix, sx.ImmutableSparseMatrix])
def test_immutable_callbacks_refused_before_calling_function(cls):
    m = cls.eye(2)
    def fail(*args): raise AssertionError('must not call')
    for name, args in [('row_op', (0, fail)), ('col_op', (0, fail)), ('zip_row_op', (0, 1, fail))]:
        with pytest.raises(TypeError): getattr(m, name)(*args)


@pytest.mark.parametrize('cls', CLASSES)
@pytest.mark.parametrize('coefficient', [sx.oo, -sx.oo, sx.zoo, sx.nan])
def test_nonfinite_coefficients_do_not_annihilate_implicit_zeros(cls, coefficient):
    m = cls([[0, 1, 0], [2, 0, 3]])
    for axis in ['row', 'col']:
        for op, kwargs in [('n->kn', {axis: 0, 'k': coefficient}),
                           ('n->n+km', {axis: 0, axis+'2': 1, 'k': coefficient})]:
            result = getattr(m, 'elementary_'+axis+'_op')(op, **kwargs)
            expected = getattr(sp.Matrix([[0, 1, 0], [2, 0, 3]]), 'elementary_'+axis+'_op')(op,
                **{k: v.to_sympy() if k == 'k' else v for k, v in kwargs.items()})
            assert entries(result) == [[sx.sympify(v) for v in row] for row in expected.tolist()]


def test_huge_sparse_permutations_and_elementary_ops_never_densify(monkeypatch):
    def fail(*args): raise AssertionError('must not densify')
    monkeypatch.setattr(sx.SparseMatrix, 'as_dense', fail)
    n = 10**50
    m = sx.SparseMatrix(n, n, {(0, 2): 3, (n-1, 3): 4})
    assert m.permute_rows([[0, n-1]]).todok() == {(n-1, 2): 3, (0, 3): 4}
    assert m.permute_cols([[2, 3]]).todok() == {(0, 3): 3, (n-1, 2): 4}
    assert m.elementary_row_op(row=0, k=2).todok() == {(0, 2): 6, (n-1, 3): 4}
    assert m.elementary_row_op('n->n+km', row=0, row2=n-1, k=2).todok() == {(0, 2): 3, (0, 3): 8, (n-1, 3): 4}
    for action in [lambda: m.elementary_row_op(row=0, k=sx.oo), lambda: m.row_op(0, fail),
                   lambda: m.col_op(0, fail), lambda: m.zip_row_op(0, n-1, fail)]:
        with pytest.raises(NotImplementedError, match='resource limit'): action()
    assert m.todok() == {(0, 2): 3, (n-1, 3): 4}


def test_permutation_protocol_checks_size_before_materializing_array():
    class Giant:
        size = 10**50
        @property
        def array_form(self): raise AssertionError('must check size first')
    with pytest.raises(NotImplementedError): sx.eye(3).permute(Giant())
    with pytest.raises(NotImplementedError): sx.eye(3).permute(itertools.repeat(0))


def test_native_arithmetic_rejects_latent_coefficients_before_evaluation_without_sympy():
    code = '''
import sys
sys.modules['sympy'] = None
import symox as s
one=s.Integer(1); y=s.Symbol('y')
raw=s._core._restore_arithmetic('Mul',s.Integer(2),[(y,one)],False)
for _ in range(80):
 raw=s._core._restore_arithmetic('Mul',one,[(raw,one),(raw,one)],False)
for cls in [s.Matrix,s.SparseMatrix]:
 for k in [raw,s.Pow(2,10**9,evaluate=False)]:
  try: cls([[1,0]]).elementary_row_op(row=0,k=k)
  except NotImplementedError: pass
  else: raise AssertionError('expected bounded rejection')
 a=cls([[1,2],[3,4]])
 assert a.permute_rows([[0,1]]).elementary_row_op('n->n+km',row=0,row2=1,k=-2)[0,0] == 1
'''
    subprocess.run([sys.executable, '-c', code], check=True, timeout=15)
