"""Exact row-lattice reduction, independently certified with Python Fractions."""
from fractions import Fraction
import pickle
import random
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx
from symox.polys.matrices.exceptions import DMError, DMRankError, DMShapeError, DMValueError


def rows(matrix):
    return [[int(matrix[i, j]) for j in range(matrix.cols)] for i in range(matrix.rows)]


def determinant(data):
    """Independent rational Gaussian determinant, distinct from native LLL."""
    data = [[Fraction(v) for v in row] for row in data]
    result = Fraction(1)
    for k in range(len(data)):
        pivot = next(i for i in range(k, len(data)) if data[i][k])
        if pivot != k:
            data[pivot], data[k] = data[k], data[pivot]
            result = -result
        value = data[k][k]
        result *= value
        for i in range(k+1, len(data)):
            factor = data[i][k]/value
            for j in range(k+1, len(data)):
                data[i][j] -= factor*data[k][j]
    return result


def certificate(original, basis, transform, delta=Fraction(3, 4)):
    m, n = original.shape
    assert basis.shape == (m, n) and transform.shape == (m, m)
    assert all(e.is_Integer is True for value in [basis, transform] for e in value)
    a, b, t = rows(original), rows(basis), rows(transform)
    assert [[sum(t[i][k]*a[k][j] for k in range(m)) for j in range(n)] for i in range(m)] == b
    assert abs(determinant(t)) == 1
    star, norms, mu = [], [], []
    for i in range(m):
        value = list(map(Fraction, b[i]))
        coefficients = []
        for j in range(i):
            coeff = sum(x*y for x, y in zip(b[i], star[j]))/norms[j]
            coefficients.append(coeff)
            value = [x-coeff*y for x, y in zip(value, star[j])]
        norm = sum(x*x for x in value)
        assert norm > 0
        assert all(abs(c) <= Fraction(1, 2) for c in coefficients)
        if i:
            assert norm >= (delta-coefficients[i-1]**2)*norms[i-1]
        star.append(value)
        norms.append(norm)
        mu.append(coefficients)


@pytest.mark.parametrize('shape', [(1, 1), (1, 4), (2, 2), (2, 4), (3, 5), (4, 6)])
@pytest.mark.parametrize('delta', [Fraction(1, 3), Fraction(1, 2), Fraction(3, 4), Fraction(99, 100)])
def test_deterministic_bases_match_sympy_and_independent_exact_certificates(shape, delta):
    m, n = shape
    rng = random.Random(1000*m+n)
    for _ in range(8):
        while True:
            data = [[rng.randrange(-30, 31) for j in range(n)] for i in range(m)]
            reference = sp.Matrix(data)
            if reference.rank() == m:
                break
        original = sx.Matrix(data)
        basis, transform = original.lll_transform(sx.Rational(delta.numerator, delta.denominator))
        expected_basis, expected_transform = reference.lll_transform(sp.Rational(delta.numerator, delta.denominator))
        assert rows(basis) == rows(expected_basis)
        assert rows(transform) == rows(expected_transform)
        assert original.lll(sx.Rational(delta.numerator, delta.denominator)) == basis
        assert rows(original) == data
        certificate(original, basis, transform, delta)


@pytest.mark.parametrize('immutable', [False, True])
def test_documented_classic_basis_mutability_pickle_and_interchange(immutable):
    cls = sx.ImmutableMatrix if immutable else sx.Matrix
    original = cls([[1, 0, 0, 0, -20160], [0, 1, 0, 0, 33768],
                    [0, 0, 1, 0, 39578], [0, 0, 0, 1, 47757]])
    basis, transform = original.lll_transform()
    expected = [[10, -3, -2, 8, -4], [3, -9, 8, 1, -11], [-3, 13, -9, -3, -9], [-12, -7, -11, 9, -1]]
    assert rows(basis) == expected
    assert rows(transform) == [r[:4] for r in expected]
    assert original.lll() == basis
    certificate(original, basis, transform)
    for value in [basis, transform]:
        assert isinstance(value, sx.ImmutableMatrix) is immutable
        assert pickle.loads(pickle.dumps(value)) == value
        assert sx.ImmutableMatrix(value.to_sympy()) == value


@pytest.mark.parametrize('coefficient,expected,transform', [(3, -1, -2), (-3, -1, 1), (1, 1, 0), (-1, -1, 0)])
def test_half_integer_ties_round_towards_positive_infinity(coefficient, expected, transform):
    original = sx.ImmutableMatrix([[2, 0], [coefficient, 2]])
    basis, t = original.lll_transform()
    assert basis == sx.ImmutableMatrix([[2, 0], [expected, 2]])
    assert t == sx.ImmutableMatrix([[1, 0], [transform, 1]])
    certificate(original, basis, t)


@pytest.mark.parametrize('delta', [.75, .3, .7, .3333333333333, sp.Float('.25123456789', 80), sp.Float('.3', 80)])
def test_float_delta_uses_sympy_rational_field_coercion(delta):
    original = sx.ImmutableMatrix([[15, 23, 11], [46, 15, 3], [9, 20, 29]])
    actual = original.lll_transform(sx.sympify(delta))
    expected = original.to_sympy().lll_transform(delta)
    assert tuple(rows(m) for m in actual) == tuple(rows(m) for m in expected)
    exact_delta = sp.polys.domains.QQ.from_sympy(sp.sympify(delta))
    certificate(original, *actual, Fraction(int(exact_delta.numerator), int(exact_delta.denominator)))


@pytest.mark.parametrize('delta', [sx.Rational(1, 4), 0, 1, -1, .25, 1., .9999999999999999])
def test_invalid_delta_range_raises_domain_matrix_value_error(delta):
    with pytest.raises(DMValueError, match='delta.*range'):
        sx.eye(2).lll(delta)
    with pytest.raises(DMValueError):
        sx.Matrix(0, 3, []).lll_transform(delta)


@pytest.mark.parametrize('delta', [True, False, sx.I, sx.sqrt(2), sx.Symbol('d'), sx.oo, sx.nan])
def test_nonrational_delta_rejected_without_guessing(delta):
    with pytest.raises(sx.CoercionFailed, match='delta.*rational'):
        sx.eye(2).lll(delta)


@pytest.mark.parametrize('entry', [sx.Rational(1, 2), sx.Symbol('x', integer=True), sx.I, sx.oo, sx.zoo, sx.nan, 1.5, True])
def test_noninteger_entry_domains_rejected(entry):
    for method in ['lll', 'lll_transform']:
        with pytest.raises(sx.CoercionFailed, match='ZZ'):
            getattr(sx.Matrix([[entry, 1]]), method)()


def test_integral_floating_entries_are_converted_by_actual_zz_coercion():
    original = sx.ImmutableMatrix([[sx.Float('2.0', 80), -3.0, 4], [1.0, sx.Float('4.0', 80), 2.]])
    assert original[0, 0].is_Integer is False
    basis, transform = original.lll_transform()
    expected = original.to_sympy().lll_transform()
    assert tuple(rows(m) for m in [basis, transform]) == tuple(rows(m) for m in expected)
    certificate(original, basis, transform)
    assert original[0, 0].is_Integer is False
    # Singleton bases need no Gram-Schmidt squaring; compact but bounded
    # integral floats and large canonical integers can retain this fast path.
    compact = sx.Float(2, 80)**10000
    single = sx.ImmutableMatrix([[compact, 1]])
    result, t = single.lll_transform()
    assert result[0, 0] == 2**10000 and t == sx.eye(1)
    certificate(single, result, t)


@pytest.mark.parametrize('shape', [(0, 0), (0, 7), (0, 10**6)])
@pytest.mark.parametrize('immutable', [False, True])
def test_empty_row_basis_preserves_width_and_empty_transform(shape, immutable):
    cls = sx.ImmutableMatrix if immutable else sx.Matrix
    original = cls(shape[0], shape[1], [])
    basis, transform = original.lll_transform()
    assert basis.shape == shape and transform.shape == (0, 0)
    assert isinstance(basis, sx.ImmutableMatrix) is immutable
    certificate(original, basis, transform)


@pytest.mark.parametrize('shape', [(1, 0), (3, 2), (3, 0)])
def test_shape_restrictions_are_explicit(shape):
    with pytest.raises(DMShapeError, match='m <= n'):
        sx.zeros(*shape).lll()


@pytest.mark.parametrize('data', [[[0]], [[0, 0, 0]], [[1, 2], [2, 4]], [[0, 0], [1, 2]], [[1, 2], [0, 0]]])
def test_dependent_rows_include_the_single_zero_row_rejected_mathematically(data):
    with pytest.raises(DMRankError, match='dependent rows'):
        sx.Matrix(data).lll_transform()


@pytest.mark.parametrize('immutable', [False, True])
def test_sparse_delegation_certificates_and_guards_before_dense_allocation(immutable, monkeypatch):
    cls = sx.ImmutableSparseMatrix if immutable else sx.SparseMatrix
    original = cls([[2, 0, 1], [3, 2, 0]])
    basis, transform = original.lll_transform()
    assert isinstance(basis, cls) and isinstance(transform, cls)
    assert original.lll() == basis
    certificate(original, basis, transform)
    assert cls(0, 10**6, {}).lll_transform()[0].shape == (0, 10**6)
    oversized = cls(100000, 100000, {})
    monkeypatch.setattr(cls, 'as_dense', lambda _: pytest.fail('densified before LLL preflight'))
    with pytest.raises(NotImplementedError, match='LLL resource limit'):
        oversized.lll()
    with pytest.raises(NotImplementedError, match='LLL resource limit'):
        cls(1, 10**50, {}).lll_transform()
    with pytest.raises(DMShapeError):
        cls(3, 2, {}).lll()
    with pytest.raises(DMValueError):
        oversized.lll(sx.Rational(1, 4))


def test_exception_names_and_bases_do_not_require_sympy():
    assert issubclass(DMValueError, DMError)
    assert issubclass(DMShapeError, DMError)
    assert issubclass(DMRankError, DMError)


def test_native_runtime_and_resource_guards_cover_held_dags_bits_and_sparse_shapes():
    subprocess.run([sys.executable, '-c', '''
import sys
sys.modules['sympy']=None
import symox as s
from symox.polys.matrices.exceptions import DMRankError
m=s.ImmutableMatrix([[1,0,0,0,-20160],[0,1,0,0,33768],[0,0,1,0,39578],[0,0,0,1,47757]])
b,t=m.lll_transform()
assert t*m==b and abs(int(t.det()))==1
try:
    s.zeros(1,3).lll()
except DMRankError:
    pass
else:
    raise AssertionError('single zero row is not a basis')
held=s.Pow(2,10**9,evaluate=False)
one=s.Integer(1)
y=s.Symbol('y')
raw=s._core._restore_arithmetic('Mul',s.Integer(2),[(y,one)],False)
for _ in range(80):
    raw=s._core._restore_arithmetic('Mul',one,[(raw,one),(raw,one)],False)
for entry in [held,raw,s.UnevaluatedExpr(2)]:
    try:
        s.ImmutableMatrix([[entry,1]]).lll()
    except s.CoercionFailed:
        pass
    else:
        raise AssertionError('held expression was evaluated during coercion')
try:
    s.eye(2).lll(held)
except s.CoercionFailed:
    pass
else:
    raise AssertionError('held delta was evaluated')
# A compact power-of-two float must be rejected before integer materialization.
f=s._core._restore_matrix_dag(1,[('s.number',True,(('mpfloat',1,10**9,80),))],0)
for matrix in [s.ImmutableMatrix([[f,1]]),s.ImmutableMatrix([[2**100000,1],[1,2]])]:
    try:
        matrix.lll()
    except NotImplementedError as error:
        assert 'LLL resource limit' in str(error)
    else:
        raise AssertionError('numeric bit/work preflight did not run')
'''], check=True, timeout=20)
