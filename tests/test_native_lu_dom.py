"""Native exact LU, triangular solves and Gauss-Jordan beyond QQ.

GF(p) and QQ<a> run ddm_ilu / ddm_ilu_solve natively; GF(p), QQ<a> and QQ_I
run dense ddm_irref and sparse sdm_irref natively. Expected values are
hard-coded from SymPy 1.14 or compared with the reference SymPy functions on
the same inputs, and every result is also checked by an exact certificate.
"""
import copy
import importlib
import random

import pytest
import sympy as sp
import symox as sx
from symox.polys.matrices import DomainMatrix
from symox.polys.matrices.exceptions import DMNonInvertibleMatrixError
from sympy.polys.matrices import DomainMatrix as RefDomainMatrix

dense = importlib.import_module('symox.polys.matrices.dense')
sdm = importlib.import_module('symox.polys.matrices.sdm')


def symox_field(name):
    return {'GF7': sx.GF(7), 'GF10007': sx.GF(10007), 'QQ_I': sx.QQ_I,
            'sqrt2': sx.QQ.algebraic_field(sx.sqrt(2)),
            'cbrt2': sx.QQ.algebraic_field(sx.cbrt(2))}[name]


def sympy_field(name):
    return {'GF7': sp.GF(7), 'GF10007': sp.GF(10007), 'QQ_I': sp.QQ_I,
            'sqrt2': sp.QQ.algebraic_field(sp.sqrt(2)),
            'cbrt2': sp.QQ.algebraic_field(sp.cbrt(2))}[name]


def entry(name, rng, zero_rate):
    if rng.random() < zero_rate:
        return 0
    if name.startswith('GF'):
        return rng.randrange(1, 7 if name == 'GF7' else 10007)
    if name == 'QQ_I':
        return sp.Rational(rng.randint(-4, 4), rng.randint(1, 3)) + rng.randint(-3, 3) * sp.I
    root = sp.sqrt(2) if name == 'sqrt2' else sp.cbrt(2)
    return sp.Rational(rng.randint(-4, 4), rng.randint(1, 3)) + rng.randint(-3, 3) * root + (
        rng.randint(-2, 2) * root**2 if name == 'cbrt2' else 0)


def pair(name, rows):
    K, R = symox_field(name), sympy_field(name)
    if name.startswith('GF'):
        a = [[K(int(v)) for v in row] for row in rows]
        b = [[R(int(v)) for v in row] for row in rows]
    else:
        a = [[K.from_sympy(sx.sympify(str(v))) for v in row] for row in rows]
        b = [[R.from_sympy(sp.sympify(v)) for v in row] for row in rows]
    m, n = len(rows), len(rows[0]) if rows else 0
    return DomainMatrix(a, (m, n), K), RefDomainMatrix(b, (m, n), R)


def text(M):
    return [str(sp.sympify(str(v))) for v in M.to_Matrix()]


def ref_text(M):
    return [str(v) for v in M.to_Matrix()]


def record(monkeypatch, name):
    native = getattr(sx._core, name)
    calls = []

    def wrapper(*args):
        result = native(*args)
        calls.append(result)
        return result
    monkeypatch.setattr(sx._core, name, wrapper)
    return calls


def test_hard_coded_sympy_gf7_lu_solve_and_rref(monkeypatch):
    calls = record(monkeypatch, '_ddm_ilu')
    F = sx.GF(7)
    M = DomainMatrix([[F(0), F(2), F(3)], [F(4), F(5), F(6)], [F(1), F(1), F(1)]], (3, 3), F)
    L, U, swaps = M.lu()
    assert L.to_Matrix().tolist() == [[1, 0, 0], [0, 1, 0], [2, -1, 1]]
    assert U.to_Matrix().tolist() == [[-3, -2, -1], [0, 2, 3], [0, 0, -1]]
    assert swaps == [(0, 1)] and calls and calls[0] == [(0, 1)]
    assert M.lu_solve(DomainMatrix([[F(1)], [F(2)], [F(3)]], (3, 1), F)).to_Matrix().tolist() == [[-1], [-3], [0]]
    R, pivots = M.rref(method='GJ')
    assert R.to_Matrix().tolist() == [[1, 0, 0], [0, 1, 0], [0, 0, 1]] and pivots == (0, 1, 2)
    singular = DomainMatrix([[F(0), F(2), F(3)], [F(4), F(5), F(6)], [F(1), F(1), F(2)]], (3, 3), F)
    with pytest.raises(DMNonInvertibleMatrixError):
        singular.lu_solve(DomainMatrix([[F(1)], [F(2)], [F(3)]], (3, 1), F))


def test_hard_coded_sympy_sqrt2_lu_solve_and_rref(monkeypatch):
    calls = record(monkeypatch, '_ddm_ilu_solve')
    A = sx.QQ.algebraic_field(sx.sqrt(2))
    r = sx.sqrt(2)
    N = DomainMatrix([[A.from_sympy(r), A(1), A(0)], [A(2), A.from_sympy(1 + r), A(3)],
                      [A(1), A(0), A.from_sympy(r)]], (3, 3), A)
    L, U, swaps = N.lu()
    assert [str(v) for v in L.to_Matrix()] == ['1', '0', '0', 'sqrt(2)', '1', '0', 'sqrt(2)/2', '-sqrt(2)/2', '1']
    assert [str(v) for v in U.to_Matrix()] == ['sqrt(2)', '1', '0', '0', '1', '3', '0', '0', '5*sqrt(2)/2']
    assert swaps == []
    x = N.lu_solve(DomainMatrix([[A(1)], [A(0)], [A(1)]], (3, 1), A))
    assert [str(v) for v in x.to_Matrix()] == ['sqrt(2)/5 + 1', '3/5 - sqrt(2)', '-1/5']
    assert calls == [True]
    assert [str(v) for v in N.rref(method='GJ')[0].to_Matrix()] == ['1', '0', '0', '0', '1', '0', '0', '0', '1']


@pytest.mark.parametrize('sparse', [False, True])
def test_hard_coded_sympy_gaussian_rref(sparse, monkeypatch):
    calls = record(monkeypatch, '_sdm_irref' if sparse else '_ddm_irref')
    Q = DomainMatrix([[sx.QQ_I(0, 1), sx.QQ_I(2, 0), sx.QQ_I(1, 1)],
                      [sx.QQ_I(1, 0), sx.QQ_I(0, 0), sx.QQ_I(3, -1)]], (2, 3), sx.QQ_I)
    R, pivots = Q.rref(method='GJ' if sparse else 'GJ_dense')
    assert [str(v) for v in R.to_Matrix()] == ['1', '0', '3 - I', '0', '1', '-I'] and pivots == (0, 1)
    assert calls and calls[0] is not None


SHAPES = [(1, 1), (2, 2), (3, 3), (4, 4), (3, 5), (5, 3), (6, 6)]


@pytest.mark.parametrize('name', ['GF7', 'GF10007', 'sqrt2', 'cbrt2'])
@pytest.mark.parametrize('shape', SHAPES)
@pytest.mark.parametrize('seed', range(3))
def test_lu_matches_reference_and_certifies_pa_equals_lu(name, shape, seed, monkeypatch):
    calls = record(monkeypatch, '_ddm_ilu')
    rng = random.Random(hash((name, shape, seed)) & 0xffff)
    rows = [[entry(name, rng, 0.35) for _ in range(shape[1])] for _ in range(shape[0])]
    M, R = pair(name, rows)
    L, U, swaps = M.to_dense().lu()
    RL, RU, rswaps = R.to_dense().lu()
    assert swaps == rswaps
    assert text(L) == ref_text(RL) and text(U) == ref_text(RU)
    assert calls and calls[0] == swaps
    permuted = M.to_dense().to_list()
    for i, j in swaps:
        permuted[i], permuted[j] = permuted[j], permuted[i]
    assert L * U == DomainMatrix(permuted, shape, M.domain)


@pytest.mark.parametrize('name', ['GF7', 'GF10007', 'sqrt2'])
@pytest.mark.parametrize('n,cols', [(1, 1), (2, 1), (3, 2), (5, 1), (6, 3)])
@pytest.mark.parametrize('seed', range(3))
def test_lu_solve_matches_reference_and_certifies_ax_equals_b(name, n, cols, seed, monkeypatch):
    calls = record(monkeypatch, '_ddm_ilu_solve')
    rng = random.Random(hash((name, n, cols, seed)) & 0xffff)
    rows = [[entry(name, rng, 0.2) for _ in range(n)] for _ in range(n)]
    rhs = [[entry(name, rng, 0.2) for _ in range(cols)] for _ in range(n)]
    M, R = pair(name, rows)
    B, RB = pair(name, rhs)
    try:
        expected = R.lu_solve(RB)
    except Exception as error:
        with pytest.raises(Exception) as raised:
            M.lu_solve(B)
        assert type(raised.value).__name__ == type(error).__name__
        assert str(raised.value) == str(error)
        return
    x = M.lu_solve(B)
    assert text(x) == ref_text(expected)
    assert M * x == B
    assert calls == [True]


@pytest.mark.parametrize('name', ['GF7', 'GF10007', 'sqrt2', 'cbrt2', 'QQ_I'])
@pytest.mark.parametrize('shape', SHAPES + [(4, 7), (7, 4)])
@pytest.mark.parametrize('sparse', [False, True])
@pytest.mark.parametrize('seed', range(2))
def test_gauss_jordan_rref_matches_reference(name, shape, sparse, seed, monkeypatch):
    calls = record(monkeypatch, '_sdm_irref' if sparse else '_ddm_irref')
    rng = random.Random(hash((name, shape, sparse, seed)) & 0xffff)
    rows = [[entry(name, rng, 0.4) for _ in range(shape[1])] for _ in range(shape[0])]
    M, R = pair(name, rows)
    method = 'GJ' if sparse else 'GJ_dense'
    reduced, pivots = M.rref(method=method)
    expected, rpivots = R.rref(method=method)
    assert pivots == rpivots and text(reduced) == ref_text(expected)
    if any(v for row in rows for v in row):
        assert any(result is not None for result in calls)


@pytest.mark.parametrize('name', ['GF7', 'sqrt2', 'QQ_I'])
def test_low_level_dense_kernels_publish_the_reference_row_objects(name):
    rng = random.Random(4242)
    rows = [[entry(name, rng, 0.3) for _ in range(4)] for _ in range(4)]
    rows[0][0] = 0
    M, R = pair(name, rows)
    a = M.to_dense().rep.to_ddm().to_list()
    r = R.to_dense().rep.to_ddm().to_list()
    original_rows = list(a)
    pivots = dense.ddm_irref(a)
    from sympy.polys.matrices.dense import ddm_irref
    assert pivots == ddm_irref(r)
    assert [[str(M.domain.to_sympy(v)) for v in row] for row in a] == \
        [[str(R.domain.to_sympy(v)) for v in row] for row in r]
    assert sorted(map(id, a)) == sorted(map(id, original_rows))


def test_partial_pivoting_and_mixed_domains_are_declined_without_mutation():
    F, G = sx.GF(7), sx.GF(11)
    a = [[F(0), F(2)], [F(3), F(4)]]
    before = copy.deepcopy(a)
    assert sx._core._ddm_irref(a, True) is None and a == before
    mixed = [[F(1), F(2)], [G(3), G(4)]]
    assert sx._core._ddm_ilu(mixed) is None
    assert sx._core._ddm_irref(mixed, False) is None
    A2 = sx.QQ.algebraic_field(sx.sqrt(2))
    A3 = sx.QQ.algebraic_field(sx.sqrt(3))
    fields = [[A2(1), A2(2)], [A3(3), A3(4)]]
    assert sx._core._ddm_ilu(fields) is None
    assert sx._core._ddm_irref(fields, False) is None


def test_gaussian_lu_and_solve_keep_the_python_wrapper_route(monkeypatch):
    # QQ_I division goes through the Python Gaussian wrapper (convert/exquo
    # callbacks): LU and triangular solves stay in Python and stay correct.
    calls = record(monkeypatch, '_ddm_ilu')
    Q = DomainMatrix([[sx.QQ_I(0, 1), sx.QQ_I(2, 0)], [sx.QQ_I(1, 0), sx.QQ_I(3, -1)]], (2, 2), sx.QQ_I)
    L, U, swaps = Q.lu()
    assert calls == [None]
    assert [str(v) for v in L.to_Matrix()] == ['1', '0', '-I', '1'] and swaps == []
    assert [str(v) for v in U.to_Matrix()] == ['I', '2', '0', '3 + I']


def test_sparse_explicit_zero_and_foreign_values_decline():
    F = sx.GF(7)
    assert sx._core._sdm_irref({0: {0: F(0), 1: F(1)}}) is None
    assert sx._core._sdm_irref({0: {0: F(1)}, 1: {0: sx.GF(11)(1)}}) is None
    assert sx._core._sdm_irref({0: {0: F(1)}, 1: {1: 3}}) is None


def test_native_interruption_free_large_prime_field_is_exact():
    p = 2**31 - 1
    F = sx.GF(p)
    rng = random.Random(7)
    rows = [[rng.randrange(p) for _ in range(5)] for _ in range(5)]
    M = DomainMatrix([[F(v) for v in row] for row in rows], (5, 5), F)
    R = RefDomainMatrix([[sp.GF(p)(v) for v in row] for row in rows], (5, 5), sp.GF(p))
    L, U, swaps = M.lu()
    RL, RU, rswaps = R.lu()
    assert swaps == rswaps and L.to_Matrix() == sx.sympify(str(RL.to_Matrix()))
    assert U.to_Matrix() == sx.sympify(str(RU.to_Matrix()))
