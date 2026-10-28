"""Differential and protocol tests for the native Smith decomposition (stock ZZ)
and the number-field integer HNF kernels.

Expected values are hard-coded from SymPy 1.14.0 (python ground types);
every Smith case is also certified by s*m*t == diag(invariants) with
unimodular s, t and the divisibility chain.
"""
import pytest

from symox import ZZ, QQ, Matrix, Poly, symbols, cyclotomic_poly, _core
from symox.polys.matrices import DomainMatrix
from symox.polys.matrices import normalforms as nf
from symox.polys.matrices import _smith_native
from symox.polys.matrices.normalforms import (smith_normal_form, smith_normal_decomp,
                                               invariant_factors)
from symox.polys.numberfields import _normalforms as nfh
from symox.polys.numberfields.basis import round_two
from symox.polys.numberfields.primes import prime_decomp

x = symbols('x')

SMITH = [
    ([[12, 6, 4], [3, 9, 6], [2, 16, 14]], (1, 10, 30), [[0, 1, -1], [1, -4, 0], [0, -2, 3]], [[1, 1, 10], [0, -1, -2], [0, 1, 3]]),
    ([[0, 0], [0, 0]], (0, 0), [[0, 1], [1, 0]], [[0, 1], [1, 0]]),
    ([[0, 0, 3]], (3,), [[1]], [[0, 0, 1], [0, 1, 0], [1, 0, 0]]),
    ([[0], [-4]], (4,), [[0, -1], [1, 0]], [[1]]),
    ([[0, 0], [0, 6], [4, 0]], (2, 12), [[0, 1, -1], [0, 2, -3], [1, 0, 0]], [[1, -3], [1, -2]]),
    ([[-6, 4, 0], [0, 0, 0], [10, -15, 21]], (1, 2, 0), [[-2, 0, -1], [-45, 0, -22], [0, 1, 0]], [[-3, 7, -42], [1, 10, -63], [0, 4, -25]]),
    ([[2, 4, 4], [-6, 6, 12], [10, -4, -16]], (2, 6, 12), [[1, 0, 0], [2, -1, -1], [3, -4, -3]], [[1, -2, 2], [0, 1, -2], [0, 0, 1]]),
    ([[0, 2, 0], [3, 0, 0], [0, 0, 5]], (1, 1, 30), [[-2, 0, 1], [13, 3, -7], [45, 10, -24]], [[1, -1, 10], [1, -5, 45], [1, -4, 36]]),
    ([[1, 2, 3, 4], [5, 6, 7, 8]], (1, 4), [[1, 0], [5, -1]], [[1, -2, 1, 2], [0, 1, -2, -3], [0, 0, 1, 0], [0, 0, 0, 1]]),
    ([[1, 2], [3, 4], [5, 6], [7, 8]], (1, 2), [[1, 0, 0, 0], [3, -1, 0, 0], [1, -2, 1, 0], [2, -3, 0, 1]], [[1, -2], [0, 1]]),
    ([[2, 0], [0, 3]], (1, 6), [[-1, 1], [-3, 2]], [[1, -3], [1, -2]]),
    ([[6, 0, 0], [0, 10, 0], [0, 0, 15]], (1, 30, 30), [[1, -1, 1], [5, -6, 6], [0, -3, 2]], [[1, -5, 0], [-1, 6, -3], [-1, 6, -2]]),
    ([[0, 0, 0], [0, 0, 0], [0, 0, 7]], (7, 0, 0), [[0, 0, 1], [0, 1, 0], [1, 0, 0]], [[0, 0, 1], [0, 1, 0], [1, 0, 0]]),
    ([[1000000000000000000000000000007, 30000000000000000000000000], [20000000000000000000000000000, -9999999999999999999999999999999]], (1, 10000000600000000000000000000068999999999999999999999999999993), [[-2857142857142857142857142857, 142857142857142857142857142851], [20000000000000000000000000000, -1000000000000000000000000000007]], [[1, 1428571514285714285714285714224142852857142857142857142857149], [0, 1]]),
    ([[1, 2, 3], [2, 4, 6], [3, 6, 9]], (1, 0, 0), [[1, 0, 0], [-3, 0, 1], [-2, 1, 0]], [[1, -3, -2], [0, 0, 1], [0, 1, 0]]),
    ([[-1]], (1,), [[-1]], [[1]]),
    ([[0]], (0,), [[1]], [[1]]),
    ([[5, -10]], (5,), [[1]], [[1, 2], [0, 1]]),
    ([[0, 34, 0, 23, 0], [0, -11, -13, -36, 0], [0, 5, 0, 0, 0], [-36, 0, 0, 0, 0], [-22, 10, 0, 0, 0]], (1, 1, 1, 2990, 0), [[0, -1, -2, 0, 0], [1, -1, -9, 0, 0], [-1269, -230, 8422, -2241, 3735], [-2540, -460, 16858, -4485, 7475], [0, 0, -36, -11, 18]], [[26, 1124, -7474, 11166155, 0], [1, 45, -299, 446706, 0], [0, -9, 59, -88146, 0], [0, 2, -13, 19422, 0], [0, 0, 0, 0, 1]]),
    ([[-1, 0, 0, 0], [0, 0, 23, 0], [0, 0, 0, 0], [0, 0, 0, 0]], (1, 23, 0, 0), [[-1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]], [[1, 0, 0, 0], [0, 0, 0, 1], [0, 1, 0, 0], [0, 0, 1, 0]]),
    ([[0, 0, 4, 0], [0, 11, 0, -1], [0, -27, 0, 36], [0, 5, 0, 0], [0, -24, 0, 0]], (1, 1, 4, 0), [[0, 5, 2, 0, 0], [0, -37, -1, 76, 0], [1, -168, -12, 304, 0], [0, 59616, 1656, -122208, 1], [0, -180, -5, 369, 0]], [[0, 0, 0, 1], [1, -67, 0, 0], [1, -1, 1, 0], [0, 1, 0, 0]]),
    ([[-25, 0, 0], [0, 0, 0]], (25, 0), [[-1, 0], [0, 1]], [[1, 0, 0], [0, 0, 1], [0, 1, 0]]),
    ([[27], [0]], (27,), [[1, 0], [0, 1]], [[1]]),
    ([[-24, 0, 0]], (24,), [[-1]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
    ([[0, -17], [3, 35], [0, 0], [-40, 0], [0, 0]], (1, 1), [[0, -13, 0, -1, 0], [247, 120, 0, 9, 0], [0, 0, 1, 0, 0], [-1400, -680, 0, -51, 0], [0, 0, 0, 0, 1]], [[1, 455], [0, 1]]),
    ([[19, 0], [20, -38], [-23, 0], [0, 0]], (1, 38), [[-1, 1, 0, 0], [-5, -1, -5, 0], [-23, 0, -19, 0], [0, 0, 0, 1]], [[1, 38], [0, 1]]),
    ([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 6], [0, 0, 0, 0]], (6, 0, 0, 0), [[0, 0, 1, 0], [0, 0, 0, 1], [0, 1, 0, 0], [1, 0, 0, 0]], [[0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0]]),
    ([[14, -5, 30, 25], [-33, 0, 0, 16], [0, 0, 0, 0]], (1, 1, 0), [[-7, -3, 0], [-33, -14, 0], [0, 0, 1]], [[1, -7, -80, 0], [0, -89, -1049, 6], [0, 0, 0, 1], [0, -14, -165, 0]]),
    ([[0, 0, -6], [0, 18, 0], [0, -17, 0], [11, -40, -1], [1, -10, 0]], (1, 1, 1), [[0, 0, 0, 0, 1], [0, 1, 1, 0, 0], [0, 70, 70, -1, 11], [0, -17, -18, 0, 0], [-1, -420, -420, 6, -66]], [[1, 10, 0], [0, 1, 0], [0, 0, 1]]),
    ([[0, 0, 0, -4, 40]], (4,), [[-1]], [[0, 0, 0, 1, 0], [0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [1, 0, 0, 0, 10], [0, 0, 0, 0, 1]]),
    ([[38, 0, 19], [29, 0, 0], [0, 9, 0]], (1, 1, 4959), [[13, -17, 0], [29, -38, 245], [261, -342, 2204]], [[1, 988, -544635], [0, 1, -551], [0, -4, 2205]]),
    ([[0, 0, 0, 0, 19], [0, 0, 0, 0, 0], [0, 0, 27, 0, 0]], (1, 513, 0), [[10, 0, 1], [189, 0, 19], [0, 1, 0]], [[0, 0, 0, 1, 0], [0, 0, 0, 0, 1], [-7, 190, 0, 0, 0], [0, 0, 1, 0, 0], [1, -27, 0, 0, 0]]),
    ([[20, -17], [0, 0], [0, 0]], (1, 0), [[1, 0, 0], [0, 0, 1], [0, 1, 0]], [[6, -17], [7, -20]]),
    ([[0, 0, 0, 39], [0, 0, 38, 0]], (1, 1482), [[1, 1], [38, 39]], [[0, 0, 0, 1], [0, 0, 1, 0], [-1, 39, 0, 0], [1, -38, 0, 0]]),
    ([[0, 0, -15, 0], [0, 0, 0, 0]], (15, 0), [[-1, 0], [0, 1]], [[0, 1, 0, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 0, 1, 0]]),
    ([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]], (0, 0, 0), [[0, 0, 1], [0, 1, 0], [1, 0, 0]], [[0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0]]),
    ([[15, 0], [0, -11], [0, 0], [0, -1], [0, 0]], (1, 15), [[0, 0, 0, -1, 0], [-1, 0, 0, -15, 0], [0, 0, 1, 0, 0], [0, -1, 0, 11, 0], [0, 0, 0, 0, 1]], [[1, -1], [1, 0]]),
    ([[39, 0, 0]], (39,), [[1]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
    ([[-5, 0, 0, 27], [0, 0, 1, 0]], (1, 1), [[1, 0], [0, 1]], [[-11, 0, 0, 27], [0, 0, 1, 0], [0, 1, 0, 0], [-2, 0, 0, 5]]),
    ([[0, 0, 0, -29]], (29,), [[-1]], [[0, 0, 0, 1], [0, 1, 0, 0], [0, 0, 1, 0], [1, 0, 0, 0]]),
    ([[0, 1], [0, -15], [14, 0]], (1, 14), [[1, 0, 0], [14, 0, -1], [15, 1, 0]], [[1, -1], [1, 0]]),
    ([[-32], [0], [0], [0], [0]], (32,), [[-1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]], [[1]]),
    ([[28, 0, 0], [39, 26, 0]], (1, 728), [[7, -5], [-39, 28]], [[1, 130, 0], [0, 1, 0], [0, 0, 1]]),
    ([[0, -27, 0], [0, 0, 0], [-38, 0, 0]], (1, 1026, 0), [[-1, 0, -5], [-38, 0, -189], [0, 1, 0]], [[1, -27, 0], [-7, 190, 0], [0, 0, 1]]),
    ([[0, 0, 0], [0, 0, 0], [0, -26, -22], [-29, 0, 36], [0, 0, 0]], (1, 2, 0), [[0, 0, 0, 1, 0], [0, 0, 1, -88, 0], [1, 0, 0, 0, 0], [0, 0, 0, 0, 1], [0, 1, 0, 0, 0]], [[-5, -72, 468], [0, 49, -319], [-4, -58, 377]]),
    ([[0, 34, 0, 33]], (1,), [[1]], [[0, 1, 0, 0], [1, 0, 0, 33], [0, 0, 1, 0], [-1, 0, 0, -34]]),
    ([[-18], [0], [21]], (3,), [[1, 0, 1], [0, 1, 0], [7, 0, 6]], [[1]]),
    ([[36, 12, 0, 0, 0], [0, 0, 0, 29, 0], [-21, 0, 0, 0, 0], [35, 0, 0, 0, 0]], (1, 1, 2436, 0), [[36, 0, 60, -1], [-70, 1, -120, 0], [-5887, 84, -10092, 0], [0, 0, -5, -3]], [[1, -432, 12528, 0, 0], [0, 1, -29, 0, 0], [0, 0, 0, 0, 1], [0, 29, -840, 0, 0], [0, 0, 0, 1, 0]]),
]

HNF = [
    ([[12, 6, 4], [3, 9, 6], [2, 16, 14]], None, [[10, 0, 2], [0, 15, 3], [0, 0, 2]]),
    ([[2, 7, 17, 29, 41], [3, 11, 19, 31, 43], [5, 13, 23, 37, 47]], None, [[1, 0, 0], [0, 2, 1], [0, 0, 1]]),
    ([[0, 0], [0, 0]], None, [[], []]),
    ([[-3, 5, 0], [6, -10, 0]], None, [[-1], [2]]),
    ([[100000000000000000003, 7], [11, 9999999999999999999]], None, [[999999999999999999929999999999999999920, 545454545454545454516363636363636363593], [0, 1]]),
    ([[-5, 11, -27], [-26, 22, 4], [-24, -7, 7]], None, [[9567, 7358, 1994], [0, 2, 0], [0, 0, 1]]),
    ([[-5, 11, -27], [-26, 22, 4], [-24, -7, 7]], 19134, [[9567, 7358, 1994], [0, 2, 0], [0, 0, 1]]),
    ([[-25, -3], [-4, -26]], None, [[319, 172], [0, 2]]),
    ([[-25, -3], [-4, -26]], 638, [[319, 172], [0, 2]]),
    ([[-3, -27, 22, 6, -23]], None, [[1]]),
    ([[-3, -27, 22, 6, -23]], 1, [[1]]),
    ([[7, -5, -27, -16, -28]], None, [[1]]),
    ([[7, -5, -27, -16, -28]], 3, [[1]]),
    ([[-4, -21, 4, -23], [6, -11, 5, 22]], None, [[1, 0], [0, 1]]),
    ([[-4, -21, 4, -23], [6, -11, 5, 22]], 3, [[1, 0], [0, 1]]),
    ([[7, 6], [10, -18]], None, [[93, 20], [0, 2]]),
    ([[7, 6], [10, -18]], 372, [[93, 20], [0, 2]]),
    ([[15, -26, 6, -27, 9]], None, [[1]]),
    ([[15, -26, 6, -27, 9]], 1, [[1]]),
    ([[19, -10, -1, 7, 29], [-1, -7, -11, -15, 20], [-19, 14, 19, -15, -25], [6, -11, 3, 1, 26]], None, [[2, 0, 1, 0], [0, 3, 1, 0], [0, 0, 1, 0], [0, 0, 0, 1]]),
    ([[19, -10, -1, 7, 29], [-1, -7, -11, -15, 20], [-19, 14, 19, -15, -25], [6, -11, 3, 1, 26]], 12, [[2, 0, 1, 0], [0, 3, 1, 0], [0, 0, 1, 0], [0, 0, 0, 1]]),
    ([[8, -26, -23, 2, -4], [-20, 18, -9, -21, 29], [1, -4, -28, 12, -26], [18, 5, 6, 20, 26]], None, [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]),
    ([[8, -26, -23, 2, -4], [-20, 18, -9, -21, 29], [1, -4, -28, 12, -26], [18, 5, 6, 20, 26]], 2, [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]),
    ([[-8, 8, 1, 7, 21], [-1, -26, 23, -25, 30], [-13, 0, 14, 12, -26]], None, [[3, 2, 0], [0, 1, 0], [0, 0, 1]]),
    ([[-8, 8, 1, 7, 21], [-1, -26, 23, -25, 30], [-13, 0, 14, 12, -26]], 3, [[3, 2, 0], [0, 1, 0], [0, 0, 1]]),
]
NF = [
    ('x**3 - 2', [[1, 0, 0], [0, 1, 0], [0, 0, 1]], 1, -108, {2: ['(2, x)'], 3: ['(3, x + 1)'], 5: ['(5, x + 2)', '(5, x**2 - 2*x - 1)'], 7: ['(7)']}),
    ('x**2 + 7', [[2, 1], [0, 1]], 2, -7, {2: ['(2, 3*x/2 + 1/2)', '(2, 3*x/2 + 3/2)'], 3: ['(3)'], 5: ['(5)'], 7: ['(7, x)']}),
    ('x**4 + x**3 + x**2 + x + 1', [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]], 1, 125, {2: ['(2)'], 3: ['(3)'], 5: ['(5, x - 1)'], 7: ['(7)']}),
    ('x**3 + x**2 - 2*x + 8', [[2, 0, 0], [0, 2, 1], [0, 0, 1]], 2, -503, {2: ['(2, x**2 + 1)', '(2, x**2/2 + 3*x/2 + 1)', '(2, 3*x**2/2 + 3*x/2)'], 3: ['(3)'], 5: ['(5, x + 1)', '(5, x**2 - 2)'], 7: ['(7)']}),
]


def _dm(rows):
    return DomainMatrix([[ZZ(v) for v in row] for row in rows], (len(rows), len(rows[0])), ZZ)


def _ints(dm):
    return [[int(v) for v in row] for row in dm.to_list()]


@pytest.mark.parametrize('rows,invs,s,t', SMITH)
def test_smith_matches_sympy_and_certificate(rows, invs, s, t):
    m = _dm(rows)
    assert _smith_native.try_smith(m.to_list(), ZZ, m.shape, True) is not None
    smf, S, T = smith_normal_decomp(m)
    assert tuple(int(v) for v in invariant_factors(m)) == invs
    assert _ints(S) == s and _ints(T) == t
    assert smith_normal_form(m) == smf.to_sparse() or smith_normal_form(m).to_dense() == smf
    assert S.domain == ZZ and T.domain == ZZ and smf.domain == ZZ
    for row in S.to_list() + T.to_list() + smf.to_list():
        assert all(type(v) is type(ZZ(1)) for v in row)
    assert S * m * T == smf
    assert abs(int(S.det())) == 1 and abs(int(T.det())) == 1
    for a, b in zip(invs, invs[1:]):
        assert (a == 0 and b == 0) or (a != 0 and b % a == 0)


def test_empty_shapes():
    for shape in ((0, 3), (3, 0), (0, 0)):
        m = DomainMatrix([[]] * shape[0] if shape[1] == 0 else [], shape, ZZ)
        assert invariant_factors(m) == ()
        smf, S, T = smith_normal_decomp(m)
        assert S.shape == (shape[0], shape[0]) and T.shape == (shape[1], shape[1])
        assert S.to_dense() == DomainMatrix.eye(shape[0], ZZ).to_dense()
        assert T.to_dense() == DomainMatrix.eye(shape[1], ZZ).to_dense()


def test_matrix_level_api():
    from symox.matrices.normalforms import smith_normal_form as snf, invariant_factors as invf
    M = Matrix([[12, 6, 4], [3, 9, 6], [2, 16, 14]])
    assert snf(M, domain=ZZ) == Matrix([[1, 0, 0], [0, 10, 0], [0, 0, 30]])
    assert invf(M, domain=ZZ) == (1, 10, 30)


def test_qq_and_non_pid_keep_reference_paths():
    m = DomainMatrix([[QQ(1, 2), QQ(3)], [QQ(4), QQ(0)]], (2, 2), QQ)
    # QQ is native since the 2026-10-29 port (tests/test_native_smith_fx.py).
    assert _smith_native.try_smith(m.to_list(), QQ, m.shape, False) == (QQ(1), QQ(1))
    assert invariant_factors(m) == (QQ(1), QQ(1))
    R = ZZ[x]
    p = DomainMatrix([[R(x), R(1)], [R(0), R(x)]], (2, 2), R)
    with pytest.raises(ValueError, match='principal ideal domain'):
        invariant_factors(p)
    with pytest.raises(ValueError, match='principal ideal domain'):
        smith_normal_decomp(p)


def test_python_int_entries_decline():
    rows = [[12, 6], [3, 9]]
    assert _smith_native.try_smith(rows, ZZ, (2, 2), True) is None


def _reference(rows, full):
    return nf._smith_normal_decomp([list(r) for r in rows], ZZ, shape=(len(rows), len(rows[0])), full=full)


def test_changed_protocol_declines_and_keeps_reference(monkeypatch):
    rows = SMITH[0][0]
    m = _dm(rows)
    original = nf.add_columns

    def add_columns(m, i, j, a, b, c, d):
        return original(m, i, j, a, b, c, d)
    monkeypatch.setattr(nf, 'add_columns', add_columns)
    assert _smith_native.try_smith(m.to_list(), ZZ, m.shape, True) is None
    assert tuple(int(v) for v in invariant_factors(m)) == SMITH[0][1]
    monkeypatch.undo()
    assert _smith_native.try_smith(m.to_list(), ZZ, m.shape, True) is not None


def test_changed_domain_method_declines(monkeypatch):
    m = _dm([[-4, 6], [2, 8]])
    native = _smith_native.try_smith(m.to_list(), ZZ, m.shape, True)
    assert native is not None
    calls = []
    original = _core.Domain.__dict__['canonical_unit']

    def canonical_unit(self, a):
        calls.append(a)
        return original(self, a)
    monkeypatch.setattr(_core.Domain, 'canonical_unit', canonical_unit)
    assert _smith_native.try_smith(m.to_list(), ZZ, m.shape, True) is None
    smf, S, T = smith_normal_decomp(m)
    assert calls
    assert (tuple(smf.to_list()[i][i] for i in range(2)), S.to_list(), T.to_list()) == native
    monkeypatch.undo()
    assert _smith_native.try_smith(m.to_list(), ZZ, m.shape, True) is not None


def test_added_class_binding_declines(monkeypatch):
    m = _dm([[2, 4], [6, 8]])
    monkeypatch.setattr(DomainMatrix, 'symox_smith_probe', 1, raising=False)
    assert _smith_native.try_smith(m.to_list(), ZZ, m.shape, False) is None
    monkeypatch.undo()
    assert _smith_native.try_smith(m.to_list(), ZZ, m.shape, False) is not None


def test_changed_code_declines():
    m = _dm([[2, 4], [6, 8]])
    function = nf.add_columns
    code = function.__code__
    try:
        function.__code__ = code.replace()
        assert _smith_native.try_smith(m.to_list(), ZZ, m.shape, False) is None
    finally:
        function.__code__ = code
    assert _smith_native.try_smith(m.to_list(), ZZ, m.shape, False) is not None


def test_large_random_against_python_reference():
    import random
    rng = random.Random(5)
    for _ in range(40):
        r, c = rng.randint(1, 6), rng.randint(1, 6)
        rows = [[rng.choice([0, rng.randint(-10**12, 10**12)]) for _ in range(c)] for _ in range(r)]
        m = _dm(rows)
        native = _smith_native.try_smith(m.to_list(), ZZ, m.shape, True)
        assert native is not None
        assert native == _reference(m.to_list(), True)
        assert _smith_native.try_smith(m.to_list(), ZZ, m.shape, False) == _reference(m.to_list(), False)


@pytest.mark.parametrize('rows,D,expected', HNF)
def test_numberfield_hnf_matches_sympy(rows, D, expected):
    m = _dm(rows)
    result = nfh.hermite_normal_form(m) if D is None else nfh.hermite_normal_form(m, D=D, check_rank=True)
    assert _ints(result) == expected
    assert result.domain == ZZ
    ints = [[int(v) for v in row] for row in rows]
    if D is None:
        assert nfh._nf_native('_nf_hnf', ints, m.shape) is not None
    else:
        assert nfh._nf_native('_nf_hnf_modular', ints, m.shape, D) is not None


def test_numberfield_hnf_fallback_and_errors(monkeypatch):
    m = _dm([[12, 6, 4], [3, 9, 6], [2, 16, 14]])
    ints = [[12, 6, 4], [3, 9, 6], [2, 16, 14]]
    original = nfh._gcdex
    monkeypatch.setattr(nfh, '_gcdex', lambda a, b: original(a, b))
    assert nfh._nf_native('_nf_hnf', ints, (3, 3)) is None
    assert _ints(nfh.hermite_normal_form(m)) == [[10, 0, 2], [0, 15, 3], [0, 0, 2]]
    monkeypatch.undo()
    with pytest.raises(nfh.DMShapeError):
        nfh.hermite_normal_form(_dm([[1], [2]]), D=ZZ(2))
    with pytest.raises(nfh.DMDomainError):
        nfh.hermite_normal_form(m, D=0)
    with pytest.raises(nfh.DMDomainError):
        nfh.hermite_normal_form(m.convert_to(QQ))
    assert nfh._nf_native('_nf_hnf', [[True, 2]], (1, 2)) is None


@pytest.mark.parametrize('T,basis,denom,disc,primes', NF)
def test_round_two_and_prime_decomp_match_sympy(T, basis, denom, disc, primes):
    from symox import sympify
    T = Poly(sympify(T), x)
    ZK, dK = round_two(T)
    assert [[int(v) for v in row] for row in ZK.matrix.to_list()] == basis
    assert int(ZK.denom) == denom and int(dK) == disc
    for p, expected in primes.items():
        assert [str(P) for P in prime_decomp(p, T=T)] == expected
