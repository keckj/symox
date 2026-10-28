"""Fraction-free kernels over GF(p), ZZ_I/QQ_I, QQ<a>, ZZ[x]/QQ[x]/QQ(x) and
composite rings: SymPy 1.14 values, exact certificates and identical results
with the native kernels disabled (the complete Python bodies)."""
import json
from pathlib import Path

import pytest

import symox as sx
from symox import _core
from symox.polys.domains import GF, QQ, ZZ, ZZ_I, QQ_I
from symox.polys.matrices import DomainMatrix
from symox.polys.matrices import dense
from symox.polys.matrices.ddm import DDM

x, y = sx.symbols('x y')
s2 = sx.sqrt(2)
I = sx.I
R = sx.Rational
A = QQ.algebraic_field(s2)

DOMAINS = {
    'ZZ_I': (ZZ_I, [[1+I, 2, I], [3, 1-I, 2+I], [I, 4, 5]]),
    'QQ_I': (QQ_I, [[R(1, 2)+I, 2, I], [3, 1-I, R(2, 3)+I], [I, 4, 5]]),
    'GF7': (GF(7), [[1, 2, 3], [4, 5, 6], [0, 1, 1]]),
    'QQs2': (A, [[1, s2, 2], [s2, 3, 1], [1, 1, s2]]),
    'ZZx': (ZZ[x], [[x, 1, 2], [x**2, x+1, 0], [1, x, 3]]),
    'QQx': (QQ[x], [[x/2, 1, 2], [x**2, x+1, 0], [1, x, R(3, 4)]]),
    'ZZxy': (ZZ[x, y], [[x, y, 2], [x*y, x+1, 0], [1, y, 3]]),
    'QQfx': (QQ.frac_field(x), [[1/x, 1, 2], [x**2, x+1, 0], [1, x, 3]]),
    'ZZIx': (ZZ_I[x], [[x+I, 1, 2], [x**2, x+1, 0], [1, I*x, 3]]),
    'QQs2x': (A[x], [[x+s2, 1, 2], [x**2, x+1, 0], [1, s2*x, 3]]),
}
VARIANTS = {
    'base': lambda r: r,
    'sing': lambda r: [r[0], r[1], [a+b for a, b in zip(r[0], r[1])]],
    'rect': lambda r: [r[0], r[1]],
    'tall': lambda r: [row[:2] for row in r],
    'zcol': lambda r: [[0]+row[1:] for row in r],
}
EXPECTED = json.loads((Path(__file__).parent/'data'/'native_fflu_dom_sympy114.json').read_text())


def matrix(name, variant):
    K, rows = DOMAINS[name]
    data = VARIANTS[variant](rows)
    return K, DomainMatrix([[K.from_sympy(sx.sympify(e)) for e in row] for row in data],
                           (len(data), len(data[0])), K)


def text(K, e):
    return str(sx.expand(K.to_sympy(e)))


def table(K, M):
    return [[text(K, e) for e in row] for row in M.to_list()]


def compute(K, M, op):
    def rref(R, den, pivots):
        return [table(K, R), text(K, den), list(pivots)]
    if op == 'fflu':
        return [table(K, X) for X in M.fflu()]
    if op == 'rref_den':
        return rref(*M.rref_den(method='FF'))
    if op == 'srref_den':
        return rref(*M.to_sparse().rref_den(method='FF'))
    if op == 'rref_den_auto':
        return rref(*M.rref_den())
    if op == 'det':
        return text(K, M.det())
    if op == 'ddet':
        return text(K, M.to_dense().det())
    if op == 'charpoly':
        return [text(K, c) for c in M.charpoly()]
    if op == 'scharpoly':
        return [text(K, c) for c in M.to_sparse().charpoly()]
    if op == 'berk':
        return [text(K, c) for c in M.to_dense().rep.charpoly()]
    if op == 'sberk':
        return [text(K, c) for c in M.to_sparse().rep.charpoly()]
    raise AssertionError(op)


@pytest.mark.parametrize('key', sorted(EXPECTED))
def test_public_results_match_sympy_114(key):
    name, variant, op = key.split('/')
    K, M = matrix(name, variant)
    assert compute(K, M, op) == EXPECTED[key]


NATIVE = ('_ddm_irref_den', '_ddm_idet', '_ddm_berk', '_sdm_rref_den', '_sdm_berk',
          '_sdm_berk_combine', '_ddm_fflu', '_ddm_fflu_p', '_ddm_fflu_l', '_ddm_fflu_d')


def disable(monkeypatch):
    for name in NATIVE:
        fallback = False if name.startswith('_ddm_fflu') else None
        monkeypatch.setattr(_core, name, lambda *args, _v=fallback, **kwargs: _v)


@pytest.mark.parametrize('name', sorted(DOMAINS))
@pytest.mark.parametrize('variant', sorted(VARIANTS))
def test_native_equals_complete_python_body(name, variant, monkeypatch):
    K, M = matrix(name, variant)
    ops = ['fflu', 'rref_den', 'srref_den', 'rref_den_auto']
    if M.shape[0] == M.shape[1]:
        ops += ['det', 'ddet', 'berk', 'sberk']
    native = [compute(K, M, op) for op in ops]
    disable(monkeypatch)
    assert [compute(K, M, op) for op in ops] == native


def rows(K, data):
    return DDM([[K.from_sympy(sx.sympify(e)) for e in row] for row in data],
               (len(data), len(data[0])), K)


def aliases(a, den):
    return [[e is den for e in row] for row in a]


@pytest.mark.parametrize('name', sorted(DOMAINS))
@pytest.mark.parametrize('variant', ['base', 'sing', 'zcol'])
def test_kernels_are_native_and_mutate_like_python(name, variant, monkeypatch):
    K, data = DOMAINS[name]
    data = VARIANTS[variant](data)
    stock = dense._ddm_stock_type()
    a = rows(K, data)
    det = _core._ddm_idet(a, K, stock)
    b = rows(K, data)
    den, pivots = _core._ddm_irref_den(b, K, stock)
    poly = _core._ddm_berk(rows(K, data), K, stock)
    assert det is not None and poly is not None and poly[0] == [K.one]
    sdm = DomainMatrix(rows(K, data).to_list(), (3, 3), K).to_sparse().rep
    assert _core._sdm_berk(dict(sdm), 3, K) == {i: c[0] for i, c in enumerate(poly) if c[0]}
    sparse = _core._sdm_rref_den({i: dict(r) for i, r in sdm.items()}, K)
    assert sparse is not None
    disable(monkeypatch)
    a2 = rows(K, data)
    assert dense.ddm_idet(a2, K) == det and a2.to_list() == a.to_list()
    b2 = rows(K, data)
    den2, pivots2 = dense.ddm_irref_den(b2, K)
    assert (den2, pivots2) == (den, pivots) and b2.to_list() == b.to_list()
    assert aliases(b2, den2) == aliases(b, den)
    assert dense.ddm_berk(rows(K, data), K) == poly
    assert sdm.charpoly() == [c[0] for c in poly]
    R, d, p = sdm.rref_den()
    assert ({i: dict(r) for i, r in R.items()}, d, list(p)) == (sparse[0], sparse[1], list(sparse[2]))


@pytest.mark.parametrize('name', sorted(DOMAINS))
def test_fflu_certificate(name):
    K, M = matrix(name, 'sing')
    P, L, D, U = M.fflu()
    # P*A = L*D**-1*U, multiplied through by prod(D) to stay in the ring.
    d = [D[i, i].element for i in range(D.shape[0])]
    assert all(d)
    total = K.one
    for di in d:
        total *= di
    adjugate = DomainMatrix.diag([K.exquo(total, di) for di in d], K)
    assert L*adjugate*U == (P*M)*total


def test_foreign_and_mixed_cells_decline():
    K = ZZ_I
    stock = dense._ddm_stock_type()
    a = rows(K, [[1, 2], [3, 4]])
    a[0][0] = GF(7)(1)
    assert _core._ddm_idet(a, K, stock) is None
    assert _core._ddm_irref_den(a, K, stock) is None
    assert _core._ddm_berk(a, K, stock) is None
    a = rows(K, [[1, 2], [3, 4]])
    a[1][1] = 4
    assert _core._ddm_idet(a, K, stock) is None
    # Inexact and expression domains keep their own Python bodies.
    for D in (sx.RR, sx.EX):
        a = rows(D, [[1, 2], [3, 4]])
        assert _core._ddm_idet(a, D, stock) is None
        assert _core._ddm_irref_den(a, D, stock) is None


def test_large_prime_field_and_wide_gaussian_determinant():
    p = 2**61 - 1
    K = GF(p)
    M = DomainMatrix([[K(3**(i*j+1) % 1000003) for j in range(6)] for i in range(6)], (6, 6), K)
    via_charpoly = M.charpoly()[-1]*(-1)**6
    assert M.det() == via_charpoly
    G = DomainMatrix([[ZZ_I.from_sympy((i+1)**j + (j-i)*I) for j in range(5)] for i in range(5)], (5, 5), ZZ_I)
    assert G.det() == G.charpoly()[-1]*(-1)**5


def test_composite_modulus_keeps_python_error_and_partial_mutation():
    # Z/6Z has zero divisors: the Python body raises NotInvertible after its
    # first writes, exactly as SymPy 1.14; the native kernels decline.
    K = GF(6)
    stock = dense._ddm_stock_type()
    a = [[K(2), K(1)], [K(3), K(4)]]
    assert _core._ddm_irref_den(a, K, stock) is None
    assert _core._ddm_idet(rows(K, [[2, 1], [3, 4]]), K, stock) is None
    with pytest.raises(Exception, match='zero divisor'):
        dense.ddm_irref_den(a, K)
    assert a == [[K(2), K(1)], [K(0), K(5)]]
