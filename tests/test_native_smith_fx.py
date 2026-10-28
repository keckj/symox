"""Differential and protocol tests for the native Smith decomposition over the
PIDs beyond ZZ: QQ, GF(p), QQ_I, QQ(x) and the univariate rings QQ[x].

Expected invariants, Smith forms and both transformation matrices are
hard-coded from SymPy 1.14.0 (serialised by value: QQ[x] by its coefficient
list, GF(p) by its residue in [0, p), QQ(x) by numerator/denominator with a
monic denominator). Every case is also certified by s*m*t == S with
invertible s, t, and random matrices are compared exactly (value, domain and
element type) with the complete Python reference body.
"""
import json
import random

import pytest

from symox import QQ, QQ_I, GF, ZZ, Poly, symbols, sympify, sqrt, fraction, cancel, _core
from symox.polys.matrices import DomainMatrix
from symox.polys.matrices import normalforms as nf
from symox.polys.matrices import _smith_native
from symox.polys.matrices.normalforms import (smith_normal_form, smith_normal_decomp,
                                               invariant_factors)

x = symbols('x')


def _domain(name):
    if name == 'QQ':
        return QQ
    if name.startswith('GF'):
        return GF(int(name[3:-1]))
    if name == 'QQ[x]':
        return QQ[x]
    if name == 'QQ_I':
        return QQ_I
    if name == 'QQ(x)':
        return QQ.frac_field(x)
    raise ValueError(name)


def _ser(K, e):
    if K.is_FiniteField:
        return int(e) % K.mod
    v = K.to_sympy(e)
    if K.is_QQ or K.is_GaussianField:
        return str(v)
    if K.is_PolynomialRing:
        return [str(c) for c in Poly(v, x, domain=QQ).all_coeffs()]
    n, d = fraction(cancel(v))
    pn, pd = Poly(n, x, domain=QQ), Poly(d, x, domain=QQ)
    lc = pd.LC()
    return [[str(c / lc) for c in pn.all_coeffs()], [str(c / lc) for c in pd.all_coeffs()]]


def _build(K, rows):
    return DomainMatrix([[K.from_sympy(sympify(v)) for v in r] for r in rows], (len(rows), len(rows[0])), K)


# (domain, entries, SymPy 1.14 invariants / Smith form / s / t / smith_normal_form)
CASES = [
    ["QQ",[[2,"3/2"],[4,3]],{"invs":["1","0"],"S":[["1","0"],["0","0"]],"s":[["1/2","0"],["-2","1"]],"t":[["1","-3/4"],["0","1"]],"snf":[["1","0"],["0","0"]]}],
    ["QQ",[[0,0],[0,5]],{"invs":["1","0"],"S":[["1","0"],["0","0"]],"s":[["0","1/5"],["1","0"]],"t":[["0","1"],["1","0"]],"snf":[["1","0"],["0","0"]]}],
    ["QQ",[[0,0],[0,0]],{"invs":["0","0"],"S":[["0","0"],["0","0"]],"s":[["0","1"],["1","0"]],"t":[["0","1"],["1","0"]],"snf":[["0","0"],["0","0"]]}],
    ["QQ",[["-1/3",2,0],[1,"7/5",4]],{"invs":["1","1"],"S":[["1","0","0"],["0","1","0"]],"s":[["-3","0"],["15/37","5/37"]],"t":[["1","6","-120/37"],["0","1","-20/37"],["0","0","1"]],"snf":[["1","0","0"],["0","1","0"]]}],
    ["QQ",[[0,1],[2,0],[3,3]],{"invs":["1","1"],"S":[["1","0"],["0","1"],["0","0"]],"s":[["0","1/2","0"],["1","0","0"],["-3","-3/2","1"]],"t":[["1","0"],["0","1"]],"snf":[["1","0"],["0","1"],["0","0"]]}],
    ["QQ",[[0,0,0],[0,0,7],[1,0,0]],{"invs":["1","1","0"],"S":[["1","0","0"],["0","1","0"],["0","0","0"]],"s":[["0","0","1"],["0","1/7","0"],["1","0","0"]],"t":[["1","0","0"],["0","0","1"],["0","1","0"]],"snf":[["1","0","0"],["0","1","0"],["0","0","0"]]}],
    ["QQ",[["-2"]],{"invs":["1"],"S":[["1"]],"s":[["-1/2"]],"t":[["1"]],"snf":[["1"]]}],
    ["QQ",[[1,2,3],[2,4,6],[1,1,1]],{"invs":["1","1","0"],"S":[["1","0","0"],["0","1","0"],["0","0","0"]],"s":[["1","0","0"],["1","0","-1"],["-2","1","0"]],"t":[["1","-2","1"],["0","1","-2"],["0","0","1"]],"snf":[["1","0","0"],["0","1","0"],["0","0","0"]]}],
    ["GF(7)",[[2,3,1],[4,6,5]],{"invs":[1,1],"S":[[1,0,0],[0,1,0]],"s":[[4,0],[4,5]],"t":[[1,3,2],[0,0,1],[0,1,0]],"snf":[[1,0,0],[0,1,0]]}],
    ["GF(2)",[[1,1],[1,1]],{"invs":[1,0],"S":[[1,0],[0,0]],"s":[[1,0],[1,1]],"t":[[1,1],[0,1]],"snf":[[1,0],[0,0]]}],
    ["GF(2)",[[0,1,1],[1,0,1],[1,1,0]],{"invs":[1,1,0],"S":[[1,0,0],[0,1,0],[0,0,0]],"s":[[0,1,0],[1,0,0],[1,1,1]],"t":[[1,0,1],[0,1,1],[0,0,1]],"snf":[[1,0,0],[0,1,0],[0,0,0]]}],
    ["GF(101)",[[50,-3],[17,99]],{"invs":[1,1],"S":[[1,0],[0,1]],"s":[[99,0],[56,67]],"t":[[1,95],[0,1]],"snf":[[1,0],[0,1]]}],
    ["GF(5)",[[0,0],[0,3]],{"invs":[1,0],"S":[[1,0],[0,0]],"s":[[0,2],[1,0]],"t":[[0,1],[1,0]],"snf":[[1,0],[0,0]]}],
    ["GF(3)",[[0,2],[1,0],[2,1]],{"invs":[1,1],"S":[[1,0],[0,1],[0,0]],"s":[[0,1,0],[2,0,0],[1,1,1]],"t":[[1,0],[0,1]],"snf":[[1,0],[0,1],[0,0]]}],
    ["GF(13)",[[0,0,0]],{"invs":[0],"S":[[0,0,0]],"s":[[1]],"t":[[0,0,1],[1,0,0],[0,1,0]],"snf":[[0,0,0]]}],
    ["QQ[x]",[["x**2-1","x+1"],["x-1","x**2+3*x"]],{"invs":[["1"],["1","3","-2","-3","1"]],"S":[[["1"],["0"]],[["0"],["1","3","-2","-3","1"]]],"s":[[["0"],["1"]],[["1"],["1/4","1","-1/2","-5/4"]]],"t":[[["-1/4","-1"],["1","3","0"]],[["1/4"],["-1","1"]]],"snf":[[["1"],["0"]],[["0"],["1","3","-2","-3","1"]]]}],
    ["QQ[x]",[["2*x","x**2+1",3],["x-1",0,"x**3"],[1,"x","x+2"]],{"invs":[["1"],["1"],["1","1","0","-4","4","-2"]],"S":[[["1"],["0"],["0"]],[["0"],["1"],["0"]],[["0"],["0"],["1","1","0","-4","4","-2"]]],"s":[[["1/2"],["-1"],["0"]],[["-1"],["1"],["1","1"]],[["-1/4","-1/2","-1/2","3/2","-1/2"],["1/4","1/2","1/2","-3/2","-1/2"],["1/4","3/4","1","-1","0","-1/2"]]],"t":[[["1"],["1/8","1/2","3/4","1/4","1/4"],["-1/2","-3/2","-1","0","0","-1"]],[["0"],["-1/4","-1/2","-5/4"],["1","1","3","-1"]],[["0"],["1/4"],["-1","1"]]],"snf":[[["1"],["0"],["0"]],[["0"],["1"],["0"]],[["0"],["0"],["1","1","0","-4","4","-2"]]]}],
    ["QQ[x]",[["x**2",0],[0,"x"]],{"invs":[["1","0"],["1","0","0"]],"S":[[["1","0"],["0"]],[["0"],["1","0","0"]]],"s":[[["0"],["1"]],[["-1"],["1","0"]]],"t":[[["1"],["-1"]],[["1"],["0"]]],"snf":[[["1","0"],["0"]],[["0"],["1","0","0"]]]}],
    ["QQ[x]",[["-x",0],[0,"-2*x**2"]],{"invs":[["1","0"],["2","0","0"]],"S":[[["1","0"],["0"]],[["0"],["2","0","0"]]],"s":[[["-1"],["0"]],[["0"],["-1"]]],"t":[[["1"],["0"]],[["0"],["1"]]],"snf":[[["1","0"],["0"]],[["0"],["2","0","0"]]]}],
    ["QQ[x]",[["x",0],[0,"x"]],{"invs":[["1","0"],["1","0"]],"S":[[["1","0"],["0"]],[["0"],["1","0"]]],"s":[[["1"],["0"]],[["0"],["1"]]],"t":[[["1"],["0"]],[["0"],["1"]]],"snf":[[["1","0"],["0"]],[["0"],["1","0"]]]}],
    ["QQ[x]",[[0,0],[0,"x-3"]],{"invs":[["1","-3"],["0"]],"S":[[["1","-3"],["0"]],[["0"],["0"]]],"s":[[["0"],["1"]],[["1"],["0"]]],"t":[[["0"],["1"]],[["1"],["0"]]],"snf":[[["1","-3"],["0"]],[["0"],["0"]]]}],
    ["QQ[x]",[[0,"-x"],["x+1",0]],{"invs":[["1"],["1","1","0"]],"S":[[["1"],["0"]],[["0"],["1","1","0"]]],"s":[[["-1"],["1"]],[["-1","-1"],["1","0"]]],"t":[[["1"],["-1","0"]],[["-1"],["1","1"]]],"snf":[[["1"],["0"]],[["0"],["1","1","0"]]]}],
    ["QQ[x]",[["-3*x**2+1","2*x"],["x","-x**3"]],{"invs":[["1"],["3","0","-1","-2","0","0"]],"S":[[["1"],["0"]],[["0"],["3","0","-1","-2","0","0"]]],"s":[[["1"],["3","0"]],[["-1","0"],["-3","0","1"]]],"t":[[["1"],["3","0","0","-2","0"]],[["0"],["1"]]],"snf":[[["1"],["0"]],[["0"],["3","0","-1","-2","0","0"]]]}],
    ["QQ[x]",[["x-1",0,0],[0,"x-1",0],[0,0,"x-2"]],{"invs":[["1"],["1","-1"],["1","-3","2"]],"S":[[["1"],["0"],["0"]],[["0"],["1","-1"],["0"]],[["0"],["0"],["1","-3","2"]]],"s":[[["0"],["1"],["1"]],[["-1"],["1","-1"],["1","-1"]],[["0"],["1","-2"],["1","-1"]]],"t":[[["1"],["-1"],["0"]],[["1"],["0"],["-1","2"]],[["-1"],["0"],["1","-1"]]],"snf":[[["1"],["0"],["0"]],[["0"],["1","-1"],["0"]],[["0"],["0"],["1","-3","2"]]]}],
    ["QQ[x]",[["x-2",-1,0],[0,"x-2",0],[0,0,"x-2"]],{"invs":[["1"],["1","-2"],["1","-4","4"]],"S":[[["1"],["0"],["0"]],[["0"],["1","-2"],["0"]],[["0"],["0"],["1","-4","4"]]],"s":[[["1"],["0"],["0"]],[["0"],["0"],["1"]],[["1","-2"],["1"],["1","-2"]]],"t":[[["0"],["-1"],["1"]],[["-1"],["-1","2"],["1","-2"]],[["0"],["1"],["0"]]],"snf":[[["1"],["0"],["0"]],[["0"],["1","-2"],["0"]],[["0"],["0"],["1","-4","4"]]]}],
    ["QQ[x]",[["x",1,0],[-1,"x",0],[0,0,"x+1"]],{"invs":[["1"],["1"],["1","1","1","1"]],"S":[[["1"],["0"],["0"]],[["0"],["1"],["0"]],[["0"],["0"],["1","1","1","1"]]],"s":[[["0"],["-1"],["0"]],[["1/2"],["1/2","0"],["1"]],[["1/2","0","-1/2"],["1/2","0","-1/2","0"],["1","0","1"]]],"t":[[["1"],["1","0"],["-1","-1","0"]],[["0"],["1"],["-1","-1"]],[["0"],["-1/2","1/2"],["1/2","0","1/2"]]],"snf":[[["1"],["0"],["0"]],[["0"],["1"],["0"]],[["0"],["0"],["1","1","1","1"]]]}],
    ["QQ[x]",[["x-1",-2,-3,-4],[0,"x-1",-2,-3],[0,0,"x-1",-2],[0,0,0,"x-1"]],{"invs":[["1"],["1"],["1"],["1","-4","6","-4","1"]],"S":[[["1"],["0"],["0"],["0"]],[["0"],["1"],["0"],["0"]],[["0"],["0"],["1"],["0"]],[["0"],["0"],["0"],["1","-4","6","-4","1"]]],"s":[[["1"],["0"],["0"],["0"]],[["1/2","-1/2"],["1"],["0"],["0"]],[["-3/16","13/16","-17/16","7/16"],["-3/8","5/4","-7/8"],["1"],["0"]],[["15/256","-1/2","463/256","-113/32","989/256","-71/32","133/256"],["15/128","-113/128","175/64","-277/64","435/128","-133/128"],["-5/16","21/16","-35/16","19/16"],["-1"]]],"t":[[["0"],["9/8"],["-5/4","-1"],["-4","-4","0"]],[["-1/2"],["3/4"],["-15/16","13/16"],["-3","2","1"]],[["0"],["3/8","-7/8"],["-5/8","11/8"],["-2","4","-2"]],[["0"],["0"],["-5/16","1","-19/16"],["-1","3","-3","1"]]],"snf":[[["1"],["0"],["0"],["0"]],[["0"],["1"],["0"],["0"]],[["0"],["0"],["1"],["0"]],[["0"],["0"],["0"],["1","-4","6","-4","1"]]]}],
    ["QQ[x]",[["1/2*x","3/4"],["-x**2","x/3"]],{"invs":[["1"],["11/12","0","0"]],"S":[[["1"],["0"]],[["0"],["11/12","0","0"]]],"s":[[["1"],["0"]],[["4/9","0"],["-1"]]],"t":[[["0"],["3/4"]],[["4/3"],["-1/2","0"]]],"snf":[[["1"],["0"]],[["0"],["11/12","0","0"]]]}],
    ["QQ[x]",[["x**2+x","x**3","x"],["x",1,"x-1"]],{"invs":[["1"],["1","0"]],"S":[[["1"],["0"],["0"]],[["0"],["1","0"],["0"]]],"s":[[["0"],["1"]],[["1"],["-1","0","0","0"]]],"t":[[["0"],["-1","0"],["-1","1","0","1"]],[["1"],["1"],["1","-1","-1"]],[["0"],["1","1"],["1","0","-1","-1"]]],"snf":[[["1"],["0"],["0"]],[["0"],["1","0"],["0"]]]}],
    ["QQ[x]",[["x"],["-x**2"],["x**3+1"]],{"invs":[["1"]],"S":[[["1"]],[["0"]],[["0"]]],"s":[[["-1","0","0"],["0"],["1"]],[["1","0"],["1"],["0"]],[["1","0","0","1"],["0"],["-1","0"]]],"t":[[["1"]]],"snf":[[["1"]],[["0"]],[["0"]]]}],
    ["QQ[x]",[["-x**4","x**2-4","x+2"]],{"invs":[["1"]],"S":[[["1"],["0"],["0"]]],"s":[[["1"]]],"t":[[["-1/16"],["1","0","-4"],["1/16","1/8"]],[["-1/16","0","-1/4"],["1","0","0","0","0"],["1/16","1/8","1/4","1/2"]],[["0"],["0"],["1"]]],"snf":[[["1"],["0"],["0"]]]}],
    ["QQ[x]",[[0,0],[0,0]],{"invs":[["0"],["0"]],"S":[[["0"],["0"]],[["0"],["0"]]],"s":[[["0"],["1"]],[["1"],["0"]]],"t":[[["0"],["1"]],[["1"],["0"]]],"snf":[[["0"],["0"]],[["0"],["0"]]]}],
    ["QQ[x]",[["-5"]],{"invs":[["5"]],"S":[[["5"]]],"s":[[["-1"]]],"t":[[["1"]]],"snf":[[["5"]]]}],
    ["QQ_I",[["1+I",2],["I","3-I"]],{"invs":["1","1"],"S":[["1","0"],["0","1"]],"s":[["1/2 - I/2","0"],["-I/4","1/4 + I/4"]],"t":[["1","-1 + I"],["0","1"]],"snf":[["1","0"],["0","1"]]}],
    ["QQ_I",[[0,"I"],["-1/2",0]],{"invs":["1","1"],"S":[["1","0"],["0","1"]],"s":[["0","-2"],["-I","0"]],"t":[["1","0"],["0","1"]],"snf":[["1","0"],["0","1"]]}],
    ["QQ_I",[["2*I",0,1],[1,"I",0]],{"invs":["1","1"],"S":[["1","0","0"],["0","1","0"]],"s":[["-I/2","0"],["1/2","-I"]],"t":[["1","0","I/2"],["0","1","-1/2"],["0","0","1"]],"snf":[["1","0","0"],["0","1","0"]]}],
    ["QQ(x)",[["1/x","x"],["x+1",1]],{"invs":[[["1"],["1"]],[["1"],["1"]]],"S":[[[["1"],["1"]],[["0"],["1"]]],[[["0"],["1"]],[["1"],["1"]]]],"s":[[[["1","0"],["1"]],[["0"],["1"]]],[[["1","1","0"],["1","1","0","-1"]],[["-1"],["1","1","0","-1"]]]],"t":[[[["1"],["1"]],[["-1","0","0"],["1"]]],[[["0"],["1"]],[["1"],["1"]]]],"snf":[[[["1"],["1"]],[["0"],["1"]]],[[["0"],["1"]],[["1"],["1"]]]]}],
    ["QQ(x)",[[0,"1/(x-1)"],["x",0]],{"invs":[[["1"],["1"]],[["1"],["1"]]],"S":[[[["1"],["1"]],[["0"],["1"]]],[[["0"],["1"]],[["1"],["1"]]]],"s":[[[["0"],["1"]],[["1"],["1","0"]]],[[["1","-1"],["1"]],[["0"],["1"]]]],"t":[[[["1"],["1"]],[["0"],["1"]]],[[["0"],["1"]],[["1"],["1"]]]],"snf":[[[["1"],["1"]],[["0"],["1"]]],[[["0"],["1"]],[["1"],["1"]]]]}]
]


@pytest.mark.parametrize('name, rows, expected', CASES)
def test_matches_sympy_and_certificate(name, rows, expected):
    K = _domain(name)
    m = _build(K, rows)
    assert _smith_native.try_smith(m.to_list(), K, m.shape, True) is not None
    assert _smith_native.try_smith(m.to_list(), K, m.shape, False) is not None
    invs = invariant_factors(m)
    assert [_ser(K, v) for v in invs] == expected['invs']
    smf, s, t = smith_normal_decomp(m)
    assert [[_ser(K, v) for v in r] for r in smf.to_list()] == expected['S']
    assert [[_ser(K, v) for v in r] for r in s.to_list()] == expected['s']
    assert [[_ser(K, v) for v in r] for r in t.to_list()] == expected['t']
    snf = smith_normal_form(m)
    assert [[_ser(K, v) for v in r] for r in snf.to_dense().to_list()] == expected['snf']
    assert s.domain == K and t.domain == K and smf.domain == K and snf.domain == K
    assert all(type(v) is _core.DomainElement for r in s.to_list() + t.to_list() for v in r)
    assert all(type(v) is _core.DomainElement for v in invs)
    assert s * m * t == smf
    for d in (s.det(), t.det()):
        # A unit: a non-zero field element, a non-zero constant of QQ[x].
        assert d != K.zero and (not K.is_PolynomialRing or K.to_sympy(d).is_number)


def _reference(rows, K, full):
    return nf._smith_normal_decomp([list(r) for r in rows], K, shape=(len(rows), len(rows[0])), full=full)


def _random_poly(rng, K):
    if rng.random() < 0.25:
        return K.zero
    terms = sum(rng.randint(-4, 4) * x**k / rng.randint(1, 3) for k in range(rng.randint(0, 3)))
    return K.from_sympy(sympify(terms))


@pytest.mark.parametrize('name', ['QQ', 'GF(2)', 'GF(7)', 'GF(65537)', 'QQ[x]', 'QQ_I', 'QQ(x)'])
def test_random_against_python_reference(name):
    # The native result equals the Python body; where the Python body raises
    # (GF(p) index quirk, QQ_I / QQ(x) zero pivot `1 / 0`) the kernel declines.
    K = _domain(name)
    rng = random.Random(11)
    used = 0
    for _ in range(25 if name not in ('QQ[x]', 'QQ(x)') else 15):
        r, c = rng.randint(1, 4), rng.randint(1, 4)
        if name == 'QQ[x]':
            rows = [[_random_poly(rng, K) for _ in range(c)] for _ in range(r)]
        elif name == 'QQ(x)':
            rows = [[K.from_sympy(sympify(rng.choice([0, rng.randint(-3, 3)])) / (x + rng.randint(0, 2))) if rng.random() < 0.5
                     else K.from_sympy(sympify(rng.choice([0, rng.randint(-3, 3) * x]))) for _ in range(c)] for _ in range(r)]
        elif name == 'QQ_I':
            rows = [[K.from_sympy(rng.choice([0, rng.randint(-5, 5) + rng.randint(-5, 5) * sympify('I')])) for _ in range(c)]
                    for _ in range(r)]
        else:
            rows = [[K.from_sympy(sympify(rng.choice([0, rng.randint(-9, 9)])) / rng.randint(1, 4) if K.is_QQ
                                  else rng.choice([0, rng.randint(-9, 9)])) for _ in range(c)] for _ in range(r)]
        m = DomainMatrix(rows, (r, c), K)
        for full in (True, False):
            native = _smith_native.try_smith(m.to_list(), K, m.shape, full)
            try:
                expected = _reference(m.to_list(), K, full)
            except ZeroDivisionError:
                assert native is None
                continue
            assert native is not None
            assert native == expected
            used += 1
    assert used >= 6


def test_characteristic_matrix_invariants():
    # x*I - A: the invariant factors are the similarity invariants of A.
    R = QQ[x]
    A = [[2, 1, 0, 0], [0, 2, 0, 0], [0, 0, 2, 0], [0, 0, 0, 3]]
    m = DomainMatrix([[R(x if i == j else 0) - R(A[i][j]) for j in range(4)] for i in range(4)], (4, 4), R)
    assert _smith_native.try_smith(m.to_list(), R, m.shape, False) is not None
    expected = [1, 1, x - 2, ((x - 2)**2 * (x - 3)).expand()]
    assert [R.to_sympy(v).expand() for v in invariant_factors(m)] == expected


def test_empty_shapes_over_fields_and_qqx():
    for K in (QQ, GF(5), QQ[x]):
        for shape in ((0, 0), (0, 3), (2, 0)):
            m = DomainMatrix([[]] * shape[0] if shape[1] == 0 else [], shape, K)
            assert invariant_factors(m) == ()
            smf, S, T = smith_normal_decomp(m)
            assert S.to_dense() == DomainMatrix.eye(shape[0], K).to_dense()
            assert T.to_dense() == DomainMatrix.eye(shape[1], K).to_dense()


def test_out_of_scope_domains_keep_reference_paths():
    # QQ<a>: SymPy 1.14 raises TypeError at 1/m[0][0]; the Python body runs.
    A = QQ.algebraic_field(sqrt(2))
    m = DomainMatrix([[A.from_sympy(sqrt(2)), A(1)], [A(2), A.from_sympy(sqrt(2))]], (2, 2), A)
    assert _smith_native.try_smith(m.to_list(), A, m.shape, True) is None
    # Non-PIDs raise the reference ValueError.
    for K in (ZZ[x], QQ[x, symbols('y')]):
        p = DomainMatrix([[K(1), K(2)], [K(3), K(4)]], (2, 2), K)
        with pytest.raises(ValueError, match='principal ideal domain'):
            invariant_factors(p)
    # Mixed-domain cells decline.
    R = QQ[x]
    rows = [[R(x), QQ(1)], [R(1), R(x)]]
    assert _smith_native.try_smith(rows, R, (2, 2), True) is None


def test_qqx_canonical_unit_patch_declines(monkeypatch):
    R = QQ[x]
    m = DomainMatrix([[R(-x), R(0)], [R(0), R(-2*x**2)]], (2, 2), R)
    native = _smith_native.try_smith(m.to_list(), R, m.shape, True)
    assert native is not None
    calls = []
    original = _core.Domain.__dict__['canonical_unit']

    def canonical_unit(self, a):
        calls.append(a)
        return original(self, a)
    monkeypatch.setattr(_core.Domain, 'canonical_unit', canonical_unit)
    assert _smith_native.try_smith(m.to_list(), R, m.shape, True) is None
    smf, S, T = smith_normal_decomp(m)
    assert calls
    assert (tuple(smf.to_list()[i][i] for i in range(2)), S.to_list(), T.to_list()) == native
    monkeypatch.undo()
    assert _smith_native.try_smith(m.to_list(), R, m.shape, True) is not None
