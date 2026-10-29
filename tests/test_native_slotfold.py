"""copyreg slot caches keep native admission; native Piecewise pairs and
piecewise_fold as SymPy 1.14 (expected values hard-coded from SymPy 1.14)."""
import copy
import copyreg
import importlib
import pickle
import sys

import pytest
import symox as sx

risch = importlib.import_module('symox.integrals.risch')
polymatrix = importlib.import_module('symox.polys.polymatrix')
exceptions = importlib.import_module('symox.polys.matrices.exceptions')
dispatch = importlib.import_module('symox._native_dispatch')
from symox.polys.matrices.ddm import DDM


def _symox_classes():
    found = {}
    for name, module in list(sys.modules.items()):
        if module is None or not (name == 'symox' or name.startswith('symox.')):
            continue
        for value in list(vars(module).values()):
            if isinstance(value, type):
                found[id(value)] = value
    return list(found.values())


@pytest.fixture
def slot_caches():
    """Yield a recorder of copied objects; afterwards remove every copyreg
    cache the test created, leaving the classes as an untouched import."""
    watched = {}

    def record(*classes):
        for cls in classes:
            for base in cls.__mro__:
                watched.setdefault(id(base), (base, '__slotnames__' in vars(base)))

    for cls in _symox_classes():
        record(cls)
    yield record
    for cls in _symox_classes():
        record(cls)
    for cls, had in watched.values():
        if not had and '__slotnames__' in vars(cls):
            try:
                delattr(cls, '__slotnames__')
            except (AttributeError, TypeError):
                pass


def _copy_all(record, *objects):
    for value in objects:
        record(type(value))
        copy.copy(value)


# -- risch polynomial_reduce ------------------------------------------------

def _risch_case():
    x, t = sx.symbols('slot_x slot_t')
    p = sx.Poly(3*t**4 - 2*t**2 + 5, t, domain='QQ')
    de = risch.DifferentialExtension(extension={'D': [sx.Poly(1, x), sx.Poly(t**2 + 1, t)]})
    return p, de


def _spy(monkeypatch, name):
    calls = []
    native = getattr(sx._core, name)

    def spy(*args):
        result = native(*args)
        calls.append(result is not None)
        return result
    monkeypatch.setattr(sx._core, name, spy)
    return calls


def test_polynomial_reduce_stays_native_after_copy_of_extension_and_poly(slot_caches, monkeypatch):
    p, de = _risch_case()
    expected = risch.polynomial_reduce(p, de)
    _copy_all(slot_caches, de, p)
    assert '__slotnames__' in vars(risch.DifferentialExtension)
    calls = _spy(monkeypatch, '_risch_polynomial_reduce')
    q, r = risch.polynomial_reduce(p, de)
    assert calls == [True]
    assert (q.as_expr(), r.as_expr()) == tuple(e.as_expr() for e in expected)
    assert p.as_dict() == (risch.derivation(q, de) + r).as_dict()


def test_polynomial_reduce_stays_native_after_pickle_of_poly(slot_caches, monkeypatch):
    p, de = _risch_case()
    expected = risch.polynomial_reduce(p, de)
    slot_caches(type(p))
    assert pickle.loads(pickle.dumps(p)) == p
    calls = _spy(monkeypatch, '_risch_polynomial_reduce')
    q, r = risch.polynomial_reduce(p, de)
    assert calls == [True]
    assert (q.as_expr(), r.as_expr()) == tuple(e.as_expr() for e in expected)


def test_a_slotnames_entry_that_is_not_the_cache_still_declines(slot_caches, monkeypatch):
    p, de = _risch_case()
    expected = risch.polynomial_reduce(p, de)
    slot_caches(risch.DifferentialExtension)
    had = vars(risch.DifferentialExtension).get('__slotnames__')
    monkeypatch.setattr(risch.DifferentialExtension, '__slotnames__', ['bogus'], raising=False)
    calls = _spy(monkeypatch, '_risch_polynomial_reduce')
    q, r = risch.polynomial_reduce(p, de)
    assert calls == []
    assert (q.as_expr(), r.as_expr()) == tuple(e.as_expr() for e in expected)
    assert had is None or isinstance(had, list)


def test_cache_removed_after_bootstrap_snapshot_is_not_a_change(slot_caches, monkeypatch):
    p, de = _risch_case()
    copy.copy(de)
    slot_caches(risch.DifferentialExtension)
    monkeypatch.delattr(risch.DifferentialExtension, '__slotnames__')
    calls = _spy(monkeypatch, '_risch_polynomial_reduce')
    risch.polynomial_reduce(p, de)
    assert calls == [True]


# -- risch splitfactor --------------------------------------------------------

def test_splitfactor_stays_native_after_copy(slot_caches):
    x, t = sx.symbols('slot_x slot_t')
    de = risch.DifferentialExtension(extension={'D': [sx.Poly(1, x), sx.Poly(t, t)]})
    p = sx.Poly(t**3 + 2*t**2 + t, t, domain='QQ')
    expected = risch.splitfactor(p, de)
    _copy_all(slot_caches, de, p)
    kernel = sx._core._risch_splitfactor
    calls = []
    previous = sys.getprofile()

    def profile(frame, event, arg):
        if event == 'c_call' and arg is kernel:
            calls.append('native')
    sys.setprofile(profile)
    try:
        actual = risch.splitfactor(p, de)
    finally:
        sys.setprofile(previous)
    assert calls == ['native']
    assert [a.as_expr() for a in actual] == [e.as_expr() for e in expected]
    assert actual[0].as_expr() == t**2 + 2*t + 1 and actual[1].as_expr() == t


# -- PolyMatrix exact ground reduction ---------------------------------------

def test_polymatrix_ground_stays_native_after_copy(slot_caches, monkeypatch):
    x = sx.Symbol('x')
    a = polymatrix.PolyMatrix(2, 3, [sx.Rational(1, 2), 2, 3, 4, 5, sx.Rational(-7, 3)], x, ring=sx.QQ)
    expected = a.rref()
    _copy_all(slot_caches, *a.flat(), sx.QQ, x)
    calls = _spy(monkeypatch, '_polymatrix_ground')
    reduced, pivots = a.rref()
    assert calls and all(calls)
    assert pivots == expected[1] == (0, 1)
    assert [e.as_expr() for e in reduced.flat()] == [e.as_expr() for e in expected[0].flat()]
    assert [e.as_expr() for e in reduced.flat()] == [1, 0, sx.Rational(-118, 33), 0, 1, sx.Rational(79, 33)]


# -- DomainMatrix LU solve / LLL exception snapshots --------------------------

def _cache_exception_classes(record, *classes):
    for cls in classes:
        for base in cls.__mro__:
            record(base)
            if base.__module__ == exceptions.__name__:
                copyreg._slotnames(base)


def test_lu_solve_stays_native_with_cached_exception_slotnames(slot_caches, monkeypatch):
    _cache_exception_classes(slot_caches, exceptions.DMNonInvertibleMatrixError)
    assert '__slotnames__' in vars(exceptions.DMError)
    q = sx.QQ
    a = DDM([[q(2), q(1)], [q(1), q(3)]], (2, 2), q)
    b = DDM([[q(3)], [q(5)]], (2, 1), q)
    calls = _spy(monkeypatch, '_ddm_ilu_solve')
    out = a.lu_solve(b)
    assert calls == [True]
    assert [[sx.Rational(int(v.numerator), int(v.denominator)) for v in row] for row in out] == [
        [sx.Rational(4, 5)], [sx.Rational(7, 5)]]


def test_lll_stays_native_with_cached_exception_slotnames(slot_caches, monkeypatch):
    _cache_exception_classes(slot_caches, exceptions.DMRankError)
    z = sx.ZZ
    a = DDM([[z(1), z(0), z(0)], [z(4), z(2), z(15)], [z(0), z(0), z(3)]], (3, 3), z)
    assert '__slotnames__' in vars(exceptions.DMRankError)
    calls = _spy(monkeypatch, '_ddm_lll_native')
    delta = sx.QQ(3, 4)
    basis, transform = a.lll_transform(delta)
    if not calls:
        # the first call may warm the ABC negative cache the kernel needs
        assert a.lll_transform(delta) == (basis, transform)
    assert calls and all(calls)
    # sympy 1.14: DM([[1,0,0],[4,2,15],[0,0,3]], ZZ).lll_transform()
    assert [[int(v) for v in row] for row in basis] == [[1, 0, 0], [0, 2, 0], [0, 0, 3]]
    assert [[int(v) for v in row] for row in transform] == [[1, 0, 0], [-4, 1, -5], [0, 0, 1]]


def test_stock_slotnames_matches_copyreg_for_slotted_classes():
    class Base:
        __slots__ = ('a', '__hidden')

    class _Child(Base):
        __slots__ = 'b'

    class Plain:
        pass
    for cls in (Base, _Child, Plain):
        expected = copyreg._slotnames(cls)
        assert dispatch._stock_slotnames(cls) == expected
        assert dispatch.copyreg_slot_cache(cls, vars(cls)['__slotnames__'])
    assert not dispatch.copyreg_slot_cache(Base, ['a'])
    assert not dispatch.copyreg_slot_cache(Plain, ('x',))
    assert dispatch.strip_slot_cache(Plain, {'__slotnames__': ['x']}, ()) is None


# -- native Piecewise pairs ---------------------------------------------------

def test_native_piecewise_pairs_are_stock_exprcondpairs():
    from symox.functions.elementary.piecewise import ExprCondPair
    x = sx.Symbol('x')
    p = sx.Piecewise((x, x > 0), (sx.sin(x), sx.And(x < -1, sx.Symbol('y') > 0)), (0, True))
    args = p.args
    assert all(type(a) is ExprCondPair for a in args)
    assert args[-1].cond is sx.true
    assert args[0] == ExprCondPair(x, x > 0) and hash(args[0]) == hash(ExprCondPair(x, x > 0))
    assert [str(a) for a in args] == ['(x, x > 0)', '(sin(x), (y > 0) & (x < -1))', '(0, True)']
    assert p.func(*args) == p


def test_piecewise_condition_with_piecewise_keeps_the_full_constructor():
    x, y = sx.symbols('x y')
    p = sx.Piecewise((1, y > 0), (0, True)).subs(y, sx.Piecewise((x, x > 0), (-x, True)))
    # sympy 1.14
    assert str(p) == 'Piecewise((1, ITE(x > 0, x > 0, x < 0)), (0, True))'


def test_patched_pair_constructor_is_still_called():
    from symox.functions.elementary.piecewise import ExprCondPair
    x = sx.Symbol('x')
    p = sx.Piecewise((x, x > 0), (0, True))
    seen = []
    original = ExprCondPair.__dict__['__new__']

    def new(cls, expr, cond):
        seen.append(cond)
        return original.__func__(cls, expr, cond)
    ExprCondPair.__new__ = new
    try:
        p.args
    finally:
        ExprCondPair.__new__ = original
    assert seen == [x > 0, sx.true]
    seen.clear()
    p.args
    assert seen == []


# -- piecewise_fold and the Sum/Product/Integral constructors (sympy 1.14) ---

x, y, n, k = sx.symbols('x y n k')
Pw = sx.Piecewise
FOLD_CASES = [
    (lambda: sx.piecewise_fold(Pw((x, x > 0), (0, True))*n),
     'Piecewise((n*x, x > 0), (0, True))'),
    (lambda: sx.Sum(Pw((x, x > 0), (0, True))*n, (n, 1, 3)),
     'Sum(n*Piecewise((x, x > 0), (0, True)), (n, 1, 3))'),
    (lambda: sx.Sum(Pw((n, n > 2), (0, True))*x, (n, 1, k)),
     'Sum(Piecewise((n*x, n > 2), (0, True)), (n, 1, k))'),
    (lambda: sx.Integral(Pw((x, x < 1), (x**2, True))*y + 1, (x, 0, 2)),
     'Integral(Piecewise((x*y + 1, x < 1), (x**2*y + 1, True)), (x, 0, 2))'),
    (lambda: sx.Product(Pw((n, n > 1), (1, True))*Pw((2, n < 3), (n, True)), (n, 1, 4)),
     'Product(Piecewise((2*n, (n > 1) & (n < 3)), (n**2, n > 1), (2, n < 3), (n, True)), (n, 1, 4))'),
    (lambda: sx.piecewise_fold(Pw((x, x < 0), (1, True)) + Pw((x**2, x < 0), (2, True))),
     'Piecewise((x**2 + x, x < 0), (3, True))'),
    (lambda: sx.piecewise_fold(Pw((Pw((x, y > 0), (1, True)), x < 1), (0, True))),
     'Piecewise((x, (y > 0) & (x < 1)), (1, x < 1), (0, True))'),
    (lambda: sx.piecewise_fold(sx.sin(Pw((x, x > 0), (0, True)))),
     'Piecewise((sin(x), x > 0), (0, True))'),
    (lambda: sx.piecewise_fold(Pw((x, x > 0), (0, True))*Pw((1, x > 0), (2, True)), evaluate=False),
     'Piecewise((x, x > 0), (0, True))'),
    (lambda: sx.Integral(Pw((1, x > y), (0, True))*Pw((x, y > 0), (0, True)), (x, 0, 1)),
     'Integral(Piecewise((x, (x > y) & (y > 0)), (0, True)), (x, 0, 1))'),
    (lambda: sx.piecewise_fold(Pw((x, x > 0), (x, x < -1), (0, True))*2, evaluate=None),
     'Piecewise((2*x, (x > 0) | (x < -1)), (0, True))'),
    (lambda: sx.piecewise_fold(Pw((x, x > 0), (0, True)) + Pw((1, y > 0), (2, True))),
     'Piecewise((x + 1, (x > 0) & (y > 0)), (x + 2, x > 0), (1, y > 0), (2, True))'),
]


@pytest.mark.parametrize('case', range(len(FOLD_CASES)))
def test_piecewise_fold_and_constructors_match_sympy(case):
    build, expected = FOLD_CASES[case]
    assert str(build()) == expected


def test_piecewise_fold_args_are_exprcondpairs_of_folded_branches():
    from symox.functions.elementary.piecewise import ExprCondPair
    folded = sx.piecewise_fold(Pw((x, x > 0), (0, True))*n)
    assert [type(a) for a in folded.args] == [ExprCondPair, ExprCondPair]
    assert folded.args == (ExprCondPair(n*x, x > 0), ExprCondPair(0, sx.true))

