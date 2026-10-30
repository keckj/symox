"""Native nth_order_reducible and Riccati_special_minus2 matches against
SymPy 1.14 and against the reference matchers."""
import pytest

import symox as S
from symox import _core
from symox.solvers.ode import _classify_native as cn
from symox.solvers.ode.single import NthOrderReducible, RiccatiSpecial, SingleODEProblem

x, a, b = S.symbols('x a b')
f = S.Function('f')
fx = f(x)
ENV = {'x': x, 'a': a, 'b': b, 'fx': fx, 'sin': S.sin, 'exp': S.exp, 'Float': S.Float}
ENV.update({'d%d' % n: fx.diff(x, n) for n in range(1, 5)})


def _eq(src):
    return eval(src, dict(ENV))


# SymPy 1.14: RiccatiSpecial(...).matches() and wilds_match().
RICCATI = [
    ('d1 - 2*fx**2 - 3*fx/x + 5/x**2', ['1', '-2', '-3', '5']),
    ('d1 + fx**2 + fx/x + 1/x**2', ['1', '1', '1', '1']),
    ('3*d1 - a*fx**2 + b/x**2', ['3', '-a', '0', 'b']),
    ('d1 - fx**2 - 2/x**2', ['1', '-1', '0', '-2']),
    ('x**2*d1 - x**2*fx**2 - 2', None),
    ('d1 - fx**2 - fx/x', ['1', '-1', '-1', '0']),
    ('d1 - fx**2', ['1', '-1', '0', '0']),
    ('a*d1 + b*fx**2 + fx/x', ['a', 'b', '1', '0']),
    ('d1 - x*fx**2 + 1/x**2', None),
    ('d1 + fx**2 - 1/x', None),
]

# SymPy 1.14: NthOrderReducible(...).matches() and smallest.
REDUCIBLE = [
    ('d2 + x*d1 - x', 1), ('d3 + d2', 2), ('d2 + fx', None), ('x*d2 - d1', 1), ('d3 - d1**2', 1),
    ('d2 + d1 + fx', None), ('d1 - x', None), ('d4 + d2**2 - x', 2), ('d2*d3 - 1', 2), ('d2 + sin(d1)', 1),
    ('d2 + exp(x)*d1**2', 1), ('d3 + d2 + fx', None),
]

# SymPy 1.14: dsolve(eq, f(x), hint=hint) and classify_ode(eq, f(x)).
DSOLVE = [
    ('d1 - 2*fx**2 - 3*fx/x + 5/x**2', 'Riccati_special_minus2',
     'Eq(f(x), -(-sqrt(14)*I*tan(C1 + sqrt(14)*I*log(x)) + 2)/(2*x))'),
    ('x*d2 - d1', 'nth_order_reducible', 'Eq(f(x), C1 + C2*x**2)'),
    ('d3 + d2', 'nth_order_reducible', 'Eq(f(x), C1 + C2*x + C3*exp(-x))'),
]
CLASSIFY = [
    ('d1 - 2*fx**2 - 3*fx/x + 5/x**2', ('1st_rational_riccati', 'Riccati_special_minus2', 'separable_reduced',
                                        'lie_group', 'separable_reduced_Integral')),
    ('x*d2 - d1', ('factorable', 'nth_linear_euler_eq_homogeneous', 'Liouville', 'nth_order_reducible',
                   '2nd_power_series_regular', 'Liouville_Integral')),
    ('d2 + x*d1 - x', ('factorable', 'nth_order_reducible')),
]


def _riccati(src):
    solver = RiccatiSpecial(SingleODEProblem(_eq(src), fx, x))
    return [str(v) for v in solver.wilds_match()] if solver.matches() else None


def _reducible(src):
    solver = NthOrderReducible(SingleODEProblem(_eq(src), fx, x))
    return solver.smallest if solver.matches() else None


@pytest.fixture
def stock(monkeypatch):
    def off():
        monkeypatch.setattr(cn, 'nth_order_reducible', lambda solver: None)
        monkeypatch.setattr(cn, 'solver_match', lambda solver: None)
    return off


def test_binding_present():
    assert hasattr(_core, 'ode_nth_order_reducible')


@pytest.mark.parametrize('src,want', RICCATI)
def test_riccati_special_matches_sympy(src, want, stock):
    assert _riccati(src) == want
    stock()
    assert _riccati(src) == want


@pytest.mark.parametrize('src,want', RICCATI)
def test_riccati_special_decided_natively(src, want):
    solver = RiccatiSpecial(SingleODEProblem(_eq(src), fx, x))
    value = cn.solver_match(solver)
    if want is None:
        assert value is False
    else:
        assert [str(value[w]) for w in solver.wilds()] == want


@pytest.mark.parametrize('src,want', REDUCIBLE)
def test_nth_order_reducible_matches_sympy(src, want, stock):
    assert _reducible(src) == want
    solver = NthOrderReducible(SingleODEProblem(_eq(src), fx, x))
    assert cn.nth_order_reducible(solver) == ((True, want) if want is not None else (False, None))
    stock()
    assert _reducible(src) == want


@pytest.mark.parametrize('src,hint,want', DSOLVE)
def test_dsolve_matches_sympy(src, hint, want):
    assert str(S.dsolve(_eq(src), fx, hint=hint)) == want


@pytest.mark.parametrize('src,want', CLASSIFY)
def test_classify_ode_matches_sympy(src, want):
    assert S.classify_ode(_eq(src), fx) == want


def test_unmodelled_declines():
    g = S.Function('g')
    solver = NthOrderReducible(SingleODEProblem(fx.diff(x, 2) + S.Integral(fx.diff(x), x), fx, x))
    assert cn.nth_order_reducible(solver) is None
    assert cn.nth_order_reducible(NthOrderReducible(SingleODEProblem(fx.diff(x, 2) + g(x)*fx.diff(x), fx, x))) == (True, 1)
