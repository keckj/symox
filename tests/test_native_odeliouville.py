"""Native Liouville and 2nd_nonlinear_autonomous_conserved matches against
SymPy 1.14 and against the reference matchers."""
import pytest

import symox as S
from symox import _core
from symox.solvers.ode import _classify_native as cn
from symox.solvers.ode.single import Liouville, SecondNonlinearAutonomousConserved, SingleODEProblem

x, a = S.symbols('x a')
f, g = S.Function('f'), S.Function('g')
fx = f(x)
ENV = {'x': x, 'a': a, 'fx': fx, 'g': g, 'sin': S.sin, 'exp': S.exp, 'log': S.log, 'Float': S.Float, 'I': S.I}
ENV.update({'d%d' % n: fx.diff(x, n) for n in range(1, 4)})


def _eq(src):
    return eval(src, dict(ENV))


# SymPy 1.14: Liouville / SecondNonlinearAutonomousConserved ._matches() with
# wilds_match(), and classify_ode(eq, f(x)).
EXPECTED = [
    ('d2 + d1**2/fx + d1/x', ['1', '1/f(x)', '1/x'], None, ('factorable', 'Liouville', 'Liouville_Integral')),
    ('d2 + fx*d1**2 + x*d1', ['1', 'f(x)', 'x'], None, ('factorable', 'Liouville', 'Liouville_Integral')),
    ('x*d2 + d1**2', None, None, ('factorable', 'nth_order_reducible')),
    ('d2 + exp(fx)*d1**2', ['1', 'exp(f(x))', '0'], None, ('factorable', 'Liouville', 'Liouville_Integral')),
    ('d2 + d1**2 + d1', ['1', '1', '1'], None, ('factorable', 'Liouville', 'nth_order_reducible', 'Liouville_Integral')),
    ('2*d2 + 3*d1**2/fx', ['2', '3/f(x)', '0'], None, ('factorable', 'Liouville', 'Liouville_Integral')),
    ('d2 + x*d1**2', None, None, ('factorable', 'nth_order_reducible')),
    ('d2 + d1/x', ['1', '0', '1/x'], None, ('factorable', 'nth_linear_euler_eq_homogeneous', 'Liouville',
                                           'nth_order_reducible', '2nd_power_series_regular', 'Liouville_Integral')),
    ('d2 - fx', None, ['-f(x)'], ('nth_linear_constant_coeff_homogeneous', '2nd_power_series_ordinary',
                                  '2nd_nonlinear_autonomous_conserved', '2nd_nonlinear_autonomous_conserved_Integral')),
    ('d2 + d1**2 + 1', None, None, ('factorable', 'nth_order_reducible')),
    ('fx*d2 + d1**2 + a*d1', None, None, ('factorable',)),
    ('d2 + sin(fx) + fx', None, ['f(x) + sin(f(x))'], ('2nd_nonlinear_autonomous_conserved',
                                                       '2nd_nonlinear_autonomous_conserved_Integral')),
    ('d2 - exp(fx) + log(fx)', None, ['-exp(f(x)) + log(f(x))'], ('2nd_nonlinear_autonomous_conserved',
                                                                  '2nd_nonlinear_autonomous_conserved_Integral')),
    ('d2 + fx**3', None, ['f(x)**3'], ('2nd_nonlinear_autonomous_conserved', '2nd_nonlinear_autonomous_conserved_Integral')),
    ('d2 + x*fx', None, None, ('2nd_linear_airy', '2nd_power_series_ordinary')),
    ('d2 + d1 + fx', None, None, ('nth_linear_constant_coeff_homogeneous', '2nd_power_series_ordinary')),
    ('d2 + 1/fx**2', None, ['f(x)**(-2)'], ('2nd_nonlinear_autonomous_conserved', '2nd_nonlinear_autonomous_conserved_Integral')),
    ('fx*d2 - 1', None, ['-1/f(x)'], ('2nd_nonlinear_autonomous_conserved', '2nd_nonlinear_autonomous_conserved_Integral')),
    ('d2**2 + fx', None, None, ()),
    ('x*d2 + fx', None, None, ('2nd_linear_bessel', '2nd_power_series_regular')),
    ('d2 + a*fx', None, ['a*f(x)'], ('nth_linear_constant_coeff_homogeneous', '2nd_power_series_ordinary',
                                     '2nd_nonlinear_autonomous_conserved', '2nd_nonlinear_autonomous_conserved_Integral')),
    ('d2 + g(fx)', None, ['g(f(x))'], ('2nd_nonlinear_autonomous_conserved', '2nd_nonlinear_autonomous_conserved_Integral')),
    ('d2 + d1*fx', None, None, ('factorable',)),
]

# SymPy 1.14: dsolve(eq, f(x), hint=hint).
DSOLVE = [
    ('d2 + d1**2/fx + d1/x', 'Liouville', '[Eq(f(x), -sqrt(C1 + C2*log(x))), Eq(f(x), sqrt(C1 + C2*log(x)))]'),
]

CLASSES = (Liouville, SecondNonlinearAutonomousConserved)


def _decide(cls, src):
    solver = cls(SingleODEProblem(_eq(src), fx, x))
    return [str(v) for v in solver.wilds_match()] if solver._matches() else None


@pytest.fixture
def stock(monkeypatch):
    def off():
        monkeypatch.setattr(cn, 'solver_match', lambda solver: None)
    return off


def test_binding_present():
    assert hasattr(_core, 'ode_second_order_patterns')


@pytest.mark.parametrize('src,liouville,conserved,hints', EXPECTED)
def test_matches_sympy(src, liouville, conserved, hints, stock):
    assert _decide(Liouville, src) == liouville
    assert _decide(SecondNonlinearAutonomousConserved, src) == conserved
    assert S.classify_ode(_eq(src), fx) == hints
    stock()
    assert _decide(Liouville, src) == liouville
    assert _decide(SecondNonlinearAutonomousConserved, src) == conserved


@pytest.mark.parametrize('src', [row[0] for row in EXPECTED])
def test_decided_natively_equals_reference(src, monkeypatch):
    for cls in CLASSES:
        native = cn.solver_match(cls(SingleODEProblem(_eq(src), fx, x)))
        if src == 'd2**2 + fx':
            # f'' only inside a power: the reference divides by `eq.coeff(f'') == 0`
            continue
        assert native is not None
        solver = cls(SingleODEProblem(_eq(src), fx, x))
        monkeypatch.setattr(cn, 'solver_match', lambda solver: None)
        solver._matches()
        monkeypatch.undo()
        reference = getattr(solver, '_wilds_match', None)
        assert (native is False and reference is None) or dict(native) == dict(reference)


@pytest.mark.parametrize('src,hint,want', DSOLVE)
def test_dsolve_matches_sympy(src, hint, want):
    assert str(S.dsolve(_eq(src), fx, hint=hint)) == want


@pytest.mark.parametrize('src', ['d2 + Float(1.5)*d1**2', 'd2 + I*fx', 'd2 + S.Integral(fx, x)'])
def test_declines(src):
    env = dict(ENV, S=S)
    solver = Liouville(SingleODEProblem(eval(src, env), fx, x))
    assert cn.solver_match(solver) is None
