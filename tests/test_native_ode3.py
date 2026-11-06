"""Native nth_algebraic, factorable and 2nd_linear_bessel solvers
(symox::ode::algebraic, symox::ode::bessel) against SymPy 1.14 and against
the reference bodies."""
import pytest

import symox as S
from symox import _core
from symox.solvers.ode import _classify_native as cn
from symox.solvers.ode import single as s
from symox.solvers.ode.single import SingleODEProblem

x, y = S.symbols('x y')
f = S.Function('f')
fx = f(x)
d1 = fx.diff(x)
d2 = fx.diff(x, 2)
d3 = fx.diff(x, 3)
ENV = {'x': x, 'y': y, 'fx': fx, 'd1': d1, 'd2': d2, 'd3': d3, 'Float': S.Float, 'exp': S.exp, 'sin': S.sin}


def _eq(src):
    return eval(src, dict(ENV))


# SymPy 1.14: (equation, hint, classify_ode(eq, f(x)), dsolve(eq, f(x), hint=hint)).
CASES = [
    ('fx*(d1**2 - 1)', 'nth_algebraic', "('factorable', 'nth_algebraic', 'lie_group', 'nth_algebraic_Integral')", '[Eq(f(x), 0), Eq(f(x), C1 - x), Eq(f(x), C1 + x)]'),
    ('d2 - d1/x', 'nth_algebraic', "('factorable', 'nth_linear_euler_eq_homogeneous', 'Liouville', 'nth_order_reducible', '2nd_power_series_regular', 'Liouville_Integral')", 'ERR ValueError'),
    ('d1*d2', 'nth_algebraic', "('factorable', 'nth_algebraic', 'nth_order_reducible', 'nth_algebraic_Integral')", 'Eq(f(x), C1 + C2*x)'),
    ('d3 + d2 - d1', None, "('factorable', 'nth_linear_constant_coeff_homogeneous', 'nth_order_reducible')", None),
    ('(d1 - x)*(d2 + 1)', 'nth_algebraic', "('factorable', 'nth_algebraic', 'nth_order_reducible', 'nth_algebraic_Integral')", '[Eq(f(x), C1 + C2*x - x**2/2), Eq(f(x), C1 + x**2/2)]'),
    ('d2*(fx - 1)', 'nth_algebraic', "('factorable', 'nth_algebraic', 'Liouville', 'nth_algebraic_Integral', 'Liouville_Integral')", 'Eq(f(x), C1 + C2*x)'),
    ('(fx**2 - 4)*(d1 + fx)', 'factorable', "('factorable', 'separable', '1st_exact', '1st_power_series', 'lie_group', 'separable_Integral', '1st_exact_Integral')", '[Eq(f(x), 2), Eq(f(x), -2), Eq(f(x), C1*exp(-x))]'),
    ('fx*(d1 - x)', 'factorable', "('factorable', 'nth_algebraic', 'separable', '1st_exact', '1st_linear', 'Bernoulli', '1st_power_series', 'lie_group', 'nth_algebraic_Integral', 'separable_Integral', '1st_exact_Integral', '1st_linear_Integral', 'Bernoulli_Integral')", '[Eq(f(x), C1 + x**2/2), Eq(f(x), 0)]'),
    ('d1**2 - 1', 'factorable', "('factorable', 'nth_algebraic', 'lie_group', 'nth_algebraic_Integral')", '[Eq(f(x), C1 - x), Eq(f(x), C1 + x)]'),
    ('(d1 - y)**2', 'factorable', "('factorable', 'nth_algebraic', 'separable', '1st_exact', '1st_linear', 'Bernoulli', '1st_power_series', 'lie_group', 'nth_linear_euler_eq_nonhomogeneous_variation_of_parameters', 'nth_algebraic_Integral', 'separable_Integral', '1st_exact_Integral', '1st_linear_Integral', 'Bernoulli_Integral', 'nth_linear_euler_eq_nonhomogeneous_variation_of_parameters_Integral')", 'ERR ValueError'),
    ('fx**2*d1', 'factorable', "('factorable', 'nth_algebraic', 'separable', '1st_exact', '1st_linear', 'Bernoulli', '1st_power_series', 'lie_group', 'nth_algebraic_Integral', 'separable_Integral', '1st_exact_Integral', '1st_linear_Integral', 'Bernoulli_Integral')", 'Eq(f(x), C1)'),
    ('x**2*d2 + x*d1 + (x**2 - 4)*fx', '2nd_linear_bessel', "('factorable', '2nd_linear_bessel', '2nd_power_series_regular')", 'Eq(f(x), C1*besselj(2, x) + C2*bessely(2, x))'),
    ('4*x**2*d2 + 8*x*d1 + (x - 3)*fx', '2nd_linear_bessel', "('factorable', '2nd_linear_bessel', '2nd_power_series_regular')", 'Eq(f(x), (C1*besselj(2, sqrt(x)) + C2*bessely(2, sqrt(x)))/sqrt(x))'),
    ('(x - 1)**2*d2 + (x - 1)*d1 + ((x - 1)**2 - 1)*fx', '2nd_linear_bessel', "('factorable', '2nd_linear_bessel', '2nd_power_series_ordinary')", 'Eq(f(x), C1*besselj(1, x - 1) + C2*bessely(1, x - 1))'),
    ('x**2*d2 + x*d1 + (y**2*x**2 - 1)*fx', '2nd_linear_bessel', "('factorable', '2nd_linear_bessel', '2nd_power_series_regular')", 'Eq(f(x), C1*besselj(1, x*y) + C2*bessely(1, x*y))'),
    ('x*d2 + d1 + x*fx', '2nd_linear_bessel', "('factorable', '2nd_linear_bessel', '2nd_power_series_regular')", 'Eq(f(x), C1*besselj(0, x) + C2*bessely(0, x))'),
    ('d2 + fx/x', '2nd_linear_bessel', "('2nd_linear_bessel', '2nd_power_series_regular')", 'Eq(f(x), sqrt(x)*(C1*besselj(1, 2*sqrt(x)) + C2*bessely(1, 2*sqrt(x))))'),
]


@pytest.mark.parametrize('src,hint,cls,sol', CASES)
def test_against_sympy(src, hint, cls, sol):
    eq = _eq(src)
    assert str(S.classify_ode(eq, fx)) == cls
    if hint is None:
        return
    try:
        got = str(S.dsolve(eq, fx, hint=hint))
    except Exception as e:  # noqa: BLE001
        got = 'ERR %s' % type(e).__name__
    assert got == sol


def _problem(src):
    eq = _eq(src)
    return SingleODEProblem(eq, fx, x)


def _reference(solver_cls, src, attr):
    """The reference body with the native kernels hidden."""
    names = ('ode_nth_algebraic', 'ode_factorable_match', 'ode_factorable_solutions',
             'ode_bessel_match', 'ode_bessel_solution', 'ode_bessel_rejected')
    mp = pytest.MonkeyPatch()
    try:
        for n in names:
            mp.setattr(_core, n, None)
        solver = solver_cls(_problem(src))
        m = solver._matches()
        value = getattr(solver, attr, None) if m else None
        sols = solver._get_general_solution() if m else None
        return m, value, sols
    finally:
        mp.undo()


def _native(solver_cls, src, attr):
    solver = solver_cls(_problem(src))
    m = solver._matches()
    value = getattr(solver, attr, None) if m else None
    sols = solver._get_general_solution() if m else None
    return m, value, sols


ALGEBRAIC = [c[0] for c in CASES]


def _same(a, b):
    """Equal, the integration constants of nth_algebraic being fresh Dummies
    on each call (compared by their printed form)."""
    return a == b or str(a) == str(b) and a[0] == b[0]


@pytest.mark.parametrize('src', ALGEBRAIC)
def test_nth_algebraic_native_equals_reference(src):
    assert _same(_native(s.NthAlgebraic, src, 'solutions')[:2], _reference(s.NthAlgebraic, src, 'solutions')[:2])


@pytest.mark.parametrize('src', ALGEBRAIC)
def test_factorable_native_equals_reference(src):
    assert _native(s.Factorable, src, 'eqs')[:2] == _reference(s.Factorable, src, 'eqs')[:2]


@pytest.mark.parametrize('src', ALGEBRAIC)
def test_bessel_native_equals_reference(src):
    assert _native(s.SecondLinearBessel, src, 'rn') == _reference(s.SecondLinearBessel, src, 'rn')


def test_native_path_is_taken(monkeypatch):
    calls = []
    for name in ('ode_nth_algebraic', 'ode_factorable_match', 'ode_factorable_solutions',
                 'ode_bessel_match', 'ode_bessel_solution'):
        native = getattr(_core, name)
        monkeypatch.setattr(_core, name, lambda *a, _n=native, _name=name: (calls.append(_name), _n(*a))[1])
    # The reference bodies would call these; the native path must not.
    for mod, name in ((s, 'collect'), (s, 'fraction')):
        monkeypatch.setattr(mod, name, lambda *a, **k: pytest.fail('reference body ran'))
    assert str(S.dsolve(_eq('fx*(d1**2 - 1)'), fx, hint='nth_algebraic')) == '[Eq(f(x), 0), Eq(f(x), C1 - x), Eq(f(x), C1 + x)]'
    assert str(S.dsolve(_eq('(fx**2 - 4)*(d1 + fx)'), fx, hint='factorable')) == '[Eq(f(x), 2), Eq(f(x), -2), Eq(f(x), C1*exp(-x))]'
    assert str(S.dsolve(_eq('4*x**2*d2 + 8*x*d1 + (x - 3)*fx'), fx, hint='2nd_linear_bessel')) == \
        'Eq(f(x), (C1*besselj(2, sqrt(x)) + C2*bessely(2, sqrt(x)))/sqrt(x))'
    for name in ('ode_nth_algebraic', 'ode_factorable_match', 'ode_factorable_solutions',
                 'ode_bessel_match', 'ode_bessel_solution'):
        assert name in calls


def test_nth_algebraic_kernel_values():
    diffx = s.NthAlgebraic._get_diffx(x)
    out = _core.ode_nth_algebraic(fx*(d1**2 - 1), fx, x, diffx(x),
                                  lambda e, v: S.solve(e, v, simplify=False), S.simplify)
    assert [str(e) for e in out][0] == 'Eq(f(x), 0)'
    assert len(out) == 3
    # No solution: an empty list (the reference returns False).
    assert _core.ode_nth_algebraic(d1 - x, fx, x, diffx(x), lambda e, v: [], S.simplify) == []
    # NotImplementedError from solve counts as no solution.
    def ni(e, v):
        raise NotImplementedError
    assert _core.ode_nth_algebraic(d1 - x, fx, x, diffx(x), ni, S.simplify) == []
    # A non-list result of solve declines.
    assert _core.ode_nth_algebraic(d1 - x, fx, x, diffx(x), lambda e, v: {}, S.simplify) is None


def test_hook_errors_propagate():
    diffx = s.NthAlgebraic._get_diffx(x)
    def bad(e, v):
        raise ValueError('boom')
    with pytest.raises(ValueError, match='boom'):
        _core.ode_nth_algebraic(d1 - x, fx, x, diffx(x), bad, S.simplify)


def test_factorable_unsolved_reports_last_equation():
    def ni(eq, func):
        raise NotImplementedError
    sols, last = _core.ode_factorable_solutions([d1 - x, d1 - y], fx, ni)
    assert sols == [] and last == d1 - y
    sols, last = _core.ode_factorable_solutions([d1 - x], fx, S.dsolve)
    assert last is None and str(sols) == '[Eq(f(x), C1 + x**2/2)]'


def test_bessel_symbolic_c4_declines_solution():
    c = S.Symbol('c')
    assert _core.ode_bessel_solution((S.Integer(1), S.Integer(1), c, S.Integer(1), S.Integer(0)), fx, x,
                                     S.Symbol('C1'), S.Symbol('C2')) is None
    assert cn.bessel_solution({'n': S.Integer(1), 'a4': S.Integer(1), 'c4': S.Float(1.0), 'd4': S.Integer(1),
                               'b4': S.Integer(0)}, fx, x, S.Symbol('C1'), S.Symbol('C2')) is None
