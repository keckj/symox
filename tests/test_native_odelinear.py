"""Differential tests of the native nth_linear_* detections
(symox::ode::linear) against SymPy 1.14 and the reference matchers.

EXPECTED holds SymPy 1.14 classify_ode(eq, f(x)) results; DSOLVE holds
SymPy 1.14 dsolve(eq, f(x), hint=hint) results."""
import pytest

import symox as S
from symox import _core
from symox.solvers.ode import _classify_native as cn
from symox.solvers.ode import single as s
from symox.solvers.ode.single import SingleODEProblem

x, y, a, b = S.symbols('x y a b')
f = S.Function('f')
fx = f(x)
d1 = fx.diff(x)
d2 = fx.diff(x, 2)
d3 = fx.diff(x, 3)
NS = {name: getattr(S, name) for name in dir(S) if not name.startswith('_')}
NS.update(x=x, y=y, a=a, b=b, f=f, fx=fx, d1=d1, d2=d2, d3=d3)
CLASSES = [s.NthLinearConstantCoeffHomogeneous, s.NthLinearConstantCoeffVariationOfParameters,
           s.NthLinearConstantCoeffUndeterminedCoefficients, s.NthLinearEulerEqHomogeneous,
           s.NthLinearEulerEqNonhomogeneousVariationOfParameters,
           s.NthLinearEulerEqNonhomogeneousUndeterminedCoefficients]
PREPROCESSED = (s.NthLinearEulerEqHomogeneous, s.NthLinearEulerEqNonhomogeneousVariationOfParameters)

EXPECTED = [
    ('d2 + fx', ('nth_linear_constant_coeff_homogeneous', '2nd_power_series_ordinary', '2nd_nonlinear_autonomous_conserved', '2nd_nonlinear_autonomous_conserved_Integral')),
    ('d2 - 3*d1 + 2*fx - exp(x)', ('factorable', 'nth_linear_constant_coeff_undetermined_coefficients', 'nth_linear_constant_coeff_variation_of_parameters', 'nth_linear_constant_coeff_variation_of_parameters_Integral')),
    ('d2 + fx - sin(x)*cos(x)', ('nth_linear_constant_coeff_variation_of_parameters', 'nth_linear_constant_coeff_variation_of_parameters_Integral')),
    ('d2 + fx - log(x)', ('nth_linear_constant_coeff_variation_of_parameters', 'nth_linear_constant_coeff_variation_of_parameters_Integral')),
    ('d3 - d1 - x**2', ('factorable', 'nth_linear_constant_coeff_undetermined_coefficients', 'nth_linear_constant_coeff_variation_of_parameters', 'nth_order_reducible', 'nth_linear_constant_coeff_variation_of_parameters_Integral')),
    ('a*d2 + b*fx - exp(y*x)', ('nth_linear_constant_coeff_undetermined_coefficients', 'nth_linear_constant_coeff_variation_of_parameters', 'nth_linear_constant_coeff_variation_of_parameters_Integral')),
    ('d2 + x*fx', ('2nd_linear_airy', '2nd_power_series_ordinary')),
    ('x**2*d2 + x*d1 + fx', ('factorable', 'nth_linear_euler_eq_homogeneous', '2nd_linear_bessel', '2nd_power_series_regular')),
    ('x**2*d2 + x*d1 + fx - log(x)', ('factorable', 'nth_linear_euler_eq_nonhomogeneous_undetermined_coefficients', 'nth_linear_euler_eq_nonhomogeneous_variation_of_parameters', 'nth_linear_euler_eq_nonhomogeneous_variation_of_parameters_Integral')),
    ('x**2*d2 - 2*x*d1 + 2*fx - x**3', ('factorable', 'nth_linear_euler_eq_nonhomogeneous_undetermined_coefficients', 'nth_linear_euler_eq_nonhomogeneous_variation_of_parameters', 'nth_linear_euler_eq_nonhomogeneous_variation_of_parameters_Integral')),
    ('x**2*d2 + x*d1 + fx - 1/log(x)', ('factorable', 'nth_linear_euler_eq_nonhomogeneous_variation_of_parameters', 'nth_linear_euler_eq_nonhomogeneous_variation_of_parameters_Integral')),
    ('x**2*d2 + fx**2', ()),
    ('d2 + fx - exp(x)*exp(2*x)', ('nth_linear_constant_coeff_undetermined_coefficients', 'nth_linear_constant_coeff_variation_of_parameters', 'nth_linear_constant_coeff_variation_of_parameters_Integral')),
    ('d2 + fx - sin(x**2)', ('nth_linear_constant_coeff_variation_of_parameters', 'nth_linear_constant_coeff_variation_of_parameters_Integral')),
    ('d2 + fx - 2**x', ('nth_linear_constant_coeff_undetermined_coefficients', 'nth_linear_constant_coeff_variation_of_parameters', 'nth_linear_constant_coeff_variation_of_parameters_Integral')),
    ('d2 - fx - x*exp(x)*sin(x)', ('nth_linear_constant_coeff_undetermined_coefficients', 'nth_linear_constant_coeff_variation_of_parameters', 'nth_linear_constant_coeff_variation_of_parameters_Integral')),
    ('d2 + 2*d1 + fx - exp(-x)', ('factorable', 'nth_linear_constant_coeff_undetermined_coefficients', 'nth_linear_constant_coeff_variation_of_parameters', 'nth_linear_constant_coeff_variation_of_parameters_Integral')),
    ('x**2*d2 + x*d1 + fx - sin(log(x))', ('factorable', 'nth_linear_euler_eq_nonhomogeneous_undetermined_coefficients', 'nth_linear_euler_eq_nonhomogeneous_variation_of_parameters', 'nth_linear_euler_eq_nonhomogeneous_variation_of_parameters_Integral')),
    ('d2 + fx - tan(x)', ('nth_linear_constant_coeff_variation_of_parameters', 'nth_linear_constant_coeff_variation_of_parameters_Integral')),
    ('x**3*d3 + x*d1 - fx - x*log(x)', ('factorable', 'nth_linear_euler_eq_nonhomogeneous_undetermined_coefficients', 'nth_linear_euler_eq_nonhomogeneous_variation_of_parameters', 'nth_linear_euler_eq_nonhomogeneous_variation_of_parameters_Integral')),
    ('d2 + fx - sin(y)*cos(y)*exp(x)', ('nth_linear_constant_coeff_variation_of_parameters', 'nth_linear_constant_coeff_variation_of_parameters_Integral')),
    ('d2 + fx - exp(x)/(x + 1)', ('nth_linear_constant_coeff_variation_of_parameters', 'nth_linear_constant_coeff_variation_of_parameters_Integral')),
    ('y*x**2*d2 + y*x*d1 + fx - x**2', ('factorable', 'nth_linear_euler_eq_nonhomogeneous_undetermined_coefficients', 'nth_linear_euler_eq_nonhomogeneous_variation_of_parameters', 'nth_linear_euler_eq_nonhomogeneous_variation_of_parameters_Integral')),
]

DSOLVE = [
    ('d2 - 3*d1 + 2*fx - exp(x)', 'nth_linear_constant_coeff_undetermined_coefficients',
     'Eq(f(x), (C1 + C2*exp(x) - x)*exp(x))'),
    ('d2 - fx - x*exp(x)*sin(x)', 'nth_linear_constant_coeff_undetermined_coefficients',
     'Eq(f(x), C2*exp(-x) + (C1 - x*sin(x)/5 - 2*x*cos(x)/5 + 14*sin(x)/25 - 2*cos(x)/25)*exp(x))'),
    ('x**2*d2 + x*d1 + fx - log(x)', 'nth_linear_euler_eq_nonhomogeneous_undetermined_coefficients',
     'Eq(f(x), C1*sin(log(x)) + C2*cos(log(x)) + log(x))'),
    ('d2 + fx - log(x)', 'nth_linear_constant_coeff_variation_of_parameters',
     'Eq(f(x), (C1 - Ci(x))*cos(x) + (C2 - Si(x))*sin(x) + log(x))'),
    ('x**2*d2 - 2*x*d1 + 2*fx', 'nth_linear_euler_eq_homogeneous', 'Eq(f(x), x*(C1 + C2*x))'),
]

native = pytest.mark.skipif(not hasattr(_core, 'ode_nth_linear_match'), reason='native kernel not built')


def _eq(src):
    return eval(src, NS)


@pytest.mark.parametrize('src, hints', EXPECTED)
def test_classify_ode_matches_sympy(src, hints):
    assert S.classify_ode(_eq(src), fx) == hints


@pytest.mark.parametrize('src, hint, sol', DSOLVE)
def test_dsolve_matches_sympy(src, hint, sol):
    assert str(S.dsolve(_eq(src), fx, hint=hint)) == sol


def _reference(cls, eq, monkeypatch):
    with monkeypatch.context() as m:
        m.setattr(cn, 'nth_linear_match', lambda solver, eq: None)
        m.setattr(cn, 'linear_coefficients', lambda eq, func, order: cn._MISSING)
        solver = cls(SingleODEProblem(eq, fx, x))
        ok = solver._matches()
        trial = set(solver.trialset) if ok and cls is s.NthLinearConstantCoeffUndeterminedCoefficients else None
        return bool(ok), (list(solver.r.items()) if ok else None), trial


@native
@pytest.mark.parametrize('src', [e[0] for e in EXPECTED])
@pytest.mark.parametrize('cls', CLASSES)
def test_native_matcher_equals_reference(src, cls, monkeypatch):
    eq = _eq(src)
    solver = cls(SingleODEProblem(eq, fx, x))
    problem = solver.ode_problem
    source = problem.eq_preprocessed if cls in PREPROCESSED else problem.eq_high_order_free
    value = cn.nth_linear_match(solver, source)
    if value is None:
        pytest.skip('declined')
    ok = solver._matches()
    trial = set(solver.trialset) if ok and cls is s.NthLinearConstantCoeffUndeterminedCoefficients else None
    got = bool(ok), (list(solver.r.items()) if ok else None), trial
    assert got == _reference(cls, eq, monkeypatch)


@native
@pytest.mark.parametrize('src', [e[0] for e in EXPECTED])
def test_native_linear_coefficients_equal_reference(src, monkeypatch):
    eq = _eq(src)
    problem = SingleODEProblem(eq, fx, x)
    for source in (problem.eq_preprocessed, problem.eq_high_order_free):
        value = cn.linear_coefficients(source, fx, problem.order)
        assert value is not cn._MISSING
        with monkeypatch.context() as m:
            m.setattr(cn, 'linear_coefficients', lambda eq, func, order: cn._MISSING)
            ref = problem.get_linear_coefficients(source, fx, problem.order)
        assert (value is None and ref is None) or list(value.items()) == list(ref.items())


@native
def test_native_decides_the_textbook_cases():
    decided = {
        'd2 - 3*d1 + 2*fx - exp(x)': (False, True, True, False, False, False),
        'd2 + fx - log(x)': (False, True, False, False, False, False),
        'd2 + fx - sin(y)*cos(y)*exp(x)': (False, True, False, False, False, False),
        'x**2*d2 + x*d1 + fx - log(x)': (False, False, False, False, True, True),
        'x**2*d2 - 2*x*d1 + 2*fx': (False, False, False, True, False, False),
    }
    for src, want in decided.items():
        eq = _eq(src)
        got = []
        for cls in CLASSES:
            solver = cls(SingleODEProblem(eq, fx, x))
            problem = solver.ode_problem
            source = problem.eq_preprocessed if cls in PREPROCESSED else problem.eq_high_order_free
            value = cn.nth_linear_match(solver, source)
            assert value is not None, (src, cls.hint)
            got.append(value[0])
        assert tuple(got) == want, src


@native
def test_trialset_is_the_reference_set():
    eq = _eq('d2 - fx - x*exp(x)*sin(x)')
    solver = s.NthLinearConstantCoeffUndeterminedCoefficients(SingleODEProblem(eq, fx, x))
    assert solver._matches()
    assert '_trialset' not in solver.__dict__
    assert solver.trialset == {x*S.exp(x)*S.sin(x), x*S.exp(x)*S.cos(x), S.exp(x)*S.sin(x), S.exp(x)*S.cos(x)}


@native
def test_nonlinear_and_unmodelled_inputs():
    assert cn.linear_coefficients(d2 + fx**2, fx, 2) is None
    assert cn.linear_coefficients(d2 + S.Integral(fx, x), fx, 2) is cn._MISSING
    eq = d2 + fx - S.Piecewise((1, x > 0), (0, True))
    solver = s.NthLinearConstantCoeffUndeterminedCoefficients(SingleODEProblem(eq, fx, x))
    assert cn.nth_linear_match(solver, eq) is None
