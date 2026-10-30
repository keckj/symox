"""Differential tests of the native first-order classify_ode matches
(symox::ode::classify) against SymPy 1.14.

The expected values are SymPy 1.14 results: classify_ode(eq, f(x)) and
solver.wilds_match() after solver.matches() for each first-order pattern
solver (None when the solver does not match)."""
import pytest

import symox as S
from symox import _core
from symox.solvers.ode import _classify_native
from symox.solvers.ode.single import (SingleODEProblem, Separable, FirstExact, FirstLinear, Bernoulli,
    AlmostLinear, HomogeneousCoeffBest, LinearCoefficients, SeparableReduced, RiccatiSpecial, RationalRiccati)

x, a, b = S.symbols('x a b')
f = S.Function('f')
fx = f(x)
df = fx.diff(x)
NS = {name: getattr(S, name) for name in dir(S) if not name.startswith('_')}
NS.update(x=x, a=a, b=b, f=f, fx=fx, df=df)
CLASSES = [Separable, FirstExact, FirstLinear, Bernoulli, AlmostLinear, HomogeneousCoeffBest,
           LinearCoefficients, SeparableReduced, RiccatiSpecial]

EXPECTED = [('df - x*fx',
  ['separable',
   '1st_exact',
   '1st_linear',
   'Bernoulli',
   'almost_linear',
   '1st_power_series',
   'lie_group',
   'separable_Integral',
   '1st_exact_Integral',
   '1st_linear_Integral',
   'Bernoulli_Integral',
   'almost_linear_Integral'],
  {'1st_exact': ['-x*f(x)*exp(Integral(-x, x))', 'exp(Integral(-x, x))'],
   '1st_homogeneous_coeff_best': None,
   '1st_linear': ['-x', '0'],
   'Bernoulli': ['-x', '0', '0'],
   'Riccati_special_minus2': None,
   'almost_linear': ['1', '-x*f(x)'],
   'linear_coefficients': None,
   'separable': ['-x*f(x)', '1'],
   'separable_reduced': None}),
 ('fx*df - x',
  ['factorable',
   'separable',
   '1st_exact',
   'Bernoulli',
   '1st_homogeneous_coeff_best',
   '1st_homogeneous_coeff_subs_indep_div_dep',
   '1st_homogeneous_coeff_subs_dep_div_indep',
   '1st_power_series',
   'lie_group',
   'separable_Integral',
   '1st_exact_Integral',
   'Bernoulli_Integral',
   '1st_homogeneous_coeff_subs_indep_div_dep_Integral',
   '1st_homogeneous_coeff_subs_dep_div_indep_Integral'],
  {'1st_exact': ['-x', 'f(x)'],
   '1st_homogeneous_coeff_best': ['-x', 'f(x)'],
   '1st_linear': None,
   'Bernoulli': ['0', 'x', '-1'],
   'Riccati_special_minus2': None,
   'almost_linear': None,
   'linear_coefficients': None,
   'separable': ['-x', 'f(x)'],
   'separable_reduced': None}),
 ('x*df + fx - x**2',
  ['factorable',
   '1st_exact',
   '1st_linear',
   'Bernoulli',
   'almost_linear',
   'lie_group',
   'nth_linear_euler_eq_nonhomogeneous_undetermined_coefficients',
   'nth_linear_euler_eq_nonhomogeneous_variation_of_parameters',
   '1st_exact_Integral',
   '1st_linear_Integral',
   'Bernoulli_Integral',
   'almost_linear_Integral',
   'nth_linear_euler_eq_nonhomogeneous_variation_of_parameters_Integral'],
  {'1st_exact': ['-x**2 + f(x)', 'x'],
   '1st_homogeneous_coeff_best': None,
   '1st_linear': ['1/x', 'x'],
   'Bernoulli': ['1/x', 'x', '0'],
   'Riccati_special_minus2': None,
   'almost_linear': ['x', '-x**2 + f(x)'],
   'linear_coefficients': None,
   'separable': None,
   'separable_reduced': None}),
 ('df + fx/x - x*fx**3',
  ['Bernoulli', 'separable_reduced', 'lie_group', 'Bernoulli_Integral', 'separable_reduced_Integral'],
  {'1st_exact': None,
   '1st_homogeneous_coeff_best': None,
   '1st_linear': None,
   'Bernoulli': ['1/x', 'x', '3'],
   'Riccati_special_minus2': None,
   'almost_linear': None,
   'linear_coefficients': None,
   'separable': None,
   'separable_reduced': ['-x*f(x)**3 + f(x)/x', '1']}),
 ('(x**2 + fx**2)*df - x*fx',
  ['factorable',
   '1st_exact',
   '1st_homogeneous_coeff_best',
   '1st_homogeneous_coeff_subs_indep_div_dep',
   '1st_homogeneous_coeff_subs_dep_div_indep',
   '1st_power_series',
   'lie_group',
   '1st_exact_Integral',
   '1st_homogeneous_coeff_subs_indep_div_dep_Integral',
   '1st_homogeneous_coeff_subs_dep_div_indep_Integral'],
  {'1st_exact': ['-x*f(x)*exp(Integral(-3/f(x), f(x)))', '(x**2 + f(x)**2)*exp(Integral(-3/f(x), f(x)))'],
   '1st_homogeneous_coeff_best': ['-x*f(x)', 'x**2 + f(x)**2'],
   '1st_linear': None,
   'Bernoulli': None,
   'Riccati_special_minus2': None,
   'almost_linear': None,
   'linear_coefficients': None,
   'separable': None,
   'separable_reduced': None}),
 ('2*x*fx + (x**2 + 3*fx**2)*df',
  ['factorable',
   '1st_exact',
   '1st_homogeneous_coeff_best',
   '1st_homogeneous_coeff_subs_indep_div_dep',
   '1st_homogeneous_coeff_subs_dep_div_indep',
   '1st_power_series',
   'lie_group',
   '1st_exact_Integral',
   '1st_homogeneous_coeff_subs_indep_div_dep_Integral',
   '1st_homogeneous_coeff_subs_dep_div_indep_Integral'],
  {'1st_exact': ['2*x*f(x)', 'x**2 + 3*f(x)**2'],
   '1st_homogeneous_coeff_best': ['2*x*f(x)', 'x**2 + 3*f(x)**2'],
   '1st_linear': None,
   'Bernoulli': None,
   'Riccati_special_minus2': None,
   'almost_linear': None,
   'linear_coefficients': None,
   'separable': None,
   'separable_reduced': None}),
 ('df + 2*x*fx - exp(-x**2)',
  ['1st_exact',
   '1st_linear',
   'Bernoulli',
   'almost_linear',
   '1st_power_series',
   'lie_group',
   '1st_exact_Integral',
   '1st_linear_Integral',
   'Bernoulli_Integral',
   'almost_linear_Integral'],
  {'1st_exact': ['(2*x*f(x) - exp(-x**2))*exp(Integral(2*x, x))', 'exp(Integral(2*x, x))'],
   '1st_homogeneous_coeff_best': None,
   '1st_linear': ['2*x', 'exp(-x**2)'],
   'Bernoulli': ['2*x', 'exp(-x**2)', '0'],
   'Riccati_special_minus2': None,
   'almost_linear': ['1', '2*x*f(x) - exp(-x**2)'],
   'linear_coefficients': None,
   'separable': None,
   'separable_reduced': None}),
 ('x*df - fx - sqrt(x**2 + fx**2)',
  ['factorable',
   '1st_homogeneous_coeff_best',
   '1st_homogeneous_coeff_subs_indep_div_dep',
   '1st_homogeneous_coeff_subs_dep_div_indep',
   'lie_group',
   '1st_homogeneous_coeff_subs_indep_div_dep_Integral',
   '1st_homogeneous_coeff_subs_dep_div_indep_Integral'],
  {'1st_exact': None,
   '1st_homogeneous_coeff_best': ['-sqrt(x**2 + f(x)**2) - f(x)', 'x'],
   '1st_linear': None,
   'Bernoulli': None,
   'Riccati_special_minus2': None,
   'almost_linear': None,
   'linear_coefficients': None,
   'separable': None,
   'separable_reduced': None}),
 ('(x + fx + 1)*df - (x - fx + 3)',
  ['factorable',
   '1st_exact',
   'linear_coefficients',
   '1st_power_series',
   'lie_group',
   '1st_exact_Integral',
   'linear_coefficients_Integral'],
  {'1st_exact': ['-x + f(x) - 3', 'x + f(x) + 1'],
   '1st_homogeneous_coeff_best': None,
   '1st_linear': None,
   'Bernoulli': None,
   'Riccati_special_minus2': None,
   'almost_linear': None,
   'linear_coefficients': ['-x + f(x) - 3', 'x + f(x) + 1'],
   'separable': None,
   'separable_reduced': None}),
 ('a*df - b*fx**2 - 1/x**2',
  ['factorable', 'Riccati_special_minus2', 'lie_group'],
  {'1st_exact': None,
   '1st_homogeneous_coeff_best': None,
   '1st_linear': None,
   'Bernoulli': None,
   'Riccati_special_minus2': ['a', '-b', '0', '-1'],
   'almost_linear': None,
   'linear_coefficients': None,
   'separable': None,
   'separable_reduced': None}),
 ('df - sin(x)*cos(fx)',
  ['separable', '1st_exact', '1st_power_series', 'lie_group', 'separable_Integral', '1st_exact_Integral'],
  {'1st_exact': ['-exp(Integral(sin(f(x))/cos(f(x)), f(x)))*sin(x)*cos(f(x))',
                 'exp(Integral(sin(f(x))/cos(f(x)), f(x)))'],
   '1st_homogeneous_coeff_best': None,
   '1st_linear': None,
   'Bernoulli': None,
   'Riccati_special_minus2': None,
   'almost_linear': None,
   'linear_coefficients': None,
   'separable': ['-sin(x)*cos(f(x))', '1'],
   'separable_reduced': None}),
 ('exp(fx)*df + 2*x*exp(fx) - x',
  ['factorable',
   'separable',
   '1st_exact',
   'almost_linear',
   '1st_power_series',
   'lie_group',
   'separable_Integral',
   '1st_exact_Integral',
   'almost_linear_Integral'],
  {'1st_exact': ['(2*x*exp(f(x)) - x)*exp(Integral(2*x, x))', 'exp(f(x))*exp(Integral(2*x, x))'],
   '1st_homogeneous_coeff_best': None,
   '1st_linear': None,
   'Bernoulli': None,
   'Riccati_special_minus2': None,
   'almost_linear': ['exp(f(x))', '2*x*exp(f(x)) - x'],
   'linear_coefficients': None,
   'separable': ['2*x*exp(f(x)) - x', 'exp(f(x))'],
   'separable_reduced': None}),
 ('cos(fx)*df + sin(fx) - x',
  ['factorable',
   '1st_exact',
   'almost_linear',
   '1st_power_series',
   'lie_group',
   '1st_exact_Integral',
   'almost_linear_Integral'],
  {'1st_exact': ['(-x + sin(f(x)))*exp(Integral(1, x))', 'exp(Integral(1, x))*cos(f(x))'],
   '1st_homogeneous_coeff_best': None,
   '1st_linear': None,
   'Bernoulli': None,
   'Riccati_special_minus2': None,
   'almost_linear': ['cos(f(x))', '-x + sin(f(x))'],
   'linear_coefficients': None,
   'separable': None,
   'separable_reduced': None}),
 ('df - a*fx',
  ['separable',
   '1st_exact',
   '1st_linear',
   'Bernoulli',
   'almost_linear',
   '1st_power_series',
   'lie_group',
   'nth_linear_constant_coeff_homogeneous',
   'separable_Integral',
   '1st_exact_Integral',
   '1st_linear_Integral',
   'Bernoulli_Integral',
   'almost_linear_Integral'],
  {'1st_exact': ['-a*f(x)*exp(Integral(-a, x))', 'exp(Integral(-a, x))'],
   '1st_homogeneous_coeff_best': None,
   '1st_linear': ['-a', '0'],
   'Bernoulli': ['-a', '0', '0'],
   'Riccati_special_minus2': None,
   'almost_linear': ['1', '-a*f(x)'],
   'linear_coefficients': None,
   'separable': ['-a*f(x)', '1'],
   'separable_reduced': None}),
 ('x**2*df - fx**2 + x*fx',
  ['factorable',
   'Bernoulli',
   '1st_rational_riccati',
   '1st_homogeneous_coeff_best',
   '1st_homogeneous_coeff_subs_indep_div_dep',
   '1st_homogeneous_coeff_subs_dep_div_indep',
   'lie_group',
   'Bernoulli_Integral',
   '1st_homogeneous_coeff_subs_indep_div_dep_Integral',
   '1st_homogeneous_coeff_subs_dep_div_indep_Integral'],
  {'1st_exact': None,
   '1st_homogeneous_coeff_best': ['x*f(x) - f(x)**2', 'x**2'],
   '1st_linear': None,
   'Bernoulli': ['1/x', 'x**(-2)', '2'],
   'Riccati_special_minus2': None,
   'almost_linear': None,
   'linear_coefficients': None,
   'separable': None,
   'separable_reduced': None}),
 ('df + fx**2 - 2/x**2',
  ['1st_rational_riccati',
   'Riccati_special_minus2',
   'separable_reduced',
   'lie_group',
   'separable_reduced_Integral'],
  {'1st_exact': None,
   '1st_homogeneous_coeff_best': None,
   '1st_linear': None,
   'Bernoulli': None,
   'Riccati_special_minus2': ['1', '1', '0', '-2'],
   'almost_linear': None,
   'linear_coefficients': None,
   'separable': None,
   'separable_reduced': ['f(x)**2 - 2/x**2', '1']}),
 ('x*df - fx*(log(x*fx) - 1)',
  ['factorable', 'separable_reduced', 'lie_group', 'separable_reduced_Integral'],
  {'1st_exact': None,
   '1st_homogeneous_coeff_best': None,
   '1st_linear': None,
   'Bernoulli': None,
   'Riccati_special_minus2': None,
   'almost_linear': None,
   'linear_coefficients': None,
   'separable': None,
   'separable_reduced': ['(1 - log(x*f(x)))*f(x)', 'x']}),
 ('df - (x + fx)**2',
  ['factorable', '1st_rational_riccati', '1st_power_series', 'lie_group'],
  {'1st_exact': None,
   '1st_homogeneous_coeff_best': None,
   '1st_linear': None,
   'Bernoulli': None,
   'Riccati_special_minus2': None,
   'almost_linear': None,
   'linear_coefficients': None,
   'separable': None,
   'separable_reduced': None}),
 ('df + fx*tan(x) - sec(x)',
  ['1st_exact',
   '1st_linear',
   'Bernoulli',
   'almost_linear',
   '1st_power_series',
   'lie_group',
   '1st_exact_Integral',
   '1st_linear_Integral',
   'Bernoulli_Integral',
   'almost_linear_Integral'],
  {'1st_exact': ['(f(x)*tan(x) - sec(x))*exp(Integral(tan(x), x))', 'exp(Integral(tan(x), x))'],
   '1st_homogeneous_coeff_best': None,
   '1st_linear': ['tan(x)', 'sec(x)'],
   'Bernoulli': ['tan(x)', 'sec(x)', '0'],
   'Riccati_special_minus2': None,
   'almost_linear': ['1', 'f(x)*tan(x) - sec(x)'],
   'linear_coefficients': None,
   'separable': None,
   'separable_reduced': None}),
 ('(2*x*fx + 1)*df + fx**2',
  ['factorable',
   '1st_exact',
   'separable_reduced',
   '1st_power_series',
   'lie_group',
   '1st_exact_Integral',
   'separable_reduced_Integral'],
  {'1st_exact': ['f(x)**2', '2*x*f(x) + 1'],
   '1st_homogeneous_coeff_best': None,
   '1st_linear': None,
   'Bernoulli': None,
   'Riccati_special_minus2': None,
   'almost_linear': None,
   'linear_coefficients': None,
   'separable': None,
   'separable_reduced': ['f(x)**2', '2*x*f(x) + 1']})]


def _eq(source):
    return eval(source, NS)


@pytest.mark.parametrize('source,hints,rows', EXPECTED, ids=[c[0] for c in EXPECTED])
def test_classify_ode_matches_sympy(source, hints, rows):
    assert list(S.classify_ode(_eq(source), fx)) == hints


@pytest.mark.parametrize('source,hints,rows', EXPECTED, ids=[c[0] for c in EXPECTED])
def test_solver_matches_match_sympy(source, hints, rows):
    eq = _eq(source)
    for cls in CLASSES:
        solver = cls(SingleODEProblem(eq, fx, x))
        expected = rows[cls.hint]
        if expected is None:
            assert not solver.matches(), cls.hint
        else:
            assert solver.matches(), cls.hint
            assert solver.wilds_match() == [S.sympify(v, locals=NS) for v in expected], cls.hint


@pytest.mark.parametrize('source,hints,rows', EXPECTED, ids=[c[0] for c in EXPECTED])
def test_native_agrees_with_reference_matcher(source, hints, rows, monkeypatch):
    eq = _eq(source)
    for cls in CLASSES:
        native = cls(SingleODEProblem(eq, fx, x))
        decided = _classify_native.solver_match(native)
        if decided is None:
            continue
        native_ok = native.matches()
        with monkeypatch.context() as m:
            m.setattr(_classify_native, 'solver_match', lambda solver: None)
            reference = cls(SingleODEProblem(eq, fx, x))
            reference_ok = reference.matches()
        assert native_ok == reference_ok, cls.hint
        if native_ok:
            assert native.wilds_match() == reference.wilds_match(), cls.hint


def test_native_kernel_decides_common_shapes():
    for source in ['df - x*fx', 'x*df + fx - x**2', '2*x*fx + (x**2 + 3*fx**2)*df']:
        problem = SingleODEProblem(_eq(source), fx, x)
        native = _classify_native._first_order(problem)
        assert native is not None
        assert native['de'] is not None
    native = _classify_native._first_order(SingleODEProblem(_eq('df + 2*x*fx - exp(-x**2)'), fx, x))
    assert native['linear'] == (2*x, S.exp(-x**2))
    assert native['bernoulli'] == (2*x, S.exp(-x**2), 0)
    native = _classify_native._first_order(SingleODEProblem(_eq('df + fx/x - x*fx**3'), fx, x))
    assert native['linear'] is False
    assert native['bernoulli'] == (1/x, x, 3)


def test_unmodelled_content_declines():
    for source in ['df - I*fx', 'df - oo*fx', 'df - 1.5*fx']:
        native = _classify_native._first_order(SingleODEProblem(_eq(source), fx, x))
        assert native is None or all(v is None for v in native.values()), source


def test_ordered_hints_follow_allhints():
    from symox.solvers.ode.ode import allhints
    assert _core.ode_ordered_hints(['separable', '1st_power_series', 'factorable', 'order']) == [
        'factorable', 'separable', '1st_power_series', 'separable_Integral']
    assert _core.ode_ordered_hints(list(allhints)) == list(allhints)


def test_split_de_power_series_match():
    collected = S.collect(S.expand(x*df + fx), df, exact=True)
    d, e = S.Wild('d', exclude=[df, fx.diff(x, 2)]), S.Wild('e', exclude=[df])
    assert _classify_native.split_de(collected, df, d, e, fx.diff(x, 2)) == {d: fx, e: x}


# SymPy 1.14 match_riccati(SingleODEProblem(eq, f(x), x).eq_expanded, f, x).
RICCATI = [
    ('df - (x**2 + 1) - fx**2', True, ['x**2 + 1', '0', '1']),
    ('df - 2*fx**2/x - 3/x**2', True, ['3/x**2', '0', '2/x']),
    ('2*df + 3*fx**2 - 1/x**2', True, ['1/(2*x**2)', '0', '-3/2']),
    ('x*df - fx**2 - 3*fx + 4', True, ['-4/x', '3/x', '1/x']),
    ('(x**2 + 1)*df + fx**2 - x', True, ['x/(x**2 + 1)', '0', '-1/(x**2 + 1)']),
    ('df - fx**2 - a*x', False, []),
    ('df - exp(x)*fx**2 - 1', False, []),
    ('df - sqrt(2)*fx**2 - x', True, ['x', '0', 'sqrt(2)']),
    ('df - fx**3 - 1', False, []),
    ('df - x*fx', False, []),
    ('df*x**2 - x**2*fx**2 - 2', True, ['2/x**2', '0', '1']),
    ('df + (x**2 + 2)*fx**2/(x - 1) - x', True, ['x', '0', '-(x**2 + 2)/(x - 1)']),
    ('df - 1.5*fx**2', False, []),
]


@pytest.mark.parametrize('source,match,funcs', RICCATI, ids=[c[0] for c in RICCATI])
def test_riccati_match_matches_sympy(source, match, funcs):
    # Compared through str: SymPy's printed forms of the coefficients.
    from symox.solvers.ode.riccati import match_riccati
    eq = SingleODEProblem(_eq(source), fx, x).eq_expanded
    native = _classify_native.riccati_match(eq, f, x)
    if native is not None:
        assert (native[0], [str(v) for v in native[1]]) == (match, funcs)
    reference = match_riccati(eq, f, x)
    assert (reference[0], [str(v) for v in reference[1]]) == (match, funcs)
    solver = RationalRiccati(SingleODEProblem(_eq(source), fx, x))
    assert solver.matches() == match
    if match:
        assert [str(v) for v in solver.wilds_match()] == funcs


def test_riccati_match_decided_natively():
    eq = SingleODEProblem(_eq('x*df - fx**2 - 3*fx + 4'), fx, x).eq_expanded
    assert _classify_native.riccati_match(eq, f, x) == (True, [-4/x, 3/x, 1/x])
    eq = SingleODEProblem(_eq('df - fx**2 - a*x'), fx, x).eq_expanded
    assert _classify_native.riccati_match(eq, f, x) == (False, [])
    for source in ['df - I*fx**2 - x', 'df - 1.5*fx**2', 'df - oo*fx**2']:
        eq = SingleODEProblem(_eq(source), fx, x).eq_expanded
        assert _classify_native.riccati_match(eq, f, x) is None, source
