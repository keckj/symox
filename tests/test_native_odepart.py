"""Differential tests of the native particular-solution kernels of the linear
ODE solvers (symox::ode::undetermined, symox::ode::variation) against
SymPy 1.14.

The expected strings are SymPy 1.14 results of
dsolve(eq, f(x), hint=hint, simplify=simp). The other tests check that the
kernels decline the inputs outside their domain (the reference body then
runs unchanged), that the native results equal the reference bodies run in
symox with the kernels removed, and that the reference's errors are kept."""
import pytest

import symox as S
from symox import _core
from symox.solvers.ode import _odesolve_native as nat
from symox.solvers.ode import nonhomogeneous as nh

x, a = S.symbols('x a')
f = S.Function('f')
fx = f(x)
exp, sin, cos, log = S.exp, S.sin, S.cos, S.log


def d(n):
    return fx.diff(x, n)


U = 'nth_linear_constant_coeff_undetermined_coefficients'
V = 'nth_linear_constant_coeff_variation_of_parameters'
EU = 'nth_linear_euler_eq_nonhomogeneous_undetermined_coefficients'
EV = 'nth_linear_euler_eq_nonhomogeneous_variation_of_parameters'
KERNELS = ('ode_uc_trialset', 'ode_uc_solve', 'ode_vop_solve')

CASES = [
    ('uc_res_sin', lambda: d(2) + fx - sin(x), U),
    ('uc_res_exp2', lambda: d(2) - 2*d(1) + fx - x*exp(x), U),
    ('uc_mixed', lambda: d(2) + 3*d(1) + 2*fx - x**2*exp(-x) - sin(2*x), U),
    ('uc_3rd', lambda: d(3) - d(1) - x*exp(x)*cos(x) - 3, U),
    ('uc_res4', lambda: d(4) + 2*d(2) + fx - x**2*sin(x), U),
    ('uc_poly', lambda: d(2) - fx - x**3 - 2*x + 1, U),
    ('uc_half', lambda: 2*d(2) - fx - exp(x/2)*cos(x/2)/3, U),
    ('uc_const', lambda: d(1) - 5, U),
    ('vop_sec', lambda: d(2) + fx - 1/cos(x), V),
    ('vop_expx', lambda: d(2) - 2*d(1) + fx - exp(x)/x, V),
    ('vop_3rd', lambda: d(3) - d(1) - exp(2*x), V),
    ('vop_lead', lambda: 3*d(2) - 3*fx - x, V),
    ('eu_uc', lambda: x**2*d(2) - 2*x*d(1) + 2*fx - log(x), EU),
    ('eu_uc3', lambda: x**2*d(2) - 2*x*d(1) + 2*fx - x**3, EU),
    ('eu_vop', lambda: x**2*d(2) - 2*x*d(1) + 2*fx - x**4, EV),
    ('eu_vop_log', lambda: x**2*d(2) + x*d(1) - fx - log(x), EV),
]

# SymPy 1.14: str(dsolve(eq, f(x), hint=hint, simplify=simp))
SYMPY = {
    'uc_res_sin/True': 'Eq(f(x), C2*sin(x) + (C1 - x/2)*cos(x))',
    'uc_res_sin/False': 'Eq(f(x), C1*sin(x) + C2*cos(x) - x*cos(x)/2)',
    'uc_res_exp2/True': 'Eq(f(x), (C1 + x*(C2 + x**2/6))*exp(x))',
    'uc_res_exp2/False': 'Eq(f(x), C1*exp(x) + C2*x*exp(x) + x**3*exp(x)/6)',
    'uc_mixed/True': 'Eq(f(x), C2*exp(-2*x) + (C1 + x**3/3 - x**2 + 2*x)*exp(-x) - sin(2*x)/20'
                     ' - 3*cos(2*x)/20)',
    'uc_mixed/False': 'Eq(f(x), C1*exp(-2*x) + C2*exp(-x) + x**3*exp(-x)/3 - x**2*exp(-x) + 2*x*exp(-x)'
                      ' - sin(2*x)/20 - 3*cos(2*x)/20)',
    'uc_3rd/True': 'Eq(f(x), C1 + C3*exp(-x) - 3*x + (C2 + x*sin(x)/10 - 3*x*cos(x)/10 + 21*sin(x)/50'
                   ' + 11*cos(x)/25)*exp(x))',
    'uc_3rd/False': 'Eq(f(x), C1*exp(-x) + C2 + C3*exp(x) + x*exp(x)*sin(x)/10 - 3*x*exp(x)*cos(x)/10'
                    ' - 3*x + 21*exp(x)*sin(x)/50 + 11*exp(x)*cos(x)/25)',
    'uc_res4/True': 'Eq(f(x), (C1 + x*(C2 - x**2/12))*cos(x) + (C3 + x*(C4 - x**3/48 + 3*x/16))*sin(x))',
    'uc_res4/False': 'Eq(f(x), C1*sin(x) + C2*cos(x) + C3*x*sin(x) + C4*x*cos(x) - x**4*sin(x)/48'
                     ' - x**3*cos(x)/12 + 3*x**2*sin(x)/16)',
    'uc_poly/True': 'Eq(f(x), C1*exp(-x) + C2*exp(x) - x**3 - 8*x + 1)',
    'uc_poly/False': 'Eq(f(x), C1*exp(-x) + C2*exp(x) - x**3 - 8*x + 1)',
    'uc_half/True': 'Eq(f(x), C1*exp(-sqrt(2)*x/2) + C2*exp(sqrt(2)*x/2) + exp(x/2)*sin(x/2)/6'
                    ' - exp(x/2)*cos(x/2)/6)',
    'uc_half/False': 'Eq(f(x), C1*exp(-sqrt(2)*x/2) + C2*exp(sqrt(2)*x/2) + exp(x/2)*sin(x/2)/6'
                     ' - exp(x/2)*cos(x/2)/6)',
    'uc_const/True': 'Eq(f(x), C1 + 5*x)',
    'uc_const/False': 'Eq(f(x), C1 + 5*x)',
    'vop_sec/True': 'Eq(f(x), (C1 + x)*sin(x) + (C2 + log(cos(x)))*cos(x))',
    'vop_sec/False': 'Eq(f(x), C1*sin(x) + C2*cos(x) + x*sin(x)/(sin(x)**2 + cos(x)**2)'
                     ' - (-log(tan(x/2) - 1) - log(tan(x/2) + 1) + log(tan(x/2)**2 + 1))*cos(x))',
    'vop_expx/True': 'Eq(f(x), (C1 + x*(C2 + log(x)))*exp(x))',
    'vop_expx/False': 'Eq(f(x), C1*exp(x) + C2*x*exp(x) + x*exp(x)*log(x) - x*exp(x))',
    'vop_3rd/True': 'Eq(f(x), C1 + C2*exp(-x) + C3*exp(x) + exp(2*x)/6)',
    'vop_3rd/False': 'Eq(f(x), C1*exp(-x) + C2 + C3*exp(x) + exp(2*x)/6)',
    'vop_lead/True': 'Eq(f(x), C1*exp(-x) + C2*exp(x) - x/3)',
    'vop_lead/False': 'Eq(f(x), C1*exp(-x) + C2*exp(x) - x/3)',
    'eu_uc/True': 'Eq(f(x), C1*x + C2*x**2 + log(x)/2 + 3/4)',
    'eu_uc/False': 'Eq(f(x), C1*x + C2*x**2 + log(x)/2 + 3/4)',
    'eu_uc3/True': 'Eq(f(x), x*(C1 + C2*x + x**2/2))',
    'eu_uc3/False': 'Eq(f(x), C1*x + C2*x**2 + x**3/2)',
    'eu_vop/True': 'Eq(f(x), x*(C1 + C2*x + x**3/6))',
    'eu_vop/False': 'Eq(f(x), C1*x + C2*x**2 + x**4/6)',
    'eu_vop_log/True': 'Eq(f(x), C1/x + C2*x - log(x))',
    'eu_vop_log/False': 'Eq(f(x), C1/x + C2*x + x**2*((-log(x)/(2*x) - 1/(2*x))/x'
                        ' - (x*log(x)/2 - x/2)/x**3))',
}

# The held variation-of-parameters integrals are SymPy's; symox evaluates
# Integral(-sin(x)/((-sin(x)**2 - cos(x)**2)*cos(x)), x) to -log(cos(x)) and
# Integral(-1/(-sin(x)**2 - cos(x)**2), x) to x, before and after the port
# (integration gap, docs/algorithmic-parity.md).
INTEGRATION_GAP = {'vop_sec/False'}

PARAMS = [pytest.param(name, eq, hint, simp, id='%s/%s' % (name, simp),
                       marks=[pytest.mark.xfail(strict=True, reason='integration gap')]
                       if '%s/%s' % (name, simp) in INTEGRATION_GAP else [])
          for name, eq, hint in CASES for simp in (True, False)]
ALL = [pytest.param(name, eq, hint, simp, id='%s/%s' % (name, simp))
       for name, eq, hint in CASES for simp in (True, False)]


@pytest.mark.parametrize('name,eq,hint,simp', PARAMS)
def test_dsolve_matches_sympy(name, eq, hint, simp):
    assert str(S.dsolve(eq(), fx, hint=hint, simplify=simp)) == SYMPY['%s/%s' % (name, simp)]


def test_held_integrals_match_sympy():
    sol = S.dsolve(d(2) + fx - 1/cos(x), fx, hint=V + '_Integral', simplify=False)
    assert str(sol) == ('Eq(f(x), C1*sin(x) + C2*cos(x) + sin(x)*Integral(-1/(-sin(x)**2 - cos(x)**2), x)'
                        ' - cos(x)*Integral(-sin(x)/((-sin(x)**2 - cos(x)**2)*cos(x)), x))')


def test_wronskian_keeps_exp_negative_argument_as_factor():
    # dotprodsimp hides exp(-x) behind a Dummy before cancel (SymPy 1.14),
    # so the Wronskian and the held integrands keep SymPy's expanded form
    w = S.wronskian([exp(-x)*sin(3*x), exp(-x)*cos(3*x)], x)
    assert str(w) == '-3*exp(-2*x)*sin(3*x)**2 - 3*exp(-2*x)*cos(3*x)**2'
    sol = S.dsolve(d(2) + 2*d(1) + 10*fx - exp(x)/x, fx, hint=V + '_Integral', simplify=False)
    assert str(sol) == (
        'Eq(f(x), C1*exp(-x)*sin(3*x) + C2*exp(-x)*cos(3*x) + exp(-x)*sin(3*x)*Integral(-cos(3*x)/'
        '(x*(-3*exp(-2*x)*sin(3*x)**2 - 3*exp(-2*x)*cos(3*x)**2)), x) - exp(-x)*cos(3*x)*Integral('
        '-sin(3*x)/(x*(-3*exp(-2*x)*sin(3*x)**2 - 3*exp(-2*x)*cos(3*x)**2)), x))')


def _stock(monkeypatch):
    for k in KERNELS:
        monkeypatch.delattr(_core, k)


@pytest.mark.parametrize('name,eq,hint,simp', ALL)
def test_native_equals_reference_body(monkeypatch, name, eq, hint, simp):
    native = S.dsolve(eq(), fx, hint=hint, simplify=simp)
    _stock(monkeypatch)
    assert S.dsolve(eq(), fx, hint=hint, simplify=simp) == native


MATCH = [exp(x), x**2*exp(-3*x)*sin(2*x), sin(x) + cos(x), x**3 + 2*x - 1, exp(x/2)*cos(x/3)/7,
         S.Integer(4), 3*x*exp(2*x) - exp(2*x)*cos(x), sin(-x), exp(x)*exp(2*x)**2]


@pytest.mark.parametrize('expr', MATCH, ids=[str(e) for e in MATCH])
def test_trialset_equals_reference_body(monkeypatch, expr):
    native = nh._undetermined_coefficients_match(expr, x)
    hom = (d(2) + 9*fx)
    native_h = nh._undetermined_coefficients_match(expr, x, func=fx, eq_homogeneous=hom)
    assert native['test'] and native_h['test']
    _stock(monkeypatch)
    assert nh._undetermined_coefficients_match(expr, x) == native
    assert nh._undetermined_coefficients_match(expr, x, func=fx, eq_homogeneous=hom) == native_h


def test_trialset_values():
    assert nh._undetermined_coefficients_match(9*x*exp(x) + exp(-x), x) == \
        {'test': True, 'trialset': {x*exp(x), exp(-x), exp(x)}}
    # resonance with the homogeneous solutions sin(x), cos(x)
    assert nh._undetermined_coefficients_match(sin(x), x, func=fx, eq_homogeneous=d(2) + fx) == \
        {'test': True, 'trialset': {x*sin(x), x*cos(x)}}


@pytest.mark.parametrize('expr', [
    log(x),              # outside the method
    x*(x + 1),           # an unexpanded product (the reference's mutable default argument)
    S.Float(1.5)*exp(x),  # inexact coefficient
    a*exp(x),            # symbolic coefficient
    exp(a*x),            # symbolic frequency
    S.sqrt(x)*exp(x),    # non-integer power of x
    exp(x**2),           # non-linear exponent
    sin(x)*cos(x),       # two trigonometric factors
    S.sqrt(2)*sin(x),    # algebraic coefficient
    sin(x)**2,           # power of a trigonometric factor
    1/x,
])
def test_trialset_kernel_declines(expr):
    assert _core.ode_uc_trialset(expr, x, x, None) is None


def test_trialset_declined_inputs_keep_reference(monkeypatch):
    exprs = [log(x), x*(x + 1), a*exp(x), sin(x)**2, 1/x, S.Float(1.5)*exp(x)]
    native = [nh._undetermined_coefficients_match(e, x) for e in exprs]
    _stock(monkeypatch)
    assert [nh._undetermined_coefficients_match(e, x) for e in exprs] == native


def test_solve_kernel_declines():
    eq = d(2) + fx - sin(x)
    terms = [x*sin(x), x*cos(x)]
    rhs = S.Symbol('C1')*sin(x) + S.Symbol('C2')*cos(x)
    assert _core.ode_uc_solve(eq, fx, x, 2, 2, rhs, terms) == rhs - x*cos(x)/2
    # number of homogeneous solutions differs from the order
    assert _core.ode_uc_solve(eq, fx, x, 2, 1, rhs, terms) is None
    # singular system: the trial set misses the resonance
    assert _core.ode_uc_solve(eq, fx, x, 2, 2, rhs, [sin(x), cos(x)]) is None
    # symbolic coefficients in the equation
    assert _core.ode_uc_solve(d(2) + a*fx - sin(x), fx, x, 2, 2, rhs, terms) is None
    # trial term that is not a canonical monomial
    assert _core.ode_uc_solve(eq, fx, x, 2, 2, rhs, [x*sin(x), x*cos(x), log(x)]) is None


def test_singular_system_keeps_reference_error(monkeypatch):
    # the reference raises NotImplementedError when solve finds no values
    match = {'list': [sin(x), cos(x)], 'sol': S.Eq(fx, S.Symbol('C1')*sin(x) + S.Symbol('C2')*cos(x))}
    eq = d(2) + fx - sin(x)
    with pytest.raises(NotImplementedError) as native:
        nh._solve_undetermined_coefficients(eq, fx, 2, match, {sin(x), cos(x)})
    _stock(monkeypatch)
    with pytest.raises(NotImplementedError) as stock:
        nh._solve_undetermined_coefficients(eq, fx, 2, match, {sin(x), cos(x)})
    assert str(native.value) == str(stock.value)


def test_dependent_solutions_keep_reference_error(monkeypatch):
    sol = S.Eq(fx, S.Symbol('C1')*exp(x) + S.Symbol('C2')*exp(x))
    r = {2: S.Integer(1), 1: S.Integer(-2), 0: S.Integer(1), -1: exp(x)/x}
    with pytest.raises(NotImplementedError) as native:
        nh._solve_variation_of_parameters(d(2) - 2*d(1) + fx - exp(x)/x, fx, [exp(x), 2*exp(x)], sol, 2, r, False)
    _stock(monkeypatch)
    with pytest.raises(NotImplementedError) as stock:
        nh._solve_variation_of_parameters(d(2) - 2*d(1) + fx - exp(x)/x, fx, [exp(x), 2*exp(x)], sol, 2, r, False)
    assert str(native.value) == str(stock.value)


def test_vop_callback_errors_propagate():
    sol = S.Eq(fx, S.Symbol('C1')*exp(x) + S.Symbol('C2')*x*exp(x))
    r = {2: S.Integer(1), 1: S.Integer(-2), 0: S.Integer(1), -1: exp(x)/x}

    def boom(e):
        raise ZeroDivisionError('callback')
    with pytest.raises(ZeroDivisionError):
        nat.variation_of_parameters(fx, [exp(x), x*exp(x)], sol, 2, r, True, boom, boom)
    # without simplification the callbacks are not called
    rhs = nat.variation_of_parameters(fx, [exp(x), x*exp(x)], sol, 2, r, False, boom, boom)
    assert rhs is not None and rhs.has(S.Integral)


def test_vop_kernel_declines():
    sol = S.Eq(fx, S.Symbol('C1')*exp(x) + S.Symbol('C2')*x*exp(x))
    r = {2: S.Integer(1), 1: S.Integer(-2), 0: S.Integer(1), -1: exp(x)/x}
    ident = lambda e: e  # noqa: E731
    # order differs from the number of solutions
    assert nat.variation_of_parameters(fx, [exp(x)], sol, 2, r, False, ident, ident) is None
    # a Piecewise integrand is folded by the reference Integral constructor
    rp = dict(r)
    rp[-1] = S.Piecewise((x, x > 0), (0, True))
    assert nat.variation_of_parameters(fx, [exp(x), x*exp(x)], sol, 2, rp, False, ident, ident) is None
    # nan forcing
    rn = dict(r)
    rn[-1] = S.nan
    assert nat.variation_of_parameters(fx, [exp(x), x*exp(x)], sol, 2, rn, False, ident, ident) is None
