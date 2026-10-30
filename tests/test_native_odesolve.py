"""Differential tests of the native characteristic-equation kernels of the
constant-coefficient and Cauchy-Euler linear ODE solvers
(symox::ode::constcoeff) against SymPy 1.14.

The expected strings are SymPy 1.14 results of
dsolve(eq, f(x), hint=hint, simplify=simp). The other tests check that the
kernels decline the inputs on which the native construction is not known to
reproduce the reference (the reference body then runs unchanged) and that
the native results equal the reference bodies run in symox."""
import pytest

import symox as S
from symox import _core
from symox.solvers.ode import _odesolve_native as nat
from symox.solvers.ode import nonhomogeneous as nh

x, a = S.symbols('x a')
k = S.Symbol('k', real=True)
f = S.Function('f')
fx = f(x)


def d(n):
    return fx.diff(x, n)


H = 'nth_linear_constant_coeff_homogeneous'
E = 'nth_linear_euler_eq_homogeneous'
U = 'nth_linear_constant_coeff_undetermined_coefficients'
V = 'nth_linear_constant_coeff_variation_of_parameters'
EU = 'nth_linear_euler_eq_nonhomogeneous_undetermined_coefficients'
EV = 'nth_linear_euler_eq_nonhomogeneous_variation_of_parameters'

CASES = [
    ('d2m5', lambda: d(2) - 5*fx, H, False),
    ('d2m5s', lambda: d(2) - 5*fx, H, True),
    ('d6', lambda: d(6) - 7*d(4) + 16*d(2) - 12*fx, H, False),
    ('d6s', lambda: d(6) - 7*d(4) + 16*d(2) - 12*fx, H, True),
    ('d3', lambda: d(3) - d(2) - 14*d(1) - 6*fx, H, False),
    ('sym_a', lambda: d(2) + a*fx, H, False),
    ('sym_a2', lambda: d(2) - a**2*fx, H, False),
    ('sym_a4', lambda: d(4) - a**4*fx, H, False),
    ('sym_k', lambda: d(2) + 2*k*d(1) + k**2*fx, H, True),
    ('quint', lambda: d(5) + 10*d(1) - 2*fx, H, True),
    ('quint_ns', lambda: d(5) + 10*d(1) - 2*fx, H, False),
    ('cplx_rep', lambda: d(4) + 2*d(2) + fx, H, True),
    ('cplx_rep_ns', lambda: d(4) + 2*d(2) + fx, H, False),
    ('mixed', lambda: d(3) - 3*d(2) + 7*d(1) - 5*fx, H, True),
    ('rat', lambda: 6*d(2) - 5*d(1) + fx, H, True),
    ('d1', lambda: d(1) - 3*fx, H, True),
    ('eu_cplx', lambda: x**2*d(2) - 3*x*d(1) + 13*fx, E, True),
    ('eu_cplx_ns', lambda: x**2*d(2) - 3*x*d(1) + 13*fx, E, False),
    ('eu_3', lambda: x**3*d(3) - x*d(1) + 5*fx, E, True),
    ('eu_rep', lambda: x**2*d(2) - 3*x*d(1) + 4*fx, E, True),
    ('eu_rat', lambda: 4*x**2*d(2) + fx, E, True),
    ('eu_irr', lambda: x**2*d(2) - 2*x*d(1) - 2*fx, E, True),
    ('eu_ln', lambda: x**2*d(2) + x*d(1) + fx, E, True),
    ('und', lambda: d(2) + 3*d(1) + 2*fx - S.exp(x), U, True),
    ('und_sqrt', lambda: d(2) - 2*fx - x**2, U, True),
    ('vop', lambda: d(2) + fx - 1/S.cos(x), V, True),
    ('vop_sqrt', lambda: d(2) - 3*fx - S.exp(x), V, True),
    ('eu_und', lambda: x**2*d(2) - 2*x*d(1) + 2*fx - x**3, EU, True),
    ('eu_var', lambda: x**2*d(2) - 2*x*d(1) + 2*fx - x**4, EV, True),
]

# SymPy 1.14: str(dsolve(eq, f(x), hint=hint, simplify=simp))
SYMPY = {
    'd2m5': 'Eq(f(x), C1*exp(sqrt(5)*x) + C2*exp(-sqrt(5)*x))',
    'd2m5s': 'Eq(f(x), C1*exp(-sqrt(5)*x) + C2*exp(sqrt(5)*x))',
    'd6': 'Eq(f(x), C1*exp(sqrt(2)*x) + C2*x*exp(sqrt(2)*x) + C3*exp(sqrt(3)*x) + C4*exp(-sqrt(2)*x)'
          ' + C5*x*exp(-sqrt(2)*x) + C6*exp(-sqrt(3)*x))',
    'd6s': 'Eq(f(x), C5*exp(-sqrt(3)*x) + C6*exp(sqrt(3)*x) + (C1 + C2*x)*exp(-sqrt(2)*x)'
           ' + (C3 + C4*x)*exp(sqrt(2)*x))',
    'd3': 'Eq(f(x), C1*exp(-3*x) + C2*exp(x*(2 + sqrt(6))) + C3*exp(x*(2 - sqrt(6))))',
    'sym_a': 'Eq(f(x), C1*exp(x*sqrt(-a)) + C2*exp(-x*sqrt(-a)))',
    'sym_a2': 'Eq(f(x), C1*exp(a*x) + C2*exp(-a*x))',
    'sym_a4': 'Eq(f(x), C1*exp(a*x) + C2*exp(-a*x) + C3*exp(I*a*x) + C4*exp(-I*a*x))',
    'sym_k': 'Eq(f(x), (C1 + C2*x)*exp(-k*x))',
    'quint': 'Eq(f(x), C5*exp(x*CRootOf(_x**5 + 10*_x - 2, 0)) + (C1*sin(x*im(CRootOf(_x**5 + 10*_x - 2, 1)))'
             ' + C2*cos(x*im(CRootOf(_x**5 + 10*_x - 2, 1))))*exp(x*re(CRootOf(_x**5 + 10*_x - 2, 1)))'
             ' + (C3*sin(x*im(CRootOf(_x**5 + 10*_x - 2, 3))) + C4*cos(x*im(CRootOf(_x**5 + 10*_x - 2, 3))))'
             '*exp(x*re(CRootOf(_x**5 + 10*_x - 2, 3))))',
    'quint_ns': 'Eq(f(x), C1*exp(x*CRootOf(_x**5 + 10*_x - 2, 0)) - C2*exp(x*re(CRootOf(_x**5 + 10*_x - 2, 1)))'
                '*sin(x*im(CRootOf(_x**5 + 10*_x - 2, 1))) + C3*exp(x*re(CRootOf(_x**5 + 10*_x - 2, 1)))'
                '*cos(x*im(CRootOf(_x**5 + 10*_x - 2, 1))) - C4*exp(x*re(CRootOf(_x**5 + 10*_x - 2, 3)))'
                '*sin(x*im(CRootOf(_x**5 + 10*_x - 2, 3))) + C5*exp(x*re(CRootOf(_x**5 + 10*_x - 2, 3)))'
                '*cos(x*im(CRootOf(_x**5 + 10*_x - 2, 3))))',
    'cplx_rep': 'Eq(f(x), (C1 + C2*x)*sin(x) + (C3 + C4*x)*cos(x))',
    'cplx_rep_ns': 'Eq(f(x), C1*sin(x) + C2*cos(x) + C3*x*sin(x) + C4*x*cos(x))',
    'mixed': 'Eq(f(x), (C1 + C2*sin(2*x) + C3*cos(2*x))*exp(x))',
    'rat': 'Eq(f(x), C1*exp(x/3) + C2*exp(x/2))',
    'd1': 'Eq(f(x), C1*exp(3*x))',
    'eu_cplx': 'Eq(f(x), x**2*(C1*sin(3*log(x)) + C2*cos(3*log(x))))',
    'eu_cplx_ns': 'Eq(f(x), x**2*(C1*sin(3*log(x)) + C2*cos(3*log(x))) + x**2*(C3*sin(3*log(x))'
                  ' + C4*cos(3*log(x))))',
    'eu_3': 'Eq(f(x), C1/x + x**2*(C2*sin(log(x)) + C3*cos(log(x))))',
    'eu_rep': 'Eq(f(x), x**2*(C1 + C2*log(x)))',
    'eu_rat': 'Eq(f(x), sqrt(x)*(C1 + C2*log(x)))',
    'eu_irr': 'Eq(f(x), C1*x**(3/2 - sqrt(17)/2) + C2*x**(3/2 + sqrt(17)/2))',
    'eu_ln': 'Eq(f(x), C1*sin(log(x)) + C2*cos(log(x)))',
    'und': 'Eq(f(x), C1*exp(-2*x) + C2*exp(-x) + exp(x)/6)',
    'und_sqrt': 'Eq(f(x), C1*exp(-sqrt(2)*x) + C2*exp(sqrt(2)*x) - x**2/2 - 1/2)',
    'vop': 'Eq(f(x), (C1 + x)*sin(x) + (C2 + log(cos(x)))*cos(x))',
    'vop_sqrt': 'Eq(f(x), C1*exp(-sqrt(3)*x) + C2*exp(sqrt(3)*x) + sqrt(3)*exp(x)/(6*(1 - sqrt(3)))'
                ' - sqrt(3)*exp(x)/(6*(1 + sqrt(3))))',
    'eu_und': 'Eq(f(x), x*(C1 + C2*x + x**2/2))',
    'eu_var': 'Eq(f(x), x*(C1 + C2*x + x**3/6))',
}


@pytest.mark.parametrize('name,eq,hint,simp', CASES, ids=[c[0] for c in CASES])
def test_dsolve_matches_sympy(name, eq, hint, simp):
    assert str(S.dsolve(eq(), fx, hint=hint, simplify=simp)) == SYMPY[name]


def test_default_hint_matches_sympy():
    assert str(S.dsolve(d(2) - 5*fx, fx)) == SYMPY['d2m5s']
    assert str(S.dsolve(d(5) + 10*d(1) - 2*fx, fx)) == SYMPY['quint']
    assert str(S.dsolve(x**2*d(2) - 3*x*d(1) + 13*fx, fx)) == SYMPY['eu_cplx']


def _stock(monkeypatch):
    monkeypatch.setattr(nat, 'const_characteristic_sols', lambda *args: None)
    monkeypatch.setattr(nat, 'simplified_sol', lambda *args: None)
    monkeypatch.setattr(nat, 'euler_characteristic_sols', lambda *args: None)


def _r(*coeffs):
    r = {i: S.sympify(c) for i, c in enumerate(coeffs)}
    r[-1] = S.Integer(0)
    return r


CHAR = [(-5, 0, 1), (1, 0, 1), (4, 4, 1), (-12, 0, 16, 0, -7, 0, 1), (1, 0, 2, 0, 1), (-5, 7, -3, 1),
        (1, -5, 6), (-6, -14, -1, 1), (-2, 10, 0, 0, 0, 1), (a, 0, 1), (-a**2, 0, 1), (k**2, 2*k, 1),
        (0, 0, 1), (8, 0, 0, 1), (4, -20, 41, -44, 26, -8, 1)]


@pytest.mark.parametrize('coeffs', CHAR, ids=[str(c) for c in CHAR])
def test_native_helpers_equal_reference_bodies(monkeypatch, coeffs):
    # The kernel reproduces the reference helper body run in symox up to the
    # order of the roots, which the kernel takes from SymPy's `ordered` (stock
    # symox orders `roots(multiple=True)` by sort key only).
    order = len(coeffs) - 1
    native = nat.const_characteristic_sols(_r(*coeffs), fx, order)
    assert native is not None
    _stock(monkeypatch)
    stock = nh._get_const_characteristic_eq_sols(_r(*coeffs), fx, order)
    # str: the RootOf polynomials are in a fresh dummy at each call
    assert sorted(map(str, native[0])) == sorted(map(str, stock[0]))
    assert sorted(map(str, native[1])) == sorted(map(str, stock[1]))


def test_const_kernel_values():
    gensols, collectterms = nh._get_const_characteristic_eq_sols(_r(1, 0, 1), fx, 2)
    assert gensols == [S.sin(x), S.cos(x)]
    # SymPy 1.14: roots [I, -I] in `ordered` order, terms prepended
    assert collectterms == [(0, 0, -1), (0, 0, 1)]
    gensols, collectterms = nh._get_const_characteristic_eq_sols(_r(-5, 0, 1), fx, 2)
    assert gensols == [S.exp(S.sqrt(5)*x), S.exp(-S.sqrt(5)*x)]
    assert collectterms == [(0, -S.sqrt(5), 0), (0, S.sqrt(5), 0)]


@pytest.mark.parametrize('coeffs', [
    (S.Float(2.0), 0, 1),       # inexact coefficients
    (S.I, 0, 1),                # complex coefficients
    (S.sqrt(2), 0, 1),          # algebraic coefficients
    (S.pi, 0, 1),               # transcendental coefficients
    (S.Function('g')(x), 0, 1),
])
def test_const_kernel_declines(coeffs):
    assert nat.const_characteristic_sols(_r(*coeffs), fx, 2) is None


def test_const_kernel_declines_bad_inputs():
    assert nat.const_characteristic_sols(_r(1, 0, 1), f(x, a), 2) is None
    assert nat.const_characteristic_sols(_r(1, 0, 1), f(x**2), 2) is None
    assert nat.simplified_sol([S.Eq(fx, x), S.Eq(fx, 1)], fx, []) is None
    assert nat.simplified_sol([S.Eq(f(a), x)], fx, []) is None


def test_euler_kernel_declines_symbolic_coefficients():
    from symox.solvers.ode.ode import get_numbered_constants
    eq = x**2*d(2) + a*fx
    match = {2: x**2, 0: a, -1: S.Integer(0)}
    assert nat.euler_characteristic_sols(eq, fx, match, get_numbered_constants) is None


def test_kernel_raw_bindings():
    assert _core.ode_const_char_sols([(0, S.Float(1.5)), (2, S.Integer(1))], x, 2) is None
    gensols, collectterms = _core.ode_const_char_sols([(0, S.Integer(-1)), (1, S.Integer(0)), (2, S.Integer(1))], x, 2)
    assert gensols == [S.exp(-x), S.exp(x)]
