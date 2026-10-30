"""Native second-order classify_ode matches (2nd_linear_airy,
2nd_linear_bessel rejection, match_2nd_hypergeometric, order-two power
series) against SymPy 1.14 and against the reference matchers."""
import pytest

import symox as S
from symox import _core
from symox.solvers.ode import _classify_native as cn
from symox.solvers.ode import hypergeometric as hg
from symox.solvers.ode import single as s
from symox.solvers.ode.single import SingleODEProblem

x, y = S.symbols('x y')
f = S.Function('f')
fx = f(x)
d1 = fx.diff(x)
d2 = fx.diff(x, 2)
ENV = {'x': x, 'y': y, 'fx': fx, 'd1': d1, 'd2': d2, 'Float': S.Float, 'exp': S.exp, 'sin': S.sin}


def _eq(src):
    return eval(src, dict(ENV))


# SymPy 1.14: classify_ode(eq, f(x)) and match_2nd_hypergeometric(eq, f(x)).
EXPECTED = [
    ('d2 - x*fx', ('2nd_linear_airy', '2nd_power_series_ordinary'), ['0', '-x']),
    ('d2 + (2 + 3*x)*fx', ('2nd_linear_airy', '2nd_power_series_ordinary'), ['0', '3*x + 2']),
    ('x**2*d2 + x*d1 + (x**2 - 4)*fx', ('factorable', '2nd_linear_bessel', '2nd_power_series_regular'),
     ['1/x', '(x**2 - 4)/x**2']),
    ('x*(1 - x)*d2 + (1 - 3*x)*d1 - fx', ('factorable', '2nd_power_series_regular'),
     ['(3*x - 1)/(x**2 - x)', '1/(x**2 - x)']),
    ('d2 + x*d1 + fx', ('factorable', '2nd_power_series_ordinary'), ['x', '1']),
    ('d2 + fx/x', ('2nd_linear_bessel', '2nd_power_series_regular'), ['0', '1/x']),
    ('x**2*d2 - x*d1 + fx', ('factorable', 'nth_linear_euler_eq_homogeneous', '2nd_linear_bessel',
                             '2nd_power_series_regular'), ['-1/x', 'x**(-2)']),
    ('2*d2 + 3*fx', ('nth_linear_constant_coeff_homogeneous', '2nd_power_series_ordinary',
                     '2nd_nonlinear_autonomous_conserved', '2nd_nonlinear_autonomous_conserved_Integral'), ['0', '3/2']),
    ('d2 + d1/x + fx/x**2', ('factorable', 'nth_linear_euler_eq_homogeneous', '2nd_linear_bessel',
                             '2nd_power_series_regular'), ['1/x', 'x**(-2)']),
    ('(1 - x**2)*d2 - 2*x*d1 + 6*fx', ('factorable', '2nd_hypergeometric', '2nd_hypergeometric_Integral',
                                       '2nd_power_series_ordinary'), ['2*x/(x**2 - 1)', '-6/(x**2 - 1)']),
    ('4*x**2*d2 + 8*x*d1 + (x - 3)*fx', ('factorable', '2nd_linear_bessel', '2nd_power_series_regular'),
     ['2/x', '(x - 3)/(4*x**2)']),
    ('d2 - y*x*fx', ('2nd_linear_airy', '2nd_power_series_ordinary'), ['0', '-x*y']),
]

# SymPy 1.14: dsolve(eq, f(x), hint=hint).
DSOLVE = [
    ('d2 - x*fx', '2nd_linear_airy', 'Eq(f(x), C1*airyai(x) + C2*airybi(x))'),
    ('d2 + (2 + 3*x)*fx', '2nd_linear_airy',
     'Eq(f(x), C1*airyai(3**(1/3)*(-x - 2/3)) + C2*airybi(3**(1/3)*(-x - 2/3)))'),
    ('x**2*d2 + x*d1 + (x**2 - 4)*fx', '2nd_linear_bessel', 'Eq(f(x), C1*besselj(2, x) + C2*bessely(2, x))'),
    ('d2 + x*d1 + fx', '2nd_power_series_ordinary',
     'Eq(f(x), C2*(x**4/8 - x**2/2 + 1) + C1*x*(1 - x**2/3) + O(x**6))'),
]

NATIVE = (cn.airy_match, cn.bessel_rejected, cn.match_2nd_hypergeometric, cn.power_series2)


@pytest.fixture
def stock(monkeypatch):
    def off():
        monkeypatch.setattr(cn, 'airy_match', lambda solver, eq: None)
        monkeypatch.setattr(cn, 'bessel_rejected', lambda solver, eq: False)
        monkeypatch.setattr(cn, 'match_2nd_hypergeometric', lambda eq, func: cn._MISSING)
        monkeypatch.setattr(cn, 'power_series2', lambda *args: cn._MISSING)
    return off


def test_bindings_present():
    for name in ('ode_match_2nd_hypergeometric', 'ode_airy_match', 'ode_bessel_rejected', 'ode_power_series2'):
        assert hasattr(_core, name)


@pytest.mark.parametrize('src,hints,ab', EXPECTED)
def test_classify_ode_matches_sympy(src, hints, ab):
    eq = _eq(src)
    assert S.classify_ode(eq, fx) == hints
    assert [str(v) for v in hg.match_2nd_hypergeometric(eq, fx)] == ab


@pytest.mark.parametrize('src,hint,want', DSOLVE)
def test_dsolve_matches_sympy(src, hint, want):
    assert str(S.dsolve(_eq(src), fx, hint=hint)) == want


def _decisions(eq):
    problem = SingleODEProblem(eq, fx, x)
    airy = s.SecondLinearAiry(problem)
    ok = airy._matches()
    bessel = s.SecondLinearBessel(SingleODEProblem(eq, fx, x))._matches()
    out = [ok, airy.rn if ok else None, bessel, hg.match_2nd_hypergeometric(problem.eq_preprocessed, fx)]
    for x0 in (0, 1, y):
        d = S.classify_ode(eq, fx, dict=True, x0=x0)
        out.append([k for k in d if k not in ('order', 'default')])
        out.append({k: v for k, v in d.items() if str(k).startswith('2nd_power_series')})
    return out


EDGE = [src for src, _, _ in EXPECTED] + [
    'd2 + Float(1.5)*x*fx', 'Float(0.0)*d2 + d1 + fx', 'x*d2 + Float(0.0)*fx', 'exp(x)*d2 + fx',
    'd2 + sin(x)*fx', 'd2 - (x + 1)*fx - 1', 'd2 + x**2*fx', '3*d2 - 6*x*fx', 'y*d2 + fx', 'd1 + fx', 'x**3*d2 + fx']


@pytest.mark.parametrize('src', EDGE)
def test_native_equals_reference(src, stock):
    eq = _eq(src)
    native = _decisions(eq)
    stock()
    assert _decisions(eq) == native


def test_native_decides_the_textbook_cases():
    def problem(src):
        return SingleODEProblem(_eq(src), fx, x)
    p = problem('d2 + (2 + 3*x)*fx')
    assert cn.airy_match(s.SecondLinearAiry(p), p.eq_high_order_free) == (True, {'b': 2, 'm': 3})
    p = problem('d2 + x**2*fx')
    assert cn.airy_match(s.SecondLinearAiry(p), p.eq_high_order_free) == (False, None)
    p = problem('2*d2 + 3*fx')
    assert cn.bessel_rejected(s.SecondLinearBessel(p), p.eq_high_order_free) is True
    p = problem('x**2*d2 + x*d1 + (x**2 - 4)*fx')
    # x-dependent leading coefficient: the factorization step stays with the reference
    assert cn.bessel_rejected(s.SecondLinearBessel(p), p.eq_high_order_free) is False
    assert [str(v) for v in cn.match_2nd_hypergeometric(_eq('d2 + x*d1 + fx'), fx)] == ['x', '1']
    assert cn.match_2nd_hypergeometric(_eq('d1 + fx'), fx) == []
    a3, b3, c3 = (S.Wild(n) for n in ('a3', 'b3', 'c3'))
    hints = cn.power_series2(_eq('d2 + x*d1 + fx'), f, x, 0, (a3, b3, c3), 5)
    assert list(hints) == ['2nd_power_series_ordinary']
    hints = cn.power_series2(_eq('x**2*d2 + x*d1 + (x**2 - 4)*fx'), f, x, 0, (a3, b3, c3), 5)
    assert list(hints) == ['2nd_power_series_regular']
    assert hints['2nd_power_series_regular']['p'] == 1


def test_declines():
    # non-polynomial wild values: the reference matches the expanded numerator again
    assert cn.match_2nd_hypergeometric(_eq('d2 + d1/x + fx/x**2'), fx) is cn._MISSING
    # a float leading coefficient: the float comparison `!= 0` is not modelled
    assert cn.match_2nd_hypergeometric(_eq('Float(1.5)*d2 + d1 + fx'), fx) is cn._MISSING
    # not a core expression
    assert cn.match_2nd_hypergeometric(S.Matrix([1]), fx) is cn._MISSING
