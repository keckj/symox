"""Differential holonomic checks; expected values hard-coded from SymPy 1.14."""
import pytest
from symox import *
from symox.holonomic import (HolonomicFunction, DifferentialOperators,
                             expr_to_holonomic, from_hyper)
from symox.holonomic.holonomicerrors import NotHolonomicError
del test

x, a, n = symbols('x a n')

# Each pair is (symox expression, str of the SymPy 1.14 result).
REFERENCE = [
    ('expr_to_holonomic((x+1)**n,x)', 'HolonomicFunction((-n) + (x + 1)*Dx, x, 0, [1])'),
    ('expr_to_holonomic(x**-2,x)', 'HolonomicFunction((2) + (x)*Dx, x, 0, {-2: [1]})'),
    ('expr_to_holonomic(exp(x)**-1,x)', 'HolonomicFunction((1) + (1)*Dx, x, 0, [1])'),
    ('expr_to_holonomic(sin(x),x)**3', 'HolonomicFunction((9) + (10)*Dx**2 + (1)*Dx**4, x, 0, [0, 0, 0, 6])'),
    ('expr_to_holonomic(sin(x),x)**0', 'HolonomicFunction((1)*Dx, x, 0, [1])'),
    ('expr_to_holonomic(sin(x),x,x0=1)**2', 'HolonomicFunction((4)*Dx + (1)*Dx**3, x, 0, [0, 0, 2])'),
    ('expr_to_holonomic(erfc(x),x)', 'HolonomicFunction((2*x)*Dx + (1)*Dx**2, x, 0, [1, -2/sqrt(pi)])'),
    ('expr_to_holonomic(erfc(2*x),x)', 'HolonomicFunction((8*x)*Dx + (1)*Dx**2, x, 0, [1, -4/sqrt(pi)])'),
    ('expr_to_holonomic(erfc(x)**2,x)', 'HolonomicFunction((8*x**2 + 2)*Dx + (6*x)*Dx**2 + (1)*Dx**3, x, 0, [1, -4/sqrt(pi), 8/pi])'),
    ('expr_to_holonomic(1.2*cos(3.1*x),x)', 'HolonomicFunction((9.61000000000000) + (1.00000000000000)*Dx**2, x, 0, [1.2, 0])'),
    ('expr_to_holonomic(1.2*cos(3.1*x),x).to_expr()', '1.2*cos(3.1*x)'),
    ('expr_to_holonomic(1.2*cos(3.1*x),x).integrate(x).to_expr()', '0.387096774193548*sin(3.1*x)'),
    ('expr_to_holonomic(1.1329138213*x,x).integrate((x,1,2))', '1.69937073195000'),
    ('expr_to_holonomic(0.5*exp(x),x)', 'HolonomicFunction((-1.00000000000000) + (1.00000000000000)*Dx, x, 0, [0.5])'),
    ('expr_to_holonomic(x,x)+expr_to_holonomic(1.2*cos(x),x)', 'HolonomicFunction((-1.00000000000000) + (1.0*x)*Dx + (-1.00000000000000)*Dx**2 + (1.0*x)*Dx**3, x, 0, {0: [1.2, 0], 1: [1]})'),
    ('expr_to_holonomic(x,x)*expr_to_holonomic(1.2*cos(x),x)', 'HolonomicFunction((1.0*x**2 + 2.0) + (-2.0*x)*Dx + (1.0*x**2)*Dx**2, x, 0, {1: [1.2]})'),
    ('expr_to_holonomic(1.5*a*sin(x),x)', 'HolonomicFunction((1.00000000000000) + (1.00000000000000)*Dx**2, x, 0, [0, 1.5*a])'),
    ('expr_to_holonomic(x**2+x,x)', 'HolonomicFunction((-2*x - 1) + (x**2 + x)*Dx, x, 0, {1: [1, 1]})'),
    ('expr_to_holonomic(x,x)', 'HolonomicFunction((-1) + (x)*Dx, x, 0, {1: [1]})'),
    ('expr_to_holonomic(cos(x),x)+expr_to_holonomic(x,x)', 'HolonomicFunction((-1) + (x)*Dx + (-1)*Dx**2 + (x)*Dx**3, x, 0, {0: [1, 0], 1: [1]})'),
    ('expr_to_holonomic(sin(x)/x,x)', 'HolonomicFunction((x) + (2)*Dx + (x)*Dx**2, x, 0, [1, 0])'),
    ('expr_to_holonomic(log(x),x)', 'HolonomicFunction((1)*Dx + (x)*Dx**2, x, 1, [0, 1])'),
    ('erf(x).rewrite(erfc)', '1 - erfc(x)'),
    ('erfi(x).rewrite(erfc)', 'I*erfc(I*x) - I'),
    ('erf2(a,x).rewrite(erfc)', 'erfc(a) - erfc(x)'),
    ('erfc(x).rewrite(erfi)', 'I*erfi(I*x) + 1'),
    ('erf(x).rewrite(uppergamma)', '(1 - erfc(sqrt(x**2)))*sqrt(x**2)/x'),
    ('erfc(x).rewrite(expint)', 'x*erfc(sqrt(x**2))/sqrt(x**2) + 1 - sqrt(x**2)/x'),
    ('erf(x).rewrite(hyper)', '2*x*hyper((1/2,), (3/2,), -x**2)/sqrt(pi)'),
    ('erfc(x).rewrite(fresnels)', '-(1 + I)*(fresnelc(x*(1 - I)/sqrt(pi)) - I*fresnels(x*(1 - I)/sqrt(pi))) + 1'),
    ('erfinv(x).rewrite(erfcinv)', 'erfcinv(1 - x)'),
    ('erfcinv(x).rewrite(erfinv)', 'erfinv(1 - x)'),
    ('fresnels(x).rewrite(erf)', '(1/4 + I/4)*(-I*erf(sqrt(pi)*x*(1/2 - I/2)) + erf(sqrt(pi)*x*(1/2 + I/2)))'),
    ('fresnelc(x).rewrite(hyper)', 'x*hyper((1/4,), (1/2, 5/4), -pi**2*x**4/16)'),
    ('hyperexpand(expr_to_holonomic(erfc(x),x).to_hyper()).rewrite(erfc)', 'erfc(x)'),
]

# SymPy normalizes these annihilators by a different constant factor
# (documented divergence); the initial conditions are identical.
SCALED = [
    ('expr_to_holonomic(2.5*x**2,x)', 'HolonomicFunction((-5.00000000000000) + (2.5*x)*Dx, x, 0, {2: [2.5]})'),
    ('expr_to_holonomic(x,x)*expr_to_holonomic(sqrt(x),x)', 'HolonomicFunction((-3/2) + (x)*Dx, x, 0, {3/2: [1]})'),
]

ERRORS = [
    ('expr_to_holonomic(sin(x)**-1,x)', NotHolonomicError, 'Negative Power on a Holonomic Function'),
    ('expr_to_holonomic(1/erf(x),x)', NotHolonomicError, 'Negative Power on a Holonomic Function'),
    ('expr_to_holonomic(sqrt(sin(x)),x)', TypeError, "unsupported operand type(s) for >>=: 'Half' and 'int'"),
    ('expr_to_holonomic(exp(x),x)**n', CoercionFailed, 'expected `Rational` object, got -n'),
    ('expr_to_holonomic(sin(x),x)**n', TypeError, 'cannot determine truth value of Relational: n < 0'),
]


def _evaluate(source):
    return eval(source, globals())


@pytest.mark.parametrize('source, expected', REFERENCE)
def test_matches_sympy(source, expected):
    assert str(_evaluate(source)) == expected


@pytest.mark.parametrize('source, expected', SCALED)
def test_matches_sympy_up_to_scale(source, expected):
    result = _evaluate(source)
    reference = expected.rsplit(', x, ', 1)
    assert str(result).rsplit(', x, ', 1)[1] == reference[1]
    ours = result.annihilator.coefficients
    theirs = [sympify(c) for c in _sympy_coefficients(reference[0])]
    ratio = cancel(theirs[-1]/ours[-1])
    assert ratio.is_number and ratio != 0
    assert [expand(ratio*c - t) for c, t in zip(ours, theirs)] == [0]*len(ours)


def _sympy_coefficients(text):
    # 'HolonomicFunction((c0) + (c1)*Dx' -> ['c0', 'c1'] for first order.
    body = text[len('HolonomicFunction('):]
    c0, c1 = body.split(' + (', 1)
    return [c0.strip('()'), c1.split(')*Dx')[0]]


@pytest.mark.parametrize('source, error, message', ERRORS)
def test_errors_match_sympy(source, error, message):
    with pytest.raises(error) as info:
        _evaluate(source)
    if error is not CoercionFailed:
        assert str(info.value) == message


def test_power_initial_conditions_of_singular_functions():
    # Intentional correction: SymPy 1.14 raises y0 keys to the power
    # (giving [1/8]); the Frobenius conditions of (3*sqrt(x))**3 are {3/2: [27]}.
    h = expr_to_holonomic(3*sqrt(x), x)**3
    assert h.y0 == {Rational(3, 2): [27]}
    assert expr_to_holonomic(x**n, x).y0 == {n: [1]}


def test_rational_power_type_error():
    with pytest.raises(TypeError, match="unsupported operand type"):
        expr_to_holonomic(sin(x), x)**Rational(1, 2)


def test_float_domain_numerics():
    h = expr_to_holonomic(1.2*cos(3.1*x), x)
    assert h.annihilator.parent.base.domain == RR
    # SymPy 1.14 RK4 (h=0.05): 0.0249626374429972
    assert abs(h.evalf(0.5) - Float('0.0249626374429972')) < 1e-13
    # SymPy 1.14: 1.69937073195000
    assert abs(expr_to_holonomic(1.1329138213*x, x).integrate((x, 1, 2)) - Float('1.69937073195')) < 1e-12


def test_extended_domain_matches_sympy_test():
    p = expr_to_holonomic(1.1329138213*x)
    q = expr_to_holonomic(1.1329138213*x, domain=RR)
    assert p == q
    assert p.annihilator.parent.base.domain == RR


def test_sum_and_product_move_conditions_to_zero():
    s1 = expr_to_holonomic(sin(x), x, x0=1)
    # SymPy 1.14: HolonomicFunction((1) + (1)*Dx**2, x, 0, [1, cos(1)**2 + sin(1)**2])
    h = s1 + expr_to_holonomic(cos(x), x)
    assert str(h.annihilator) == '(1) + (1)*Dx**2' and h.x0 == 0
    assert h.y0[0] == 1 and simplify(h.y0[1] - 1) == 0
    assert str(expr_to_holonomic(sin(x), x, x0=1)**2) == 'HolonomicFunction((4)*Dx + (1)*Dx**3, x, 0, [0, 0, 2])'


def test_symox_extension_product_with_polynomial_keeps_jet():
    R, D = DifferentialOperators(QQ.old_poly_ring(x), 'D')
    sine = HolonomicFunction(D**2 + 1, x, 0, [0, 1])
    assert (sine*x).y0 == [0, 0, 2]


def test_integrate_from_removable_singularity():
    # SymPy 1.14: HolonomicFunction((4)*Dx + (4*x)*Dx**2 + (3)*Dx**3 + (x)*Dx**4,
    # x, 0, [0, 0, 1, 0]) and log(x)/2 - Ci(2*x)/2 + EulerGamma/2 + log(2)/2
    # (the base returned nan here).
    h = expr_to_holonomic(sin(x)**2/x).integrate((x, 0, x))
    assert str(h) == 'HolonomicFunction((4)*Dx + (4*x)*Dx**2 + (3)*Dx**3 + (x)*Dx**4, x, 0, [0, 0, 1, 0])'
    assert h.to_expr() == log(x)/2 - Ci(2*x)/2 + EulerGamma/2 + log(2)/2
    # SymPy 1.14: -Ci(4)/2 + EulerGamma/2 + log(2) (quad: 1.05224586195418)
    value = expr_to_holonomic(sin(x)**2/x).integrate((x, 0, 2))
    assert value == -Ci(4)/2 + EulerGamma/2 + log(2)
    assert abs(value.evalf() - Float('1.05224586195418')) < 1e-13
    # SymPy 1.14: log(x)/2 - Ci(2*x)/2 + Ci(2)/2
    h = expr_to_holonomic(sin(x)**2/x, x0=1).integrate((x, 1, x))
    assert h.to_expr() == log(x)/2 - Ci(2*x)/2 + Ci(2)/2
    assert not h.to_expr().has(nan)


def test_erfc_to_expr_matches_sympy():
    # SymPy 1.14: 1 - erf(x), 1 - erf(2*x)
    assert expr_to_holonomic(erfc(x), x).to_expr() == 1 - erf(x)
    assert expr_to_holonomic(erfc(2*x), x).to_expr() == 1 - erf(2*x)


def test_power_correction_series():
    # Intentional correction: f = (sqrt(x)*exp(x))**3 vanishes at 0, so
    # SymPy's ordinary [1/8] is wrong; f/x**(3/2) = 1 + 3*x + ... gives {3/2: [1]}.
    h = expr_to_holonomic(sqrt(x)*exp(x), x)**3
    assert h.y0 == {Rational(3, 2): [1]}
    assert series(exp(3*x), x, 0, 1).removeO() == 1
