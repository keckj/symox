"""Regressions found by running SymPy 1.14's own integration test files
against symox (wave 25, integsuite). Expected values are SymPy 1.14's,
except the intentional correction checked by an exact certificate."""
import mpmath
from symox import (Integral, Ne, Piecewise, Rational, Symbol, cos, diff, floor,
                   integrate, log, pi, sqrt, symbols)

x, y, n = symbols('x y n')
theta = Symbol('theta')


def test_issue_21831_definite_values():
    # the manualintegrate stage does its inner integrals with final=False:
    # the floor terms of atan(c*tan(a)) are added once, in the final round
    assert integrate(cos(3*theta)/(5 - 4*cos(theta)), (theta, 0, 2*pi)) == pi/12
    assert integrate(cos(2*theta)/(5 - 4*cos(theta)), (theta, 0, 2*pi)) == pi/6
    assert integrate(cos(theta)/(5 - 4*cos(theta)), (theta, 0, 2*pi)) == pi/3


def test_issue_21831_antiderivative_continuous():
    f = cos(2*theta)/(5 - 4*cos(theta))
    F = integrate(f, theta)
    assert F.count(floor) == 3  # as SymPy 1.14: one per atan occurrence
    with mpmath.workdps(30):
        eps = Rational(1, 10**12)
        left, right = (F.subs(theta, v).evalf(30) for v in (pi - eps, pi + eps))
        assert abs(left - right) < 1e-9
        for v in (Rational(3, 10), Rational(17, 10), Rational(21, 5)):
            assert abs((diff(F, theta) - f).subs(theta, v).evalf(30)) < 1e-20


def test_doit_final_false_has_no_floor_terms():
    i = Integral(1/(2 + cos(x)), x)
    assert i.doit(final=False).count(floor) == 0
    assert i.doit().count(floor) == 1
    assert integrate(1/(2 + cos(x)), (x, 0, 2*pi)) == 2*sqrt(3)*pi/3


def test_exponential_symbolic_slope_guard():
    assert integrate(x**(n*y), y) == Piecewise(
        (x**(n*y)/(n*log(x)), Ne(n*log(x), 0)), (y, True))
    assert integrate(2**(n*y), y) == Piecewise(
        (2**(n*y)/(n*log(2)), Ne(n*log(2), 0)), (y, True))
    # SymPy 1.14 distributes the 3 into the branches; the same value
    assert integrate(3*x**(n*y), y) == 3*Piecewise(
        (x**(n*y)/(n*log(x)), Ne(n*log(x), 0)), (y, True))
    assert integrate(x**(2*y + 1), y) == Piecewise(
        (x**(2*y + 1)/(2*log(x)), Ne(log(x), 0)), (y, True))
    assert integrate(x**(y/n), y) == Piecewise(
        (n*x**(y/n)/log(x), Ne(log(x), 0)), (y, True))
    r = integrate(x**(-n*y), y)
    assert r == Piecewise((-1/(n*x**(n*y)*log(x)), Ne(n*log(x), 0)), (y, True))
    assert integrate(2**(3*y), y) == 2**(3*y)/(3*log(2))
    assert integrate(3*2**(n*y), y) == 3*Piecewise(
        (2**(n*y)/(n*log(2)), Ne(n*log(2), 0)), (y, True))
    # the slope is a multiple of the coefficient: SymPy also drops the guard
    assert integrate(n*2**(n*y), y) == 2**(n*y)/log(2)


def test_coefficient_not_multiple_of_guard_keeps_branch():
    # SymPy 1.14: Piecewise((x**(n*y)/log(x), Ne(log(x), 0)), (n*y, True));
    # symox keeps the equivalent n*Piecewise(..., Ne(n*log(x), 0)) form, not
    # the former x**(n*y)/log(x), which is undefined at x = 1
    r = integrate(n*x**(n*y), y)
    assert r.subs(x, 1) == n*y
    assert diff(r, y).subs({x: 2, n: 5}) == 5*2**(5*y)
    assert r.subs(n, 0) == 0


def test_exponential_symbolic_slope_offset_correction():
    # intentional correction: SymPy 1.14 gives (y, True) as the second branch,
    # but where n = 0 the integrand is x**3, whose antiderivative is x**3*y
    r = integrate(x**(n*y + 3), y)
    assert r == Piecewise((x**(n*y + 3)/(n*log(x)), Ne(n*log(x), 0)), (x**3*y, True))
    assert diff(r.args[0][0], y) == x**(n*y + 3)
    assert diff(r.args[1][0], y) == (x**(n*y + 3)).subs(n, 0)
    assert r.subs(x, 1) == y
