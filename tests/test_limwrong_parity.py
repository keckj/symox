"""Limits that were wrong before wave 12 (limwrong); expected values are
SymPy 1.14's (sympy/series/tests/test_limits.py and test_gruntz.py)."""
import pytest

from symox import *
from symox.calculus.accumulationbounds import AccumBounds
from symox.series.gruntz import gruntz

del test

x, z = symbols('x z')


def test_log_cancellation_with_positive_symbols():
    # issue 21598: log(1/b) + log(b) is a hidden zero leading coefficient
    a, b, c = symbols('a b c', positive=True)
    A = log(a/b) - (log(a) - log(b))
    assert A.limit(a, oo) == 0
    assert (A*c).limit(a, oo) == 0


def test_log_manualintegrate_value():
    tau, x = symbols('tau x', positive=True)
    expr = tau**2*((tau - 1)*(tau + 1)*log(x + 1)/(tau**2 + 1)**2 + 1/((tau**2
            + 1)*(x + 1)) - (-2*tau*atan(x/tau) + (tau**2/2 - S(1)/2)*log(tau**2
            + x**2))/(tau**2 + 1)**2)
    assert limit(expr, x, oo) == pi*tau**3/(tau**2 + 1)**2


def test_frac_power_bounds():
    # AccumBounds(0, 1)**t with t -> 0+ keeps the zero endpoint
    assert limit(frac(x)**(1/x), x, oo) == AccumBounds(0, 1)
    assert limit(frac(x)**(1/x), x, -oo) == AccumBounds(1, oo)
    assert limit((tan(x)**2)**(2/x), x, oo) == AccumBounds(0, oo)
    assert limit((tan(x)**2)**(2/x), x, -oo) == AccumBounds(0, oo)


def test_issue_6366_and_6843():
    n = Symbol('n', integer=True, positive=True)
    r = (n + 1)*x**(n + 1)/(x**(n + 1) - 1) - x/(x - 1)
    assert limit(r, x, 1).cancel() == n/2
    assert gruntz(r, x, 1).simplify() == n/2


def test_issue_6364_hidden_trig_zero():
    a = Symbol('a')
    e = z/(1 - sqrt(1 + z)*sin(a)**2 - sqrt(1 - z)*cos(a)**2)
    r = limit(e, z, 0)
    # SymPy 1.14: 1/(cos(a)**2 - 1/2); symox keeps an equivalent form
    assert r != 0 and not r.has(Limit)
    assert simplify(r - 1/(cos(a)**2 - S.Half)) == 0


def test_issue_18399_zoo():
    assert limit((1 - S(1)/2*x)**(3*x), x, oo) is zoo
    assert limit((-x)**x, x, oo) is zoo


def test_issue_22220_branch_cut():
    e1 = sqrt(30)*atan(sqrt(30)*tan(x/2)/6)/30
    e2 = sqrt(30)*I*(-log(sqrt(2)*tan(x/2) - 2*sqrt(15)*I/5) +
                     +log(sqrt(2)*tan(x/2) + 2*sqrt(15)*I/5))/60
    assert limit(e1, x, -pi) == -sqrt(30)*pi/60
    assert limit(e2, x, -pi) == -sqrt(30)*pi/30
    assert limit(e1, x, -pi, '-') == sqrt(30)*pi/60
    assert limit(e2, x, -pi, '-') == 0
    expr = -log(tan(x/2) - I) + log(tan(x/2) + I)
    assert limit(expr, x, pi, '+') == 2*I*pi
    assert limit(expr, x, pi, '-') == 0


def test_acoth_branch_point():
    assert limit(acoth(exp(1/x)), x, 0, '+') == 0
    assert limit(acoth(exp(1/x)), x, 0, '-') == -I*pi/2


def test_issue_22334_not_wrong():
    k, n = symbols('k, n', positive=True)
    r = limit((n + 1)**k/((n + 1)**(k + 1) - n**(k + 1)), n, oo)
    # was 0; SymPy gives 1/(k + 1), symox now stays unevaluated
    assert r == 1/(k + 1) or isinstance(r, Limit)


@pytest.mark.xfail(strict=True, reason="symbolic powers of the mrv variable")
def test_issue_22334_value():
    k, n = symbols('k, n', positive=True)
    assert limit((n + 1)**k/((n + 1)**(k + 1) - n**(k + 1)), n, oo) == 1/(k + 1)
