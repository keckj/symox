"""Integration-corpus residue (wave 10); expected values from SymPy 1.14."""
from symox import *
del test
from symox.solvers.inequalities import solve_univariate_inequality as sui
from symox.integrals.meijerint import meijerint_definite

x, t, s, a = symbols('x t s a')


def _ineqs(k):
    return [
        (k**2 > 4, Union(Interval.open(-oo, -2), Interval.open(2, oo))),
        (k**2 > pi, Union(Interval.open(-oo, -sqrt(pi)), Interval.open(sqrt(pi), oo))),
        (k**2 > 1/(4*pi**2),
         Union(Interval.open(-oo, -1/(2*pi)), Interval.open(1/(2*pi), oo))),
        (k**2 < pi, Interval.open(-sqrt(pi), sqrt(pi))),
        ((k - 1)*(k + sqrt(2)) > 0,
         Union(Interval.open(-oo, -sqrt(2)), Interval.open(1, oo))),
        (k**2 >= E, Union(Interval(-oo, -exp(S.Half)), Interval(exp(S.Half), oo))),
        (1/(k**2 - pi) > 0,
         Union(Interval.open(-oo, -sqrt(pi)), Interval.open(sqrt(pi), oo))),
    ]


def test_univariate_inequality_ignores_sign_assumptions():
    # SymPy keeps the real generator but finds roots over the reals: the sign
    # assumptions of the symbol do not restrict the returned set.
    for assumptions in ({}, {'real': True}, {'positive': True}, {'negative': True}):
        k = Symbol('k', **assumptions)
        for ineq, expected in _ineqs(k):
            assert sui(ineq, k, relational=False) == expected, (assumptions, ineq)
    k = Symbol('k', positive=True)
    assert sui(k**3 > 2, k, relational=False) == Interval.open(2**Rational(1, 3), oo)


def test_fourier_sinc_positive_frequency():
    k = Symbol('k', positive=True)
    assert fourier_transform(sin(x)/x, x, k) == \
        Piecewise((0, k**2 > 1/(4*pi**2)), (pi, True))
    assert reduce_inequalities(k**2 > 1/(4*pi**2), k) == (1/(2*pi) < k)


def test_laplace_power_open_condition_uses_integration():
    assert laplace_transform(t**a, t, s) == (gamma(a + 1)/(s*s**a), 0, re(a) > -1)
    assert laplace_transform(exp(-t)*t**a, t, s) == \
        (gamma(a + 1)/((s + 1)*(s + 1)**a), -1, re(a) > -1)
    ar = Symbol('a', real=True)
    assert laplace_transform(t**ar, t, s) == (gamma(ar + 1)/(s*s**ar), 0, ar > -1)
    assert laplace_transform(sqrt(t), t, s, noconds=True) == sqrt(pi)/(2*s**Rational(3, 2))
    assert laplace_transform(t**3, t, s, noconds=True) == 6/s**4


def test_meijerint_definite_split_at_symbol_root():
    assert meijerint_definite(1/(x*(1 + x)), x, 1, oo) == (log(2), True)
    assert meijerint_definite(1/(x*(2 + x)), x, 1, oo) == (log(3)/2, True)
    assert meijerint_definite(1/(x**2*(1 + x)), x, 1, oo) == (1 - log(2), True)
    assert meijerint_definite(1/(x*(1 + x)**2), x, 1, oo) == (-S.Half + log(2), True)
    assert meijerint_definite(x**a, x, 1, oo) == (-1/(a + 1), re(a) < -1)
    assert integrate(1/(x*(1 + x)), (x, 1, oo)) == log(2)
    assert integrate(1/(x**2*(1 + x)), (x, 1, oo)) == 1 - log(2)
    assert integrate(x**a, (x, 1, oo)) == \
        Piecewise((-1/(a + 1), re(a) < -1), (Integral(x**a, (x, 1, oo)), True))
