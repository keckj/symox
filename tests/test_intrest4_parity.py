"""Regression tests of task intrest4 (wave 13): integration corpus leftovers.

Expected values are SymPy 1.14's, except the intentional corrections, whose
values are checked numerically against quadrature.
"""
import mpmath

from symox import *
from symox.integrals.meijerint import meijerint_indefinite

del test

x, t, s = symbols('x t s')
a = Symbol('a', positive=True)
b = Symbol('b')
r = Symbol('r', real=True)


def _quad_heaviside(shift, sign, lo, hi):
    # integral of Heaviside(sign*(x - shift)) over (lo, hi)
    pts = sorted({lo, hi, min(max(shift, lo), hi)})
    return float(mpmath.quad(lambda u: 1 if sign*(u - shift) > 0 else 0, pts))


# Intentional corrections: the reference's Meijer G antiderivatives of a
# shifted Heaviside are valid for x > 0 only (SymPy gives 2 for the first).

def test_heaviside_shift_numeric_interval_crossing_zero():
    assert integrate(Heaviside(x - S(1)/2), (x, -1, 1)) == S(1)/2
    assert integrate(x*Heaviside(x - S(1)/2), (x, -1, 1)) == S(3)/8


def test_heaviside_shift_positive_parameter():
    v = integrate(Heaviside(x - a), (x, -1, 1))
    assert v == (1 - a)*Heaviside(1 - a)
    for f, sign, shift in [(Heaviside(a - x), -1, 1), (Heaviside(x + a), 1, -1)]:
        v = integrate(f, (x, -1, 1))
        assert not v.has(Integral, meijerg)
        for val in [S(3)/10, S(7)/10, 2, 3]:
            assert abs(float(v.subs(a, val)) - _quad_heaviside(shift*float(val), sign, -1, 1)) < 1e-9


def test_heaviside_shift_unassumed_and_real_parameter():
    for p in [b, r]:
        v = integrate(Heaviside(x - p), (x, -1, 1))
        assert v == (1 - p)*Heaviside(1 - p) - (-p - 1)*Heaviside(-p - 1)
        for val in [-3, -S(1)/2, S(1)/2, 3]:
            assert abs(float(v.subs(p, val)) - _quad_heaviside(float(val), 1, -1, 1)) < 1e-9


def test_heaviside_shift_indefinite():
    assert integrate(Heaviside(x - S(1)/2), x) == (x - S(1)/2)*Heaviside(x - S(1)/2)
    assert integrate(Heaviside(x - a), x) == (x - a)*Heaviside(x - a)
    assert str(integrate(Heaviside(2*x + 4), x)) == '(x + 2)*Heaviside(2*x + 4)'
    # the public Meijer G function keeps the reference's result
    assert str(meijerint_indefinite(Heaviside(x - S(1)/2), x)) == (
        'Piecewise((0, Abs(x) < 1/2), (x*meijerg(((0, 1), ()), ((), (-1, 0)), 2*x), True))')


def test_heaviside_shift_nonnegative_interval_keeps_reference_form():
    # on (0, 1) the Meijer G form is valid and is the reference's
    assert str(integrate(Heaviside(x - a), (x, 0, 1))) == (
        'meijerg(((0, 1), ()), ((), (-1, 0)), 1/a) - '
        'meijerg(((0, 1), ()), ((), (-1, 0)), Min(1, a)/a)*Min(1, a)')
    assert str(integrate(Heaviside(x - b), (x, 0, 1))) == (
        'Piecewise((Integral(Piecewise((0, x/Abs(b) < 1), (meijerg(((-1, 1), ()), '
        '((), (-1, 0)), x/b), True)), (x, 0, 1)), (b > 0) & (b < oo)), '
        '(Integral(Heaviside(-b + x), (x, 0, 1)), True))')
    assert integrate(exp(-x)*Heaviside(x - 1), (x, -2, 2)) == -exp(-2) + exp(-1)
    assert integrate(Heaviside(x)*x, (x, -1, 1)) == S(1)/2


# `_simpc` of `_laplace_trig_ltex`

def test_laplace_heaviside_sin_shift():
    assert str(laplace_transform(Heaviside(t - a)*sin(t), t, s, noconds=True)) == (
        '(s*cos(a - pi/2) + cos(a))*exp(-a*s)/(s**2 + 1)')
    assert str(laplace_transform(Heaviside(t - a)*sin(t), t, s)) == (
        '((s*cos(a - pi/2) + cos(a))*exp(-a*s)/(s**2 + 1), 0, True)')


def test_laplace_trig_coefficients_rewritten_in_cos():
    assert str(laplace_transform(sin(t + r), t, s, noconds=True)) == (
        '(s*cos(r - pi/2) + cos(r))/(s**2 + 1)')
    assert str(laplace_transform(sin(t + b), t, s, noconds=True)) == (
        '(s*(I*(-I*sin(b) + cos(b))/2 - I*(I*sin(b) + cos(b))/2) + cos(b))/(s**2 + 1)')
    assert str(laplace_transform(cos(2*t + r)*sin(t), t, s, noconds=True)) == (
        '(s*cos(r - pi/2)/2 + 3*cos(r)/2)/(s**2 + 9) + '
        '(-s*cos(r - pi/2)/2 - cos(r)/2)/(s**2 + 1)')
    assert laplace_transform(sin(t), t, s, noconds=True) == 1/(s**2 + 1)
    assert laplace_transform(cos(3*t), t, s, noconds=True) == s/(s**2 + 9)


# Corpus DIFF entries: SymPy's -exp(-s)*Ei(-s) is wrong.

def test_laplace_one_over_t_plus_one_value():
    v = laplace_transform(1/(t + 1), t, s, noconds=True)
    assert v == -exp(s)*Ei(-s)
    for sv in [S(1)/2, 1, 3]:
        ref = mpmath.quad(lambda u: mpmath.exp(-float(sv)*u)/(u + 1), [0, mpmath.inf])
        assert abs(complex(v.subs(s, sv).evalf()) - complex(ref)) < 1e-10


# Corpus MISSING entries kept on purpose: SymPy's unconditional values are
# wrong where the pole reaches the interval.

def test_rational_parameter_poles_stay_held():
    c = Symbol('c', real=True)
    for f, lim in [(1/(x + c), (x, 0, 1)), (1/(x - c), (x, 0, 1)),
                   (1/(1 + c*x), (x, 0, 1)), (1/x, (x, c, 1))]:
        assert isinstance(integrate(f, lim), Integral)
    assert isinstance(integrate(1/(x - a), (x, 0, 1)), Integral)
    assert integrate(1/(x + a), (x, 0, 1)) == -log(a) + log(a + 1)
