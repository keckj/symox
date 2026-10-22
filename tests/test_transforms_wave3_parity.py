"""Transform and integration gaps closed against the reference SymPy 1.14.

The expected strings are the outputs of SymPy 1.14 on the same inputs,
except where a test documents a remaining difference.
"""
import pytest

import symox as sx
from symox import (Float, Piecewise, Rational, Symbol, besselj, exp, integrate, sin, symbols)
from symox import cosine_transform, inverse_fourier_transform, sine_transform
from symox.integrals.heurisch import heurisch, heurisch_wrapper

t, w, k, x, s = symbols('t w k x s')
k10 = Symbol('k10')


# --- Piecewise results: the reference simplify of a Piecewise is
# piecewise_simplify, whose univariate conditions become real intervals ---

PIECEWISE_CASES = [
    (lambda: sine_transform(besselj(0, t), t, w, noconds=False),
     '(Piecewise((0, (w > -1) & (w < 1)), (sqrt(2)/(sqrt(pi)*w*sqrt(1 - 1/w**2)), (w > 1) | (w < -1)),'
     ' (sqrt(2)*meijerg(((1/2,), (1/2,)), ((1/2,), (0,)), w**2), True)), Eq(Abs(arg(w)), 0) & Ne(w**2, 1))'),
    (lambda: cosine_transform(besselj(0, t), t, w, noconds=False),
     '(Piecewise((sqrt(2)/(sqrt(pi)*sqrt(1 - w**2)), (w > -1) & (w < 1)), (0, (w > 1) | (w < -1)),'
     ' (sqrt(2)*meijerg(((1/2,), (1/2,)), ((0,), (1/2,)), w**2), True)), Eq(Abs(arg(w)), 0) & Ne(w**2, 1))'),
    (lambda: sine_transform(sin(t)/t, t, w, noconds=False),
     '(Piecewise((sqrt(2)*atanh(1/w)/sqrt(pi), (w > 1) | (w < -1)), (sqrt(2)*atanh(w)/sqrt(pi), (w > -1) & (w < 1)),'
     ' (sqrt(2)*sqrt(pi)*meijerg(((1/2,), (1,)), ((1/2,), (0,)), w**(-2))/2, True)),'
     ' Eq(Abs(arg(w)), 0) & Ne(w**(-2), 1))'),
    (lambda: inverse_fourier_transform(sin(k)/k, k, x, noconds=False),
     '(Piecewise((0, (x > 1/(2*pi)) | (x < -1/(2*pi))), (pi, True)),'
     ' ((Abs(2*arg(x) + pi) < pi) | (Eq(Abs(2*arg(x) + pi), pi) & Ne(1/(pi**2*x**2), 4)))'
     ' & ((Abs(2*arg(x) - pi) < pi) | (Ne(1/(pi**2*x**2), 4) & Eq(Abs(2*arg(x) - pi), pi))))'),
]


@pytest.mark.parametrize('compute, expected', PIECEWISE_CASES)
def test_piecewise_conditions_as_real_intervals(compute, expected):
    assert str(compute()) == expected


def test_cosine_transform_sinc_intervals():
    # the outer conditions are the reference ones; the reference keeps a
    # nested Piecewise in the Meijer G piece, which is folded here
    v, c = cosine_transform(sin(t)/t, t, w, noconds=False)
    assert str(v.args[0]) == '(0, (w > 1) | (w < -1))'
    assert str(v.args[1]) == '(sqrt(2)*sqrt(pi)/2, (w > -1) & (w < 1))'
    assert str(c) == 'Eq(Abs(arg(w)), 0) & Ne(w**(-2), 1)'


def test_unsimplified_transform_keeps_conditions():
    v, _ = sine_transform(besselj(0, t), t, w, noconds=False, simplify=False)
    assert v.has(sx.Abs)


# --- test_issue_15498: heurisch_wrapper with a Float rate. The linear
# system of heurisch is solved over RR(k10) like the reference, so the
# generic antiderivative has the denominator k10 - c and the degenerate
# rate is a separate case. ---

C = Float('1.66666666666667')


def _value(e, subs):
    return complex(sx.N(e.subs(subs)))


@pytest.mark.parametrize('call', [
    lambda f: integrate(f, s),
    lambda f: integrate(f, s, heurisch=True),
    lambda f: heurisch_wrapper(f, s),
])
def test_float_rate_degenerate_case(call):
    r = call(exp(k10*s)*exp(-C*s))
    assert isinstance(r, Piecewise)
    (generic, cond), (special, other) = r.args
    assert str(cond) == 'Ne(k10, 1.66666666666667)'
    assert other == True  # noqa: E712
    assert not special.has(k10)
    assert not cond.has(s)
    for kv, sv in [(Rational(1, 3), Rational(7, 10)), (Rational(5, 2), Rational(-3, 10))]:
        ref = exp((kv - C)*sv)/(kv - C)
        assert abs(_value(generic, {k10: kv, s: sv}) - complex(sx.N(ref))) < 1e-12*abs(complex(sx.N(ref)))


def test_float_rate_heurisch_is_first_order():
    h = heurisch(exp(k10*s)*exp(-C*s), s)
    _, d = sx.fraction(sx.together(h))
    assert sx.degree(d.subs(exp(C*s), 1), k10) == 1


@pytest.mark.parametrize('call, expected', [
    ('integrate(exp(k10*s)*exp(-Rational(5, 3)*s), s)',
     'Piecewise((3*exp(k10*s)/(3*k10*exp(5*s/3) - 5*exp(5*s/3)), Ne(k10, 5/3)), (s, True))'),
])
def test_exact_rate_unchanged(call, expected):
    assert str(eval(call)) == expected
