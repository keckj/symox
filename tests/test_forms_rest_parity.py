"""Output forms against hardcoded SymPy 1.14 results: the float coefficients
kept by cancel over RR (the reference dmp_cancel clears no denominators over
RR and computes the gcd over QQ after to_rational(limit=True))."""

import pytest

import symox as sx


x, y, a, s, k10 = sx.symbols('x y a s k10')


CANCEL_CASES = [
    (lambda: sx.cancel(2.0*x/(4.0*x + 2.0)), '2.0*x/(4.0*x + 2.0)'),
    (lambda: sx.cancel((1.0*x + 0.5)/(2.0*x + 1.0)), '0.500000000000000'),
    (lambda: sx.cancel(x/(0.3*x + 0.1)), '1.0*x/(0.3*x + 0.1)'),
    (lambda: sx.cancel((0.25*x**2 - 0.25)/(0.5*x - 0.5)), '0.5*x + 0.5'),
    (lambda: sx.cancel(1.0*sx.exp(k10*s)/(1.0*k10*sx.exp(a*s) - 0.6*sx.exp(a*s))),
     '1.0*exp(k10*s)/(1.0*k10*exp(a*s) - 0.6*exp(a*s))'),
    (lambda: sx.cancel((x + 1.5)/(3.0*x + 0.1)), '(1.0*x + 1.5)/(3.0*x + 0.1)'),
    (lambda: sx.cancel(0.1*x**2/(0.2*x)), '0.5*x'),
    (lambda: sx.together(1/(0.5*x) + 1/(1.5*y)), '(0.666666666666667*x + 2.0*y)/(x*y)'),
]


@pytest.mark.parametrize('case', range(len(CANCEL_CASES)))
def test_cancel_keeps_float_coefficients(case):
    compute, expected = CANCEL_CASES[case]
    assert str(compute()) == expected


def test_float_rate_integral_generic_piece():
    result = sx.integrate(sx.exp(k10*s)*sx.exp(-1.66666666666667*s), s)
    generic, condition = result.args[0]
    assert str(generic) == ('1.0*exp(k10*s)/(1.0*k10*exp(1.66666666666667*s)'
                            ' - 1.66666666666667*exp(1.66666666666667*s))')
    assert str(condition) == 'Ne(k10, 1.66666666666667)'
