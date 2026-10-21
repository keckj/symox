"""Inverse Laplace transforms of the irrational functions of Bateman 5.3.

The expected values are the reference outputs, except for the rules whose
reference closed form fails L[f](s) == F(s); those use the correct value and
are listed in docs/algorithmic-parity.md.
"""
import pytest

import symox as sx


s, t = sx.symbols('s t')
a, b = sx.symbols('a b', positive=True)
H = sx.Heaviside(t)
sqrt, exp, erf, erfc, pi = sx.sqrt, sx.exp, sx.erf, sx.erfc, sx.pi


def ilt(f):
    return sx.inverse_laplace_transform(f, s, t)


@pytest.mark.parametrize('f, expected', [
    (1/(sqrt(s) + a), H*(1/sqrt(pi)/sqrt(t) - a*exp(a**2*t)*erfc(a*sqrt(t)))),
    (1/(sqrt(s) + a)**3, H*(2/sqrt(pi)*(a**2*t + 1)*sqrt(t) - a*t*exp(a**2*t)*(2*a**2*t + 3)*erfc(a*sqrt(t)))),
    (1/(s*(sqrt(s) + a)), H*(1 - exp(a**2*t)*erfc(a*sqrt(t)))/a),
    (1/(sqrt(s)*(sqrt(s) + a)), H*exp(a**2*t)*erfc(a*sqrt(t))),
    (1/((s - a)*sqrt(s)), H*exp(a*t)*erf(sqrt(a)*sqrt(t))/sqrt(a)),
    (sqrt(s + a)/(s + b), H*(exp(-a*t)/sqrt(t)/sqrt(pi) + sqrt(a - b)*exp(-b*t)*erf(sqrt(a - b)*sqrt(t)))),
    (1/sqrt(s + b)/(s + a), H*(1/sqrt(b - a)*exp(-a*t)*erf(sqrt(b - a)*sqrt(t)))),
    ((sqrt(s) - a)/(s*(sqrt(s) + a)), (2*exp(a**2*t)*erfc(a*sqrt(t)) - 1)*H),
    ((sqrt(s) - a)**3/(s*(sqrt(s) + a)**3),
     H*(2*(8*a**4*t**2 + 8*a**2*t + 1)*exp(a**2*t)*erfc(a*sqrt(t)) - 8/sqrt(pi)*a*sqrt(t)*(2*a**2*t + 1) - 1)),
    (1/sqrt(s**2 + a**2), H*sx.besselj(0, a*t)),
    (1/sqrt(s**2 - a**2), H*sx.besseli(0, a*t)),
    (exp(-a*s)/sqrt(s**2 + 1), sx.Heaviside(t - a)*sx.besselj(0, a - t)),
])
def test_reference_values(f, expected):
    assert ilt(f) == expected


@pytest.mark.parametrize('f, expected', [
    # the reference has 1 - ... + (1 - 2*a**2*t)*exp(a**2*t)*(erf(a*sqrt(t)) - 1)
    (1/(sqrt(s) + a)**2, (-2*a*sqrt(t)/sqrt(pi) + (2*a**2*t + 1)*exp(a**2*t)*erfc(a*sqrt(t)))*H),
    # the reference has a**3*t**(5/2) for a*t**(3/2)
    (1/(sqrt(s) + a)**4, (-2*a*t**sx.Rational(3, 2)*(2*a**2*t + 5)/(3*sqrt(pi))
                          + t*(4*a**4*t**2 + 12*a**2*t + 3)*exp(a**2*t)*erfc(a*sqrt(t))/3)*H),
    # the reference has erfc for erf in the middle term
    (sqrt(s)/((sqrt(s) + a)*(s - b)), (a**2*exp(a**2*t)*erfc(a*sqrt(t)) + a*sqrt(b)*exp(b*t)*erf(sqrt(b)*sqrt(t))
                                       - b*exp(b*t))*H/(a**2 - b)),
    # the reference ignores the leading coefficient of the cubed pole
    (sqrt(5*s)/(2*sqrt(s) + 3)**3,
     sqrt(5)*(-3*sqrt(t)*(9*t/4 + 2)/(8*sqrt(pi))
              + (81*t**2/8 + 45*t/4 + 1)*exp(9*t/4)*erfc(3*sqrt(t)/2)/8)*H),
    (1/(sqrt(5*s)*(2*sqrt(s) + 3)**3),
     sqrt(5)*(-3*t**sx.Rational(3, 2)/(8*sqrt(pi)) + t*(9*t/2 + 1)*exp(9*t/4)*erfc(3*sqrt(t)/2)/8)*H/5),
])
def test_corrected_rules(f, expected):
    assert ilt(f) == expected


@pytest.mark.parametrize('f', [
    # rule 5.3.2 needs sqrt(a*s); the reference applies it to sqrt(s + 1)
    1/(s*(s - 1)*sqrt(s + 1)),
    # rule 5.3.14 needs sqrt(a*s) in the numerator, not sqrt(s**2)
    sqrt(s**2)/(sqrt(s) + 1)**3,
    # the pole order is undecidable; the reference raises TypeError
    1/((sqrt(s) + a)*(sqrt(s) + sx.Symbol('c', positive=True))),
    # no rule: the shift has the wrong sign
    1/(sqrt(s) - 1),
])
def test_unmatched_shapes_stay_unevaluated(f):
    assert ilt(f).has(sx.InverseLaplaceTransform)


def test_transforms_invert_numerically():
    import mpmath
    cases = [1/(2*sqrt(s) + 3)**2, 1/(2*sqrt(s) + 3)**4, sqrt(5*s)/((3*sqrt(s) + 2)*(2*s - 3)),
             sqrt(5*s)/(2*sqrt(s) + 3)**3, 1/(sqrt(5*s)*(2*sqrt(s) + 3)**3), (2*sqrt(s) - 3)**2/(5*s*(sqrt(s) + sx.Rational(3, 2))**2)]
    for F in cases:
        f = ilt(F).subs(H, 1)
        for point in (4, 7):
            value = mpmath.quad(lambda u: complex(sx.N(f.subs(t, sx.Float(str(u), 20)), 20))*mpmath.exp(-point*u),
                                [0, 1, 10, 40])
            assert abs(value - complex(sx.N(F.subs(s, point), 20))) < 1e-8
