"""Fourier, sine, cosine and Hankel transforms and their inverses.

The transforms are computed like the reference: through the definite
integral, whose first Piecewise condition is the returned convergence
condition. The expected strings are the outputs of SymPy 1.14.
"""
import pytest

import symox as sx
from symox import (Function, Heaviside, I, Rational, S, Symbol, besselj, cos, exp, gamma, log, meijerg, oo,
                   pi, sin, sqrt, symbols)
from symox import (cosine_transform, fourier_transform, hankel_transform, inverse_cosine_transform,
                   inverse_fourier_transform, inverse_hankel_transform, inverse_sine_transform, laplace_transform,
                   sine_transform, CosineTransform, InverseFourierTransform, InverseHankelTransform, SineTransform)
from symox.integrals.transforms import IntegralTransformError


t, w, x, s = symbols('t w x s')
k, r, nu, m = symbols('k r nu m')
a = Symbol('a')
ap = Symbol('a', positive=True)
bp = Symbol('b', positive=True)
kr = Symbol('k', real=True)
posk = Symbol('posk', positive=True)
negk = Symbol('k', negative=True)
f = Function('f')


REFERENCE = [
    ('sit_undef',
     lambda: sine_transform(f(t), t, w, noconds=False),
     'SineTransform(f(t), t, w)'),
    ('isit_undef',
     lambda: inverse_sine_transform(f(w), w, t, noconds=False),
     'InverseSineTransform(f(w), w, t)'),
    ('sit_isqrt',
     lambda: sine_transform(1/sqrt(t), t, w, noconds=False),
     '(1/sqrt(w), Eq(Abs(arg(w)), 0))'),
    ('isit_isqrt',
     lambda: inverse_sine_transform(1/sqrt(w), w, t, noconds=False),
     '(1/sqrt(t), Eq(Abs(arg(t)), 0))'),
    ('sit_isqrt3',
     lambda: sine_transform((1/sqrt(t))**3, t, w, noconds=False),
     '(2*sqrt(w), Eq(Abs(arg(w)), 0))'),
    ('sit_tpow',
     lambda: sine_transform(t**(-a), t, w, noconds=False),
     '(2**(1/2 - a)*w**(a - 1)*gamma(1 - a/2)/gamma(a/2 + 1/2), (re(a) > 0) & (re(a) < 2) & Eq(Abs(arg(w)), 0))'),
    ('isit_tpow',
     lambda: inverse_sine_transform(2**(-a + S(1)/2)*w**(a - 1)*gamma(-a/2 + 1)/gamma(a/2 + S.Half), w, t, noconds=False),
     '(t**(-a), (re(a) > -1) & (re(a) < 1) & Eq(Abs(arg(t)), 0))'),
    ('sit_exp',
     lambda: sine_transform(exp(-a*t), t, w, noconds=False),
     '(sqrt(2)*w/(sqrt(pi)*(a**2 + w**2)), Eq(Abs(arg(w)), 0) & (Abs(arg(a)) < pi/2))'),
    ('isit_exp',
     lambda: inverse_sine_transform(sqrt(2)*w/(sqrt(pi)*(a**2 + w**2)), w, t, noconds=False),
     '(exp(-a*t), Eq(Abs(arg(t)), 0) & (Abs(arg(a)) < pi/2))'),
    ('sit_gauss',
     lambda: sine_transform(t*exp(-a*t**2), t, w, noconds=False),
     '(sqrt(2)*w*exp(-w**2/(4*a))/(4*a**(3/2)), (Abs(arg(a)) < pi/2) | (Eq(Abs(arg(w)), 0) & (Abs(arg(a)) < pi/2)))'),
    ('isit_gauss',
     lambda: inverse_sine_transform(sqrt(2)*w*exp(-w**2/(4*a))/(4*a**Rational(3, 2)), w, t, noconds=False),
     '(t*exp(-a*t**2), (Abs(arg(a)) < pi/2) | (Eq(Abs(arg(t)), 0) & (Abs(arg(a)) < pi/2)))'),
    ('sit_exp_pos',
     lambda: sine_transform(exp(-ap*t), t, w, noconds=False),
     '(sqrt(2)*w/(sqrt(pi)*(a**2 + w**2)), Eq(Abs(arg(w)), 0))'),
    ('sit_rat',
     lambda: sine_transform(t/(1 + t**2), t, w, noconds=False),
     '(sqrt(2)*sqrt(pi)*exp(-w)/2, Eq(Abs(arg(w)), 0))'),
    ('sit_rat_a',
     lambda: sine_transform(t/(a**2 + t**2), t, w, noconds=False),
     '(sqrt(2)*sqrt(pi)*exp(-a*w)/2, Eq(Abs(arg(w)), 0) & (Abs(arg(a)) < pi/2))'),
    ('sit_exp_t',
     lambda: sine_transform(exp(t), t, w, noconds=False),
     'SineTransform(exp(t), t, w)'),
    ('sit_1',
     lambda: sine_transform(S(1), t, w, noconds=False),
     'SineTransform(1, t, w)'),
    ('sit_heav',
     lambda: sine_transform(Heaviside(1 - t), t, w, noconds=False),
     '(sqrt(2)*(1 - cos(w))/(sqrt(pi)*w), Eq(Abs(arg(w)), 0))'),
    ('sit_noconds',
     lambda: sine_transform(exp(-a*t), t, w),
     'sqrt(2)*w/(sqrt(pi)*(a**2 + w**2))'),
    ('sit_nosimp',
     lambda: sine_transform(exp(-a*t), t, w, simplify=False, noconds=False),
     '(sqrt(2)/(sqrt(pi)*w*(a**2/w**2 + 1)), Eq(Abs(arg(w)), 0) & (Abs(arg(a)) < pi/2))'),
    ('cot_undef',
     lambda: cosine_transform(f(t), t, w, noconds=False),
     'CosineTransform(f(t), t, w)'),
    ('icot_undef',
     lambda: inverse_cosine_transform(f(w), w, t, noconds=False),
     'InverseCosineTransform(f(w), w, t)'),
    ('cot_isqrt',
     lambda: cosine_transform(1/sqrt(t), t, w, noconds=False),
     '(1/sqrt(w), Eq(Abs(arg(w)), 0))'),
    ('icot_isqrt',
     lambda: inverse_cosine_transform(1/sqrt(w), w, t, noconds=False),
     '(1/sqrt(t), Eq(Abs(arg(t)), 0))'),
    ('cot_lor',
     lambda: cosine_transform(1/(a**2 + t**2), t, w, noconds=False),
     '(sqrt(2)*sqrt(pi)*exp(-a*w)/(2*a), Eq(Abs(arg(w)), 0) & (Abs(arg(a)) < pi/2))'),
    ('cot_tpow',
     lambda: cosine_transform(t**(-a), t, w, noconds=False),
     '(2**(1/2 - a)*w**(a - 1)*gamma(1/2 - a/2)/gamma(a/2), (re(a) > 0) & (re(a) < 1) & Eq(Abs(arg(w)), 0))'),
    ('icot_tpow',
     lambda: inverse_cosine_transform(2**(-a + S(1)/2)*w**(a - 1)*gamma(-a/2 + S.Half)/gamma(a/2), w, t, noconds=False),
     '(t**(-a), (re(a) > 0) & (re(a) < 1) & Eq(Abs(arg(t)), 0))'),
    ('cot_exp',
     lambda: cosine_transform(exp(-a*t), t, w, noconds=False),
     '(sqrt(2)*a/(sqrt(pi)*(a**2 + w**2)), Eq(Abs(arg(w)), 0) & (Abs(arg(a)) < pi/2))'),
    ('icot_exp',
     lambda: inverse_cosine_transform(sqrt(2)*a/(sqrt(pi)*(a**2 + w**2)), w, t, noconds=False),
     '(exp(-a*t), Eq(Abs(arg(t)), 0) & (Abs(arg(a)) < pi/2))'),
    ('cot_expcos',
     lambda: cosine_transform(exp(-a*sqrt(t))*cos(a*sqrt(t)), t, w, noconds=False),
     '(a*exp(-a**2/(2*w))/(2*w**(3/2)), Eq(Abs(arg(w)), 0) & (Abs(arg(a)) < pi/4))'),
    ('icot_mg',
     lambda: inverse_cosine_transform(sqrt(2)*meijerg(((S.Half, 0), ()), ((S.Half, 0, 0), (S.Half,)), a**2*w**2/4)/(2*pi), w, t, noconds=False),
     '(1/(a + t), Eq(Abs(arg(t)), 0) & (Abs(arg(a**2)) < 2*pi))'),
    ('cot_isq',
     lambda: cosine_transform(1/sqrt(a**2 + t**2), t, w, noconds=False),
     '(sqrt(2)*meijerg(((1/2,), ()), ((0, 0), (1/2,)), a**2*w**2/4)/(2*sqrt(pi)), Eq(Abs(arg(w)), 0) & (Abs(arg(a)) < pi/2))'),
    ('icot_mg2',
     lambda: inverse_cosine_transform(sqrt(2)*meijerg(((S.Half,), ()), ((0, 0), (S.Half,)), a**2*w**2/4)/(2*sqrt(pi)), w, t, noconds=False),
     '(1/(t*sqrt(a**2/t**2 + 1)), (Eq(Abs(arg(t)), 0) & (Abs(arg(a**2)) < pi)) | (Eq(Abs(arg(t)), 0) & Eq(Abs(arg(a**2)), pi) & Ne(a**2/t**2, -1)))'),
    ('cot_gauss',
     lambda: cosine_transform(exp(-t**2), t, w, noconds=False),
     '(sqrt(2)*exp(-w**2/4)/2, True)'),
    ('cot_gauss_a',
     lambda: cosine_transform(exp(-a*t**2), t, w, noconds=False),
     '(sqrt(2)*exp(-w**2/(4*a))/(2*sqrt(a)), (Abs(arg(a)) < pi/2) | (Eq(Abs(arg(w)), 0) & (Abs(arg(a)) < pi/2)))'),
    ('cot_rat1',
     lambda: cosine_transform(1/(1 + t**2), t, w, noconds=False),
     '(sqrt(2)*sqrt(pi)*exp(-w)/2, Eq(Abs(arg(w)), 0))'),
    ('cot_exp_pos',
     lambda: cosine_transform(exp(-ap*t), t, w, noconds=False),
     '(sqrt(2)*a/(sqrt(pi)*(a**2 + w**2)), Eq(Abs(arg(w)), 0))'),
    ('cot_1',
     lambda: cosine_transform(S(1), t, w, noconds=False),
     'CosineTransform(1, t, w)'),
    ('cot_heav',
     lambda: cosine_transform(Heaviside(1 - t), t, w, noconds=False),
     '(sqrt(2)*sin(w)/(sqrt(pi)*w), Eq(Abs(arg(w)), 0))'),
    ('cot_noconds',
     lambda: cosine_transform(exp(-a*t), t, w),
     'sqrt(2)*a/(sqrt(pi)*(a**2 + w**2))'),
    ('ft_undef',
     lambda: fourier_transform(f(x), x, kr, noconds=False),
     'FourierTransform(f(x), x, k)'),
    ('ift_undef',
     lambda: inverse_fourier_transform(f(kr), kr, x, noconds=False),
     'InverseFourierTransform(f(k), k, x)'),
    ('ft_box',
     lambda: fourier_transform(Heaviside(1 - abs(2*ap*x)), x, kr, noconds=False),
     '(sin(pi*k/a)/(pi*k), k > 0)'),
    ('ft_exph',
     lambda: fourier_transform(exp(-ap*x)*Heaviside(x), x, kr, noconds=False),
     '(1/(a + 2*I*pi*k), Abs(arg(k) + pi/2) <= pi/2)'),
    ('ift_pos',
     lambda: inverse_fourier_transform(1/(ap + 2*pi*I*x), x, posk, noconds=False),
     '(exp(-a*posk), True)'),
    ('ift_negpos',
     lambda: inverse_fourier_transform(1/(ap + 2*pi*I*x), x, -posk, noconds=False),
     '(0, True)'),
    ('ift_neg',
     lambda: inverse_fourier_transform(1/(ap + 2*pi*I*x), x, negk, noconds=False),
     '(0, True)'),
    ('ift_k',
     lambda: inverse_fourier_transform(1/(ap + 2*pi*I*x), x, kr, noconds=False),
     'InverseFourierTransform(1/(a + 2*I*pi*x), x, k)'),
    ('ft_gauss',
     lambda: fourier_transform(exp(-ap*x**2), x, kr, noconds=False),
     '(sqrt(pi)*exp(-pi**2*k**2/a)/sqrt(a), True)'),
    ('ift_gauss',
     lambda: inverse_fourier_transform(sqrt(pi/ap)*exp(-(pi*kr)**2/ap), kr, x, noconds=False),
     '(exp(-a*x**2), True)'),
    ('ift_gauss_k',
     lambda: inverse_fourier_transform(exp(-pi*k**2), k, x, noconds=False),
     '(exp(-pi*x**2), True)'),
    ('ift_gauss_k1',
     lambda: inverse_fourier_transform(exp(-k**2), k, x, noconds=False),
     '(sqrt(pi)*exp(-pi**2*x**2), True)'),
    ('ft_absexp',
     lambda: fourier_transform(exp(-ap*abs(x)), x, kr, noconds=False),
     '(2*a/(a**2 + 4*pi**2*k**2), k > 0)'),
    ('ift_lor',
     lambda: inverse_fourier_transform(2*ap/(ap**2 + 4*pi**2*kr**2), kr, x, noconds=False),
     '(exp(-a*x), x > 0)'),
    ('ift_lor1',
     lambda: inverse_fourier_transform(1/(1 + k**2), k, x, noconds=False),
     '(pi*exp(-2*pi*x), x > 0)'),
    ('ift_exp',
     lambda: inverse_fourier_transform(exp(kr), kr, x, noconds=False),
     'InverseFourierTransform(exp(k), k, x)'),
    ('ft_exp_real',
     lambda: fourier_transform(exp(x), x, kr, noconds=False),
     'FourierTransform(exp(x), x, k)'),
    ('ift_noconds',
     lambda: inverse_fourier_transform(exp(-pi*k**2), k, x),
     'exp(-pi*x**2)'),
    ('ht_inv',
     lambda: hankel_transform(1/r, r, k, 0, noconds=False),
     '(1/k, Eq(Abs(arg(k)), 0))'),
    ('iht_inv',
     lambda: inverse_hankel_transform(1/k, k, r, 0, noconds=False),
     '(1/r, Eq(Abs(arg(r)), 0))'),
    ('ht_pow0',
     lambda: hankel_transform(1/r**m, r, k, 0, noconds=False),
     '(2**(1 - m)*k**(m - 2)*gamma(1 - m/2)/gamma(m/2), (re(m) > 1/2) & (re(m) < 2) & Eq(Abs(arg(k)), 0))'),
    ('iht_pow0',
     lambda: inverse_hankel_transform(2**(-m + 1)*k**(m - 2)*gamma(-m/2 + 1)/gamma(m/2), k, r, 0, noconds=False),
     '(r**(-m), (re(m) > 0) & (re(m) < 3/2) & Eq(Abs(arg(r)), 0))'),
    ('iht_pownu',
     lambda: inverse_hankel_transform(2**(-m + 1)*k**(m - 2)*gamma(-m/2 + nu/2 + 1)/gamma(m/2 + nu/2), k, r, nu, noconds=False),
     '(r**(-m), (re(m) < 3/2) & Eq(Abs(arg(r)), 0) & (re(m)/2 + re(nu)/2 > 0))'),
    ('ht_exp',
     lambda: hankel_transform(r**nu*exp(-a*r), r, k, nu, noconds=False),
     '(2**(nu + 1)*a*k**(-nu - 3)*(a**2/k**2 + 1)**(-nu - 3/2)*gamma(nu + 3/2)/sqrt(pi), ((re(nu) > -3/2) & (re(nu) > -1) & Eq(Abs(arg(k)), 0) & (Abs(arg(a)) < pi/2)) | ((re(nu) > -3/2) & (re(nu) > -1) & (re(nu) < -1/2) & Eq(Abs(arg(k)), 0) & Eq(Abs(arg(a)), pi/2) & ((re(nu) < -3/2) | Ne(a**2/k**2, -1))))'),
    ('iht_exp',
     lambda: inverse_hankel_transform(2**(nu + 1)*a*k**(-nu - 3)*(a**2/k**2 + 1)**(-nu - Rational(3, 2))*gamma(nu + Rational(3, 2))/sqrt(pi), k, r, nu, noconds=False),
     '(r**nu*exp(-a*r), (re(nu) > -5/2) & (re(nu) > -1) & Eq(Abs(arg(r)), 0) & (Abs(arg(a)) < pi/2))'),
    ('ht_undef',
     lambda: hankel_transform(f(r), r, k, nu, noconds=False),
     'ERR IntegralTransformError: Hankel Transform could not be computed: could not compute integral.'),
    ('iht_undef',
     lambda: inverse_hankel_transform(f(k), k, r, nu, noconds=False),
     'ERR IntegralTransformError: Inverse Hankel Transform could not be computed: could not compute integral.'),
    ('ht_gauss',
     lambda: hankel_transform(exp(-r**2), r, k, 0, noconds=False),
     '(exp(-k**2/4)/2, True)'),
    ('iht_gauss',
     lambda: inverse_hankel_transform(exp(-k**2/4)/2, k, r, 0, noconds=False),
     '(exp(-r**2), True)'),
    ('ht_exp0',
     lambda: hankel_transform(exp(-a*r), r, k, 0, noconds=False),
     '(a/(k**3*(a**2/k**2 + 1)**(3/2)), Eq(Abs(arg(k)), 0) & (Abs(arg(a)) < pi/2))'),
    ('ht_lor',
     lambda: hankel_transform(1/(r**2 + 1), r, k, 0, noconds=False),
     '(meijerg(((0,), ()), ((0, 0), (0,)), k**2/4)/2, Eq(Abs(arg(k)), 0))'),
    ('ht_noconds',
     lambda: hankel_transform(1/r, r, k, 0),
     '1/k'),
    ('iht_noconds',
     lambda: inverse_hankel_transform(1/k, k, r, 0),
     '1/r'),
    ('lt_sinsqrt',
     lambda: laplace_transform(sin(sqrt(t)), t, s),
     '(sqrt(pi)*exp(-1/(4*s))/(2*s**(3/2)), 0, True)'),
    ('ft_doc',
     lambda: fourier_transform(exp(-x**2), x, k, noconds=False),
     '(sqrt(pi)*exp(-pi**2*k**2), True)'),
    ('ift_doc',
     lambda: inverse_fourier_transform(sqrt(pi)*exp(-(pi*k)**2), k, x, noconds=False),
     '(exp(-x**2), True)'),
    ('ft_nosimp',
     lambda: fourier_transform(exp(-x**2), x, k, simplify=False, noconds=False),
     '(sqrt(pi)*(2 - erfc(I*pi*k))*exp(-pi**2*k**2)/2 + sqrt(pi)*exp(-pi**2*k**2)*erfc(I*pi*k)/2, True)'),
    ('ft_needeval',
     lambda: fourier_transform(f(x), x, k, needeval=True),
     'ERR IntegralTransformError: Fourier Transform could not be computed: needeval.'),
    ('sit_needeval',
     lambda: sine_transform(f(t), t, w, needeval=True),
     'ERR IntegralTransformError: Sine Transform could not be computed: needeval.'),
    ('cot_needeval',
     lambda: cosine_transform(exp(t), t, w, needeval=True),
     'ERR IntegralTransformError: Cosine Transform could not be computed: needeval.'),
    ('sit_coeff_undef',
     lambda: sine_transform(3*a*f(t), t, w, noconds=False),
     '3*a*SineTransform(f(t), t, w)'),
    ('isit_doc',
     lambda: inverse_sine_transform(2**((1 - 2*a)/2)*k**(a - 1)*gamma(-a/2 + 1)/gamma((a + 1)/2), k, x, noconds=False),
     '(x**(-a), (re(a) > -1) & (re(a) < 1) & Eq(Abs(arg(x)), 0))'),
    ('ht_gauss_a',
     lambda: hankel_transform(exp(-ap*r**2), r, k, 0, noconds=False),
     '(exp(-k**2/(4*a))/(2*a), True)'),
    ('ht_gauss_nu',
     lambda: hankel_transform(r**nu*exp(-r**2), r, k, nu, noconds=False),
     '(k**nu*exp(-k**2/4)/(2*2**nu), (re(nu) > -1) | ((re(nu) > -1) & Eq(Abs(arg(k)), 0)))'),
    ('iht_noconds2',
     lambda: inverse_hankel_transform(exp(-k**2/4)/2, k, r, 0),
     'exp(-r**2)'),
    ('ht_nosimp',
     lambda: hankel_transform(1/r, r, k, 0, simplify=False, noconds=False),
     '(1/k, Eq(Abs(arg(k)), 0))'),
    ('ht_sum',
     lambda: hankel_transform(1/r + exp(-r), r, k, 0, noconds=False),
     '((k**2 + 1)**(-3/2) + 1/k, Eq(Abs(arg(k)), 0))'),
    ('cot_gauss_noconds',
     lambda: cosine_transform(exp(-t**2), t, w),
     'sqrt(2)*exp(-w**2/4)/2'),
    ('ift_exp_abs',
     lambda: inverse_fourier_transform(exp(-abs(k)), k, x, noconds=False),
     '(2/(4*pi**2*x**2 + 1), x > 0)'),
    ('ft_xgauss',
     lambda: fourier_transform(x*exp(-x**2), x, k, noconds=False),
     '(-I*pi**(3/2)*k*exp(-pi**2*k**2), True)'),
    ('icot_lor',
     lambda: inverse_cosine_transform(1/(1 + w**2), w, t, noconds=False),
     '(sqrt(2)*sqrt(pi)*exp(-t)/2, Eq(Abs(arg(t)), 0))'),
    ('isit_frac',
     lambda: inverse_sine_transform(w/(1 + w**2), w, t, noconds=False),
     '(sqrt(2)*sqrt(pi)*exp(-t)/2, Eq(Abs(arg(t)), 0))'),
    ('held_sit_doit',
     lambda: SineTransform(exp(-a*t), t, w).doit(noconds=False),
     '(sqrt(2)*w/(sqrt(pi)*(a**2 + w**2)), Eq(Abs(arg(w)), 0) & (Abs(arg(a)) < pi/2))'),
    ('held_iht_doit',
     lambda: InverseHankelTransform(1/k, k, r, 0).doit(noconds=False),
     '(1/r, Eq(Abs(arg(r)), 0))'),
    ('held_ift_doit',
     lambda: InverseFourierTransform(exp(-pi*k**2), k, x).doit(),
     'exp(-pi*x**2)'),
]


@pytest.mark.parametrize('name, compute, expected', REFERENCE, ids=[c[0] for c in REFERENCE])
def test_transform_reference_output(name, compute, expected):
    if expected.startswith('ERR '):
        kind, msg = expected[4:].split(': ', 1)
        with pytest.raises(Exception) as err:
            compute()
        assert (type(err.value).__name__, str(err.value)) == (kind, msg)
        return
    assert str(compute()) == expected


@pytest.mark.parametrize('compute', [
    lambda: fourier_transform(f(x), x, k, needeval=True),
    lambda: sine_transform(f(t), t, w, needeval=True),
    lambda: cosine_transform(exp(t), t, w, needeval=True),
])
def test_needeval_raises(compute):
    with pytest.raises(IntegralTransformError):
        compute()


def _close(u, v, subs):
    du = complex(sx.N(u.subs(subs)))
    dv = complex(sx.N(v.subs(subs)))
    return abs(du - dv) <= 1e-10*max(1, abs(dv))


# Values equal to the reference with a different printed form; the
# conditions are the reference ones.
@pytest.mark.parametrize('compute, value, cond, subs', [
    (lambda: sine_transform(log(t)/t, t, w, noconds=False),
     sqrt(2)*sqrt(pi)*(-log(w**2) - 2*sx.EulerGamma)/4, 'Eq(Abs(arg(w)), 0)', {w: Rational(13, 10)}),
    (lambda: sine_transform(exp(-t) + 1/sqrt(t), t, w, noconds=False),
     (sqrt(2)*w**Rational(3, 2) + sqrt(pi)*(w**2 + 1))/(sqrt(pi)*sqrt(w)*(w**2 + 1)), 'Eq(Abs(arg(w)), 0)',
     {w: Rational(7, 10)}),
    (lambda: cosine_transform(1/(a + t), t, w, noconds=False),
     sqrt(2)*((pi - 2*sx.Si(a*w))*sin(a*w)/2 - cos(a*w)*sx.Ci(a*w))/sqrt(pi),
     'Eq(Abs(arg(w)), 0) & (Abs(arg(a)) < pi)', {a: 2, w: Rational(3, 2)}),
    (lambda: fourier_transform(x*exp(-ap*x)*Heaviside(x), x, kr, noconds=False),
     (ap + 2*I*pi*kr)**(-2), 'Abs(arg(k) + pi/2) <= pi/2', {ap: 2, kr: Rational(3, 10)}),
    (lambda: hankel_transform(1/r**m, r, k, nu, noconds=False),
     2*k**(m - 2)*gamma(-m/2 + nu/2 + 1)/(2**m*gamma(m/2 + nu/2)),
     '(re(m) > 1/2) & Eq(Abs(arg(k)), 0) & (re(m)/2 - re(nu)/2 < 1)', {m: Rational(3, 4), nu: 2, k: 3}),
    (lambda: fourier_transform(Heaviside(x)*Heaviside(1 - x), x, k, noconds=False),
     -I/(2*pi*k) + I*exp(-2*I*pi*k)/(2*pi*k), 'Abs(arg(k) - pi/2) <= pi/2', {k: Rational(3, 10)}),
])
def test_transform_equivalent_forms(compute, value, cond, subs):
    v, c = compute()
    assert str(c) == cond
    assert _close(v, value, subs)


def test_laplace_cos_sqrt_is_evaluated():
    F, plane, cond = laplace_transform(cos(sqrt(t)), t, s)
    assert not F.has(sx.LaplaceTransform)
    assert (plane, cond) == (0, True)
    ref = -4*I*(I*sqrt(s)*exp(1/(4*s))/4 - sqrt(pi)*sx.erf(I/(2*sqrt(s)))/8)*exp(-1/(4*s))/s**Rational(3, 2)
    assert _close(F, ref, {s: 2})
    assert _close(laplace_transform(cos(sqrt(t)), t, s, noconds=True), ref, {s: Rational(1, 3)})


def test_fourier_integral_not_computed_uses_table_value():
    # the native definite integral of exp(-a*x)*sin(b*x)*Heaviside(x) is not
    # computed; the closed form table gives the reference value with the
    # condition True instead of the reference half plane condition
    v, c = fourier_transform(exp(-ap*x)*sin(bp*x)*Heaviside(x), x, kr, noconds=False)
    assert _close(v, bp/(bp**2 + (ap + 2*I*pi*kr)**2), {ap: 2, bp: 3, kr: Rational(1, 5)})
    assert c == True  # noqa: E712


def test_divergent_hankel_raises():
    # the reference returns a meaningless Piecewise for this divergent integral
    with pytest.raises(IntegralTransformError):
        hankel_transform(S(1), r, k, 0, noconds=False)
