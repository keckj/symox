"""Transforms computed through the definite integral when no table applies.

The expected values are the reference outputs, except where noted: the
reference drops the binomial coefficients and the argument scale of
derivatives of DiracDelta, which is documented in docs/algorithmic-parity.md.
"""
import pytest

import symox as sx


t, s, x, k = sx.symbols('t s x k')
a = sx.Symbol('a')


def same_value(actual, expected):
    if actual == expected:
        return True
    difference = sx.simplify(sx.gammasimp(actual - expected))
    if difference == 0:
        return True
    for point in (sx.Rational(3, 7), sx.Rational(11, 5)):
        value = sx.N(difference.subs({s: point, k: point}), 30)
        if abs(complex(value)) > 1e-20:
            return False
    return True


@pytest.mark.parametrize('f, expected', [
    (sx.sin(sx.sqrt(t)), (sx.sqrt(sx.pi)*sx.exp(-1/(4*s))/(2*s**sx.Rational(3, 2)), 0, sx.true)),
    (sx.besselj(0, sx.sqrt(t)), (sx.exp(-1/(4*s))/s, 0, sx.true)),
    (sx.erfc(sx.sqrt(t)), (1/s - 1/(s*sx.sqrt(s + 1)), 0, sx.true)),
    (sx.exp(-a*t**2), (sx.sqrt(sx.pi)*sx.exp(s**2/(4*a))*sx.erfc(s/(2*sx.sqrt(a)))/(2*sx.sqrt(a)), 0,
                       sx.Abs(sx.arg(a)) <= sx.pi/2)),
])
def test_laplace_through_the_integral(f, expected):
    assert sx.laplace_transform(f, t, s) == expected


def test_laplace_of_cos_sqrt_matches_the_reference_value():
    value, plane, cond = sx.laplace_transform(sx.cos(sx.sqrt(t)), t, s)
    reference = (-4*sx.I*(sx.I*sx.sqrt(s)*sx.exp(1/(4*s))/4 - sx.sqrt(sx.pi)*sx.erf(sx.I/(2*sx.sqrt(s)))/8)
                 * sx.exp(-1/(4*s))/s**sx.Rational(3, 2))
    assert (plane, cond) == (0, sx.true)
    assert same_value(value, reference)


def test_laplace_without_closed_form_stays_unevaluated():
    value, plane, cond = sx.laplace_transform(sx.atan(t), t, s)
    assert value == sx.LaplaceTransform(sx.atan(t), t, s)
    assert (plane, cond) == (-sx.oo, sx.true)


def test_fourier_transform_through_the_integral():
    g = sx.meijerg
    expected = sx.sqrt(2)*(
        g(((sx.Rational(3, 4),), ()), ((0, sx.Rational(1, 4), sx.S.Half, sx.Rational(3, 4)), ()),
          sx.pi**4*k**4*sx.exp_polar(-2*sx.I*sx.pi)/16)
        + g(((sx.Rational(3, 4),), ()), ((0, sx.Rational(1, 4), sx.S.Half, sx.Rational(3, 4)), ()),
            sx.pi**4*k**4*sx.exp_polar(2*sx.I*sx.pi)/16))/(8*sx.pi**sx.Rational(3, 2))
    assert sx.fourier_transform(sx.exp(-x**4), x, k, noconds=False) == (expected, sx.true)


def test_fourier_transform_of_non_integrable_function_stays_unevaluated():
    assert sx.fourier_transform(x, x, k) == sx.FourierTransform(x, x, k)


def test_hankel_transform_through_the_integral():
    value, cond = sx.hankel_transform(sx.exp(-a*x)/x, x, k, 0, noconds=False)
    assert value == 1/(k*sx.sqrt(a**2/k**2 + 1))
    assert cond == sx.Or(sx.And(sx.Eq(sx.Abs(sx.arg(k)), 0), sx.Abs(sx.arg(a)) < sx.pi/2),
                         sx.And(sx.Eq(sx.Abs(sx.arg(k)), 0), sx.Eq(sx.Abs(sx.arg(a)), sx.pi/2),
                                sx.Ne(a**2/k**2, -1)))
    value, cond = sx.hankel_transform(x**2*sx.exp(-x**2), x, k, 0, noconds=False)
    assert same_value(value, (4 - k**2)*sx.exp(-k**2/4)/8)
    assert cond == sx.true


@pytest.mark.parametrize('f, expected', [
    (sx.exp(-x)*sx.log(x), (sx.gamma(s)*sx.polygamma(0, s), (0, sx.oo), sx.true)),
    (sx.cos(x)*sx.exp(-x), (2**s*sx.gamma(s/4)*sx.gamma(s/2 + sx.S.Half)/(4*sx.gamma(sx.S.Half - s/4)),
                            (0, sx.oo), sx.true)),
    (sx.exp(-x)/(1 + x), (sx.E*sx.gamma(s)*sx.uppergamma(1 - s, 1), (0, sx.oo), sx.true)),
    (sx.log(x)/(1 + x), (-sx.lerchphi(sx.exp_polar(sx.I*sx.pi), 2, s) + sx.lerchphi(sx.exp_polar(sx.I*sx.pi), 2, 1 - s),
                         (0, 1), sx.true)),
])
def test_mellin_transform_through_the_integral(f, expected):
    assert sx.mellin_transform(f, x, s) == expected


def test_mellin_transform_equivalent_form():
    value, strip, cond = sx.mellin_transform(sx.sinh(x)*sx.exp(-2*x), x, s)
    assert (strip, cond) == ((0, sx.oo), sx.true)
    assert same_value(value, (3**s - 1)*sx.gamma(s)/(2*3**s))


def test_mellin_strip_moves_by_the_real_part_of_the_power():
    assert sx.mellin_transform(x**a*sx.exp(-x), x, s) == (sx.gamma(a + s), (-sx.re(a), sx.oo), sx.true)


def test_derivatives_of_dirac_delta_against_smooth_weights():
    D = sx.DiracDelta
    t0 = sx.Symbol('t0', negative=True)
    assert sx.integrate(D(t, 1)*sx.exp(-s*t), (t, 0, sx.oo)) == s/2
    assert sx.integrate(3*D(t, 1)*sx.exp(-s*t), (t, t0, sx.oo)) == 3*s
    assert sx.integrate(D(t, 2)*sx.exp(-s*t), (t, -sx.oo, sx.oo)) == s**2
    assert sx.integrate(D(x, 1), x) == D(x)
    assert sx.integrate(D(2*x - 1, 1), x) == D(2*x - 1)/2
    assert sx.integrate(x*D(x, 1), x) == -sx.Heaviside(x)
    assert sx.integrate(x**2*D(x, 1), x) == 0
    assert sx.integrate(sx.exp(x)*D(x - 1, 1), (x, -1, 2)) == -sx.E


def test_derivatives_of_dirac_delta_keep_the_distributional_coefficients():
    D = sx.DiracDelta
    # x*delta''(x) = -2*delta'(x); the reference returns -DiracDelta(x)
    assert sx.integrate(x*D(x, 2), x) == -2*D(x)
    # delta'(2*x - 1) = delta'(x - 1/2)/4; the reference drops the scale
    assert sx.integrate(sx.exp(x)*D(2*x - 1, 1), (x, -1, 2)) == -sx.exp(sx.S.Half)/4


def test_derivatives_of_dirac_delta_at_an_endpoint_stay_unevaluated():
    D = sx.DiracDelta
    for f in (D(x, 1), x*D(x, 2), D(x, 1)*sx.Heaviside(x)):
        assert isinstance(sx.integrate(f, (x, 0, sx.oo)), sx.Integral)
