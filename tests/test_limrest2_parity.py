"""limrest2 (wave 14): Gruntz limits through digamma/loggamma/gamma at
logarithmically growing arguments, gamma differences and imaginary
coefficients (c0*oo).  Expected values are those of SymPy 1.14."""

from symox import *
from symox.series.gruntz import gruntz

del test

x = Symbol("x")


def test_polygamma_loggamma_nested():
    assert gruntz(digamma(digamma(x)), x, oo) == oo
    assert gruntz(digamma(digamma(digamma(x))), x, oo) == oo
    assert gruntz(loggamma(loggamma(x)), x, oo) == oo
    assert limit(loggamma(loggamma(x)), x, oo) == oo
    assert gruntz(exp(exp(digamma(digamma(x)))) / x, x, oo) == exp(-S.Half)
    assert gruntz(exp(exp(digamma(log(x)))) / x, x, oo) == exp(-S.Half)
    assert limit(exp(exp(digamma(log(x)))) / x, x, oo) == exp(-S.Half)


def test_gamma_differences():
    assert gruntz(exp(x) * (gamma(x + exp(-x)) - gamma(x)), x, oo) is oo
    assert gruntz(gamma(x + 1) / sqrt(2 * pi)
                  - exp(-x) * (x**(x + S.Half) + x**(x - S.Half) / 12), x, oo) is oo
    assert gruntz(exp(gamma(x - exp(-x)) * exp(1 / x)) - exp(gamma(x)), x, oo) is oo


def test_imaginary_coefficient():
    assert limit(x * Shi(x) * I, x, oo) == oo * I
    assert limit(I * Shi(x), x, oo) == oo * I
    assert limit(I * x * exp(x), x, oo) == oo * I
    assert limit(x * Shi(x), x, oo) == oo


def test_bessel_at_imaginary_infinity():
    assert limit(besselj(1, x), x, I * oo) == oo * I
    assert limit(besselj(1, x), x, -I * oo) == -oo * I
    assert limit(besselk(1, x), x, -oo) == -oo * I


def test_unsupported_stay_unevaluated():
    # SymPy returns the meaningless Ci(1/_w)/_w**2 here; the limit does not
    # exist (x**2*Ci(x) ~ x*sin(x)).
    assert isinstance(limit(x**2 * Ci(x), x, oo), Limit)
