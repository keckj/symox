"""limrest3 (wave 21): limits left unevaluated after limrest2.

Expected values are those of SymPy 1.14 (sympy/series/tests/test_limits.py,
test_gruntz.py and functions/combinatorial/tests/test_comb_numbers.py),
checked numerically below where the limit is a constant.
"""

import mpmath

from symox import *
from symox.series.gruntz import gruntz

del test

x = Symbol("x")
n = Symbol("n")


def test_ei_erf_dominant_cancellation():
    # Taylor expansion of the auxiliary _eis/_erfs at a slower argument
    e = (Ei(x - exp(-exp(x))) - Ei(x)) * exp(-x) * exp(exp(x)) * x
    assert gruntz(e, x, oo) == -1
    assert limit(e, x, oo) == -1
    e = (erf(x - exp(-exp(x))) - erf(x)) * exp(exp(x)) * exp(x**2)
    assert gruntz(e, x, oo) == -2 / sqrt(pi)
    assert limit(e, x, oo) == -2 / sqrt(pi)


def test_ei_erf_small_instances():
    # the same expansion with a polynomially small shift
    assert limit((Ei(x + 1/x) - Ei(x)) * x**2 * exp(-x), x, oo) == 1
    assert limit((erf(x + 1/x) - erf(x)) * x * exp(x**2), x, oo) == -exp(-2)/sqrt(pi) + 1/sqrt(pi)
    assert limit((erfi(x + 1/x) - erfi(x)) * x * exp(-x**2), x, oo) == -1/sqrt(pi) + exp(2)/sqrt(pi)


def test_uppergamma_lowergamma_order_at_infinity():
    assert limit(uppergamma(x, 1) / gamma(x), x, oo) == 1
    assert limit(lowergamma(x, 1), x, oo) == 0
    assert limit(x * lowergamma(x, 1), x, oo) == exp(-1)
    assert gruntz(x * lowergamma(x, 1), x, oo) == exp(-1)


def test_lowergamma_order_numeric():
    # x*lowergamma(x, 1) -> 1/e and x*lowergamma(x, 2)/2**x -> exp(-2)
    with mpmath.workdps(30):
        v = 10**6 * mpmath.gammainc(10**6, 0, 1)
        assert abs(v - mpmath.exp(-1)) < 1e-5
    assert limit(x * lowergamma(x, 2) / 2**x, x, oo) == exp(-2)


def test_subfactorial_at_infinity():
    assert limit(subfactorial(x), x, oo) is oo


def test_stirling_remainder():
    e = ((x + S(1) / 2) * log(x) - x + log(2 * pi) / 2 - log(factorial(x)) + S(1) / (12 * x)) * x**3
    assert limit(e, x, oo) == S(1) / 360
    e = (x + S(1) / 2) * log(x) - x + log(2 * pi) / 2 - loggamma(x + 1)
    assert limit(e, x, oo) == 0
    assert limit(e * x, x, oo) == -S(1) / 12


def test_harmonic_value_at_infinity():
    m = Symbol("m", positive=True)
    y = Symbol("y", positive=True)
    z = Symbol("z", negative=True)
    assert limit(harmonic(n, m + 1), n, oo) == zeta(m + 1)
    assert harmonic(x, 1).limit(x, 0) == 0
    assert harmonic(x, y).limit(x, 0) == 0
    assert limit(harmonic(x, y), x, oo) == harmonic(oo, y)
    assert limit(harmonic(x, y + 1), x, oo) == zeta(y + 1)
    assert limit(harmonic(x, y - 1), x, oo) == harmonic(oo, y - 1)
    assert limit(harmonic(x, z), x, oo) == Limit(harmonic(x, z), x, oo, dir='-')
    assert limit(harmonic(x, z + 1), x, oo) == oo
    assert limit(harmonic(x, z + 2), x, oo) == harmonic(oo, z + 2)
    assert limit(harmonic(x, z - 1), x, oo) == Limit(harmonic(x, z - 1), x, oo, dir='-')
    # an integer order still goes through the polygamma form
    assert limit(harmonic(x, 3), x, oo) == zeta(3)


def test_harmonic_symbolic_order_not_real_stays_held():
    w = Symbol("w")
    assert isinstance(limit(harmonic(x + I, 3 + w), x, oo), Limit)


def test_atan_parametric_limit_and_integral():
    a = Symbol("a")
    ar = Symbol("a", real=True)
    y = Symbol("y", positive=True)
    assert limit(a * (log(x) - log(a**2 * x**2 + 1) / 2) - atan(a * x) / x, x, oo) == -a * log(a**2) / 2
    assert limit(ar * (log(x) - log(ar**2 * x**2 + 1) / 2) - atan(ar * x) / x, x, oo) == -ar * log(ar**2) / 2
    assert limit(atan(a * x) / x, x, oo) == 0
    expected = -a * log(a**2) / 2 + a * log(a**2 + 1) / 2 + atan(a)
    assert integrate(atan(a * x) / x**2, (x, 1, oo)) == expected
    expected = -ar * log(ar**2) / 2 + ar * log(ar**2 + 1) / 2 + atan(ar)
    r = integrate(atan(ar * x) / x**2, (x, 1, oo))
    assert r == expected
    for v in (-3, -S(1) / 2, S(1) / 3, 2):
        num = mpmath.quad(lambda t: mpmath.atan(v * t) / t**2, [1, mpmath.inf])
        assert abs(r.subs(ar, v).evalf(20) - num) < 1e-12
    r = integrate(atan(y * x) / x**2, (x, 1, oo))
    assert (r - (-y * log(y) + y * log(y**2 + 1) / 2 + atan(y))).expand(force=True) == 0


def test_acoth_log_singularity_two_sided():
    # acoth(u) with u -> +-1: the reference log form
    # (log(1 + 1/u) - log(1 - 1/u))/2, real part dominant
    for d in ('+', '-', '+-'):
        assert limit(acoth(sin(x) / x), x, 0, d) is oo
        assert limit(acoth(1 - x**2), x, 0, d) is oo
        assert limit(acoth(1 - x), x, 0, d) is oo
        assert limit(acoth(x - 1), x, 0, d) is -oo
        assert limit(acoth(cos(x)), x, 0, d) is oo
        assert limit(acoth(sin(x) / x) - log(x), x, 0, d) is oo
        assert limit(acoth(1 + x) + log(x) / 2, x, 0, d) == log(2) / 2
        assert limit(x * acoth(1 - x), x, 0, d) == 0
    e = acoth(1 - x) + log(x) / 2
    assert limit(e, x, 0, '+') == log(2) / 2 - I * pi / 2
    assert limit(e, x, 0, '-') == log(2) / 2 + I * pi / 2
    with mpmath.workdps(30):
        t = mpmath.mpf(10)**-12
        v = mpmath.acoth(1 - t) + mpmath.log(t) / 2
        assert abs(v - (mpmath.log(2) / 2 - 1j * mpmath.pi / 2)) < 1e-10
    try:
        limit(e, x, 0, '+-')
        assert False
    except ValueError:
        pass
