"""Special function API boundaries and genuine arbitrary precision values."""

import mpmath as mp
import pytest

import symox as sx


@pytest.mark.parametrize(
    "name,args",
    [
        ("zeta", ()),
        ("zeta", (2, 3, 4)),
        ("beta", ()),
        ("beta", (1, 2, 3)),
        ("harmonic", (1, 2, 3)),
        ("elliptic_e", ()),
        ("elliptic_pi", (1,)),
        ("LambertW", (1, 2, 3)),
        ("besselj", (1,)),
    ],
)
def test_invalid_arity_raises_typeerror(name, args):
    with pytest.raises(TypeError):
        getattr(sx, name)(*args)
    with pytest.raises(TypeError):
        sx._core._numeric(name, *args)


@pytest.mark.parametrize(
    "name,ratios,reference",
    [
        ("factorial", [(1, 3)], mp.factorial),
        ("binomial", [(5, 3), (1, 7)], mp.binomial),
        ("RisingFactorial", [(1, 3), (2, 7)], mp.rf),
        ("FallingFactorial", [(4, 3), (2, 7)], mp.ff),
        ("beta", [(1, 3), (2, 7)], mp.beta),
        ("beta", [(1, 3)], lambda x: mp.beta(x, x)),
        ("catalan", [(1, 3)], lambda x: 4**x * mp.gamma(x + mp.mpf("0.5")) / (mp.sqrt(mp.pi) * mp.gamma(x + 2))),
        ("fibonacci", [(1, 3)], mp.fibonacci),
        ("lucas", [(1, 3)], lambda x: mp.phi**x + mp.cospi(x) * mp.phi**(-x)),
    ],
)
def test_special_arbitrary_precision(name, ratios, reference):
    expr = getattr(sx, name)(*[sx.Rational(n, d) for n, d in ratios])
    value = sx.special_nstr(expr, 60)
    assert value is not None
    with mp.workdps(80):
        want = reference(*[mp.mpf(n) / d for n, d in ratios])
        assert abs(mp.mpf(value) - want) < mp.mpf("1e-58") * max(1, abs(want))


def test_arbitrary_precision_bessel_digits_are_certified():
    # The digits come from a multiprecision evaluation, not a 53-bit fallback.
    value = sx.special_nstr(sx.besselj(0, sx.Rational(1, 3)), 60)
    with mp.workdps(80):
        want = mp.besselj(0, mp.mpf(1) / 3)
        assert abs(mp.mpf(value) - want) < mp.mpf("1e-58") * abs(want)


def test_binary_function_cannot_be_used_as_unary_numeric_helper():
    with pytest.raises(TypeError):
        sx._core._numeric_many("besselj", [1.0, 2.0])
