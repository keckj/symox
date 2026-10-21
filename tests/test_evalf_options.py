"""evalf strict, maxn, chop and quad options follow SymPy's error behavior."""

import pytest
import sympy as sp

import symox as sm


def cancelling(module, magnitude):
    return module.Add(module.exp(magnitude), 1, -module.exp(magnitude), evaluate=False)


EXACT_ZEROS = [
    lambda m: m.sin(1) ** 2 + m.cos(1) ** 2 - 1,
    lambda m: m.Add(m.Mul(m.sqrt(2), m.sqrt(3), evaluate=False), -m.sqrt(6), evaluate=False),
]


@pytest.mark.parametrize("build", EXACT_ZEROS + [lambda m: cancelling(m, 300), lambda m: cancelling(m, 1000)])
def test_strict_raises_where_the_reference_cannot_distinguish_from_zero(build):
    with pytest.raises(sp.PrecisionExhausted):
        build(sp).evalf(15, strict=True)
    with pytest.raises(sm.PrecisionExhausted) as error:
        build(sm).evalf(15, strict=True)
    assert str(error.value).startswith("Failed to distinguish the expression: \n\n")
    assert str(error.value).endswith("from zero. Try simplifying the input, using chop=True, "
                                     "or providing a higher maxn for evalf")


@pytest.mark.parametrize("magnitude", [300, 1000])
def test_a_larger_maxn_resolves_the_cancellation(magnitude):
    expected = cancelling(sp, magnitude).evalf(15, strict=True, maxn=1000)
    actual = cancelling(sm, magnitude).evalf(15, strict=True, maxn=1000)
    assert actual.is_Float and actual == sm.Float(1) and expected == sp.Float(1)
    assert actual._prec == expected._prec == 53


@pytest.mark.parametrize("build", EXACT_ZEROS)
def test_chop_accepts_an_exact_zero_under_strict(build):
    assert build(sp).evalf(15, strict=True, chop=True) == 0
    assert build(sm).evalf(15, strict=True, chop=True) == 0


def test_strict_reports_the_failing_subexpression():
    x = sm.Symbol("x")
    with pytest.raises(sm.PrecisionExhausted, match="exp"):
        sm.sin(x).evalf(15, strict=True, subs={x: cancelling(sm, 1000)})
    with pytest.raises(sm.PrecisionExhausted):
        (sm.I * (sm.sin(1) ** 2 + sm.cos(1) ** 2 - 1) + 2).evalf(15, strict=True)


@pytest.mark.parametrize("digits", [30, 40])
@pytest.mark.parametrize("maxn", [10, 20, 50])
def test_maxn_returns_the_accuracy_it_reached(digits, maxn):
    approximation = sp.Rational(int(sp.pi.evalf(digits + 5) * 10 ** digits), 10 ** digits)
    reference = sp.Add(sp.pi, -approximation, evaluate=False)
    actual = sm.Add(sm.pi, -sm.sympify(approximation), evaluate=False).evalf(15, maxn=maxn)
    expected = reference.evalf(60)
    assert actual.is_Float and 1 <= actual._prec <= 53
    assert abs(sm.to_sympy(actual) - expected) <= abs(expected) * sp.Rational(2) ** (1 - actual._prec)


def test_maxn_without_options_keeps_the_native_budget():
    # SymPy's default maxn returns low-accuracy noise here; the native default
    # budget resolves the sum, and an explicit strict request still raises.
    assert cancelling(sm, 300).evalf(15) == sm.Float(1)
    with pytest.raises(sm.PrecisionExhausted):
        cancelling(sm, 300).evalf(15, strict=True)


def test_quad_is_ignored_without_an_integral_and_refused_for_one():
    assert sm.pi.evalf(15, quad=None, verbose=False) == sm.pi.evalf(15)
    assert sm.pi.evalf(15, quad="osc") == sp.pi.evalf(15, quad="osc")
    x = sm.Symbol("x")
    with pytest.raises(NotImplementedError):
        sm.Integral(sm.sin(x) / x, (x, 1, sm.oo)).evalf(15, quad="osc")


def test_strict_succeeds_when_the_target_is_reached():
    assert sm.pi.evalf(50, strict=True) == sp.pi.evalf(50, strict=True)
    assert sm.Integer(5).evalf(15, strict=True) == sp.Integer(5).evalf(15, strict=True)
