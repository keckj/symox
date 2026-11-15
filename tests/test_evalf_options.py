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
    if actual._prec == 1 and abs(sm.to_sympy(actual)) > 2 * abs(expected):
        # No correct bit: as SymPy, a precision 1 Float bounding the residue.
        assert reference.evalf(15, maxn=maxn)._prec == 1
        return
    assert abs(sm.to_sympy(actual) - expected) <= abs(expected) * sp.Rational(2) ** (1 - actual._prec)


def test_maxn_without_options_keeps_the_native_budget():
    # SymPy's default maxn returns low-accuracy noise here; the native default
    # budget resolves the sum, and an explicit strict request still raises.
    assert cancelling(sm, 300).evalf(15) == sm.Float(1)
    with pytest.raises(sm.PrecisionExhausted):
        cancelling(sm, 300).evalf(15, strict=True)


def test_quad_is_ignored_without_an_integral_and_used_for_one():
    assert sm.pi.evalf(15, quad=None, verbose=False) == sm.pi.evalf(15)
    assert sm.pi.evalf(15, quad="osc") == sp.pi.evalf(15, quad="osc")
    x = sm.Symbol("x")
    X = sp.Symbol("x")
    for n in (15, 30, 50):
        assert str(sm.Integral(sm.sin(x) / x**2, (x, 1, sm.oo)).evalf(n, quad="osc")) == \
            str(sp.Integral(sp.sin(X) / X**2, (X, 1, sp.oo)).evalf(n, quad="osc"))
    with pytest.raises(ValueError, match="^quadosc requires an infinite integration interval$"):
        sm.Integral(sm.sin(x) / x, (x, 1, 2)).evalf(quad="osc")
    with pytest.raises(ValueError, match="is required for oscillatory quadrature$"):
        sm.Integral(sm.exp(-x), (x, 1, sm.oo)).evalf(quad="osc")


@pytest.mark.parametrize("build", [
    lambda m, x: m.Integral(m.sin(x) / m.log(x), (x, 2, 3)),
    lambda m, x: m.Integral(m.cos(x**3) * m.sqrt(x), (x, 0, m.pi)),
    lambda m, x: m.Integral(1 / (1 + x**4 + m.exp(x)), (x, 0, m.oo)),
    lambda m, x: m.Integral(m.exp(m.I * x) / (1 + x + m.sin(x)**2), (x, 0, 1)),
])
def test_integral_without_closed_form_uses_quadrature(build):
    x, X = sm.Symbol("x"), sp.Symbol("x")
    for n in (15, 30, 50):
        assert str(build(sm, x).evalf(n)) == str(build(sp, X).evalf(n))
    assert str(sm.N(build(sm, x), 20)) == str(sp.N(build(sp, X), 20))


def test_integral_quadrature_keeps_symbolic_and_multiple_integrals_held():
    x, y = sm.symbols("x y")
    held = sm.Integral(sm.sin(x) / sm.log(x + y), (x, 2, 3))
    assert held.evalf() == held
    nested = sm.Integral(sm.sin(x * y) / sm.log(x + y), (x, 2, 3), (y, 2, 3))
    assert nested.evalf() == nested


def test_integral_quadrature_where_the_reference_loses_accuracy():
    # SymPy's quadrature returns 0.50 (sin(1) - Ci(1) through doit here)
    # and 20 correct digits (a Bessel closed form through doit here).
    x = sm.Symbol("x")
    assert str(sm.Integral(sm.sin(1 / x), (x, 0, 1)).evalf(30)) == "0.504067061906928371989856117741"
    assert str(sm.Integral(sm.exp(-x**2) * sm.cos(x**3), (x, 0, sm.oo)).evalf(50)) == \
        "0.69405413348439055027960469448964639985266233467441"


def test_strict_succeeds_when_the_target_is_reached():
    assert sm.pi.evalf(50, strict=True) == sp.pi.evalf(50, strict=True)
    assert sm.Integer(5).evalf(15, strict=True) == sp.Integer(5).evalf(15, strict=True)
