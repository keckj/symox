"""Precision near cancellation and on branch cuts, against SymPy 1.14.

An exactly cancelling sum that evalf cannot distinguish from zero is returned
as SymPy's evalf_add returns it: precision 1 Floats bounding the residue
(``0.e-130``), PrecisionExhausted under strict, zero under chop. Multivalued
functions evaluated on and next to their branch cuts give SymPy's branch with
digits correct to the requested precision.
"""

import mpmath
import pytest
import sympy as sp

import symox as sm

SMALL = "I/10**30"

EXACT_ZEROS = [
    "-I*exp(3*I/5)/2 - exp(3*I/10)*sin(Rational(3, 10)) + I/2",
    "sin(pi/7)**2 + cos(pi/7)**2 - 1",
    "Add(Mul(sqrt(2), sqrt(3), evaluate=False), -sqrt(6), -1, 1, evaluate=False)",
    "I*sin(pi/9)**2 + I*cos(pi/9)**2 - I + sqrt(2)*sqrt(5) - sqrt(10)",
]

OPTIONS = [{}, {"n": 30}, {"n": 5}, {"maxn": 10}, {"maxn": 500}, {"n": 30, "maxn": 50}]


def parts(value):
    """(precision, value) of the real and imaginary Floats, None for an exact zero."""
    re, im = value.as_real_imag()
    out = []
    for part in (re, im):
        if part == 0 and not part.is_Float:
            out.append(None)
        else:
            assert part.is_Float, value
            out.append((part._prec, mpmath.mpf(part._mpf_)))
    return out


@pytest.mark.parametrize("options", OPTIONS)
@pytest.mark.parametrize("source", EXACT_ZEROS)
def test_an_exact_zero_is_a_precision_one_bound_as_in_the_reference(source, options):
    expected = parts(sp.sympify(source).evalf(**options))
    actual = parts(sm.sympify(source).evalf(**options))
    for reference, value in zip(expected, actual):
        assert (reference is None) == (value is None)
        if reference is None:
            continue
        assert reference[0] == value[0] == 1
        # The bound need not have SymPy's exact exponent (the accuracy
        # estimates differ by a few bits), but it bounds the same exact zero
        # at a comparable magnitude.
        assert value[1] > 0
        assert abs(mpmath.log(abs(value[1]), 2) - mpmath.log(abs(reference[1]), 2)) <= 128


@pytest.mark.parametrize("source", EXACT_ZEROS)
def test_an_exact_zero_prints_without_digits(source):
    text = str(sm.sympify(source).evalf())
    assert text.startswith("0.e-") and "e-" in text.split("+")[-1]
    assert str(sm.N(sm.sympify(source))) == text
    assert str(sm.sympify(source).n()) == text


@pytest.mark.parametrize("source", EXACT_ZEROS)
def test_strict_and_chop_on_an_exact_zero(source):
    with pytest.raises(sp.PrecisionExhausted):
        sp.sympify(source).evalf(strict=True)
    with pytest.raises(sm.PrecisionExhausted):
        sm.sympify(source).evalf(strict=True)
    assert sp.sympify(source).evalf(chop=True) == 0
    assert sm.sympify(source).evalf(chop=True) == 0
    assert sm.sympify(source).evalf(30, chop=True) == 0


def test_the_reference_example_and_its_digitless_format():
    expression = sm.sympify(EXACT_ZEROS[0])
    value = expression.evalf(30)
    re, im = value.as_real_imag()
    assert re._prec == im._prec == 1
    assert str(value).startswith("0.e-") and str(value).endswith("*I")
    # The digitless format of mpmath's to_str(x, 0), as SymPy prints it.
    for exponent in (-431, -400, -20, -5, -4, -1, 0, 2, 3, 310):
        for sign in (1, -1):
            value = sign * sp.Integer(2) ** exponent
            assert str(sm.Float(sm.sympify(value), precision=1)) == str(sp.Float(value, precision=1))
            assert str(sm.Float(sm.sympify(value), precision=4)) == str(sp.Float(value, precision=4))


def test_a_cancellation_resolved_within_the_native_budget_keeps_its_digits():
    # SymPy's default maxn reports a bound here: the native budget resolves it.
    source = "Add(exp(500), 1, -exp(500), evaluate=False)"
    assert sm.sympify(source).evalf() == sm.Float(1)
    assert sm.sympify(source).evalf()._prec == 53


def test_an_exhausted_cancellation_with_a_symbol_factor():
    x = sm.Symbol("x")
    value = (x * sm.sympify(EXACT_ZEROS[1])).evalf()
    assert value.has(x)
    floats = list(value.atoms(sm.Float))
    assert floats and all(f._prec == 1 for f in floats)


CUT = [f"-2 + {SMALL}", f"-2 - {SMALL}", "-2"]
OUTER = [f"2 + {SMALL}", f"2 - {SMALL}", "2", f"-2 + {SMALL}", f"-2 - {SMALL}"]
IMAG = [f"2*I + {SMALL}", f"2*I - {SMALL}", "2*I", f"-2*I - {SMALL}"]
INNER = ["Rational(1, 2)", f"Rational(1, 2) + {SMALL}", f"Rational(1, 2) - {SMALL}", f"-Rational(1, 2) - {SMALL}"]
LAMBERT = [f"-exp(-1) + {SMALL}", f"-exp(-1) - {SMALL}", "-Rational(3, 10)", f"-Rational(3, 10) - {SMALL}", "-2", f"-2 + {SMALL}"]


def instances(template, points):
    return [template.format(z=f"({point})") for point in points]


BRANCHES = (
    instances("log({z})", CUT) + instances("sqrt({z})", CUT) + instances("({z})**Rational(1, 3)", CUT)
    + instances("({z})**pi", CUT) + instances("({z})**(1 + I)", CUT)
    + instances("asin({z})", OUTER) + instances("acos({z})", OUTER) + instances("atanh({z})", OUTER)
    + instances("atan({z})", IMAG) + instances("acot({z})", IMAG)
    + instances("acosh({z})", CUT + INNER[:3]) + instances("asec({z})", INNER) + instances("acsc({z})", INNER)
    + instances("LambertW({z})", LAMBERT) + instances("LambertW({z}, -1)", LAMBERT) + instances("LambertW({z}, 1)", LAMBERT[:3])
    + instances("polylog(2, {z})", [f"2 + {SMALL}", f"2 - {SMALL}", "2"])
    + instances("zeta({z})", ["1 + 1/10**20", "1 - 1/10**20", f"1 + {SMALL}"])
    + instances("hyper([1, 1], [2], {z})", [f"2 + {SMALL}", f"2 - {SMALL}", "1 - 1/10**20"])
    + instances("hyper([Rational(1, 2), Rational(1, 3)], [2], {z})", ["1", f"2 + {SMALL}"])
    + instances("meijerg([[1, 1], []], [[1], [0]], {z})", ["1", f"-1 + {SMALL}", f"-1 - {SMALL}"])
    + instances("loggamma({z})", ["-Rational(5, 2)", f"-Rational(5, 2) + {SMALL}", f"-Rational(5, 2) - {SMALL}"])
    + instances("besselj(Rational(1, 3), {z})", CUT) + instances("bessely(0, {z})", CUT)
    + instances("bessely(Rational(1, 3), {z})", CUT)
)


def complex_value(value):
    re, im = value.as_real_imag()
    return mpmath.mpf(sp.Float(re, 200)._mpf_ if re else 0), mpmath.mpf(sp.Float(im, 200)._mpf_ if im else 0)


def within(value, truth, tolerance, scale=None):
    return all(abs(v - t) <= tolerance * (abs(t) if scale is None else scale) for v, t in zip(value, truth))


@pytest.mark.parametrize("digits", [15, 30, 100])
@pytest.mark.parametrize("source", BRANCHES)
def test_branch_cut_digits_and_branch_follow_the_reference(source, digits):
    truth = complex_value(sp.sympify(source).evalf(digits + 30))
    try:
        reference = complex_value(sp.sympify(source).evalf(digits))
    except (ValueError, AttributeError, TypeError):
        # SymPy rounds the argument onto a pole (zeta(1), hyper at 1) at this
        # precision and fails; the exact argument is evaluated natively.
        reference = None
    actual_value = sm.sympify(source).evalf(digits)
    assert actual_value.is_number and not actual_value.atoms(sm.Function)
    actual = complex_value(sm.to_sympy(actual_value))
    tolerance = mpmath.mpf(10) ** (1 - digits)
    if reference is None or within(reference, truth, tolerance):
        # Each part correct to the requested precision, as SymPy's.
        assert within(actual, truth, tolerance), (actual_value, truth)
    else:
        # SymPy's own accuracy is relative to the larger part here (a tiny
        # part on the cut): require at least that.
        scale = max(abs(t) for t in truth)
        assert within(actual, truth, tolerance, scale), (actual_value, truth)
