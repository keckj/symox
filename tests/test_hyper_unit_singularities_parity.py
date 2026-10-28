"""Exact convergence excess and signed pFq(1) singularities.

The native certificate uses the exact gamma-ratio sign of the asymptotic
series coefficient. Terminating and forbidden parameters are kept on their
existing paths. SymPy's erroneous large finite values at a pole are not
reference values; positive-excess Gauss values are checked at high precision.
"""

import pytest
import sympy as sp
import symox as sx


def _hyper(module, ap, bq, z=1):
    return module.hyper([module.sympify(a) for a in ap],
                        [module.sympify(b) for b in bq], z)


DIVERGENT = [
    (["1/3"], [], 1),
    (["1/4"], [], 1),
    (["1/3", "1/3"], ["1/2"], 1),
    (["1/3", "1/3", "1/3"], ["1/4", "1/4"], 1),
    (["1/7", "2/7", "4/7"], ["1/5", "4/5"], 1),
    (["-1/3", "-1/3"], ["-2/3"], -1),
    (["-2/3", "1/3"], ["-1/3"], 1),
    (["-3/2", "1"], ["-3/4"], -1),
    (["-4/3", "2/3"], ["-2/3"], -1),
]


@pytest.mark.parametrize("ap,bq,sign", DIVERGENT)
@pytest.mark.parametrize("floating", [False, True])
@pytest.mark.parametrize("digits", [15, 50])
def test_exact_signed_pfq_unit_divergence(ap, bq, sign, floating, digits):
    z = sx.Float(1, 70) if floating else sx.S.One
    expression = _hyper(sx, ap, bq, z)
    for options in ({}, {"strict": True}, {"maxn": 100}):
        assert expression.evalf(digits, **options) == sign*sx.oo


@pytest.mark.parametrize("parameter", ["1", "2", "1/2", "3/2"])
def test_native_one_f_zero_poles_replace_reference_zero_division(parameter):
    # mpmath's special integer/half-integer power route raises for 0**(-a).
    # The positive defining series proves divergence to +oo instead.
    assert _hyper(sx, [parameter], []).evalf(30) == sx.oo


def test_certified_unit_poles_do_not_call_the_numeric_host(monkeypatch):
    import symox._mpmath_special as host

    def unavailable(*args):
        raise AssertionError("a certified unit pole must not call mpmath")

    monkeypatch.setitem(host._GENERIC, "hyper", unavailable)
    for ap, bq, sign in DIVERGENT:
        assert _hyper(sx, ap, bq).evalf(43) == sign*sx.oo


def test_exact_balance_does_not_accept_the_hosts_large_finite_numbers():
    # Independent rounding of thirds/quarters made the old backend return
    # large finite values on some members of this exactly balanced family.
    for i in range(12):
        a, b = sx.Rational(1, 3) + sx.Rational(i, 1001), sx.Rational(1, 4)
        assert sx.hyper([a, b], [a + b], 1).evalf(30) == sx.oo


@pytest.mark.parametrize("source,expected", [
    ("h + 1", "oo"),
    ("-2*h", "-oo"),
    ("1/h", "0"),
    ("h**2", "oo"),
    ("h**-2", "0"),
    ("sqrt(h)", "oo"),
    ("log(h)", "oo"),
    ("exp(-h)", "0"),
    ("Mul(0, h, evaluate=False)", "nan"),
    ("Pow(h, 0, evaluate=False)", "1"),
    ("-hn", "oo"),
    ("sqrt(hn)", "I*oo"),
])
def test_certified_unit_poles_propagate_through_scalar_compositions(source, expected, monkeypatch):
    import symox._mpmath_special as host

    def unavailable(*args):
        raise AssertionError("certified compositions must not call mpmath")

    h = _hyper(sx, ["1/3"], [])
    hn = _hyper(sx, ["-1/3", "-1/3"], ["-2/3"])
    expression = eval(source, dict(vars(sx), h=h, hn=hn))
    monkeypatch.setitem(host._GENERIC, "hyper", unavailable)
    actual = sx.to_sympy(expression.evalf(40))
    reference = sp.sympify(eval(expected, vars(sp)))
    if reference.is_finite is True:
        assert (actual - reference).is_zero is True
    else:
        assert actual == reference


@pytest.mark.parametrize("ap,bq,sign", [(["0.5", "0.5"], ["1.0"], 1),
                                        (["0.5", "0.75"], ["1.0"], 1),
                                        (["-1.5", "1.0"], ["-0.75"], -1)])
def test_float_parameters_use_their_exact_binary_values(ap, bq, sign):
    expression = sx.hyper([sx.Float(a, 70) for a in ap],
                          [sx.Float(b, 70) for b in bq], sx.Float(1, 70))
    assert expression.evalf(50) == sign*sx.oo


@pytest.mark.parametrize("ap,bq,sign", DIVERGENT[:4])
def test_unit_poles_through_subs_and_evalf_subs(ap, bq, sign):
    z = sx.Symbol("z")
    expression = _hyper(sx, ap, bq, z)
    assert expression.subs(z, sx.Float(1)).evalf(30) == sign*sx.oo
    assert expression.evalf(30, subs={z: 1}, strict=True) == sign*sx.oo


@pytest.mark.parametrize("ap,bq", [(["1/3", "1/4"], ["7/6"]),
                                   (["1/3", "2/3"], ["2"]),
                                   (["2", "1/3"], ["4"])])
@pytest.mark.parametrize("digits", [30, 60])
def test_gauss_positive_excess_matches_high_precision_reference(ap, bq, digits):
    actual = sx.to_sympy(_hyper(sx, ap, bq).evalf(digits, strict=True, maxn=150))
    reference = _hyper(sp, ap, bq).evalf(digits + 30)
    assert actual.is_finite is True
    assert abs((actual - reference).evalf(digits + 30)) < \
        sp.Rational(1, 10**(digits - 2))*abs(reference)


@pytest.mark.parametrize("digits", [15, 30, 60])
def test_gauss_tiny_positive_excess_is_not_rounded_onto_a_pole(digits):
    actual = sx.to_sympy(sx.hyper([sx.Rational(1, 3)]*2,
                         [sx.Rational(2, 3) + sx.Rational(1, 10**40)], 1).evalf(digits))
    reference = sp.hyper([sp.Rational(1, 3)]*2,
                         [sp.Rational(2, 3) + sp.Rational(1, 10**40)], 1).evalf(100)
    assert actual.is_finite is True
    assert abs((actual - reference).evalf(100)) < \
        sp.Rational(1, 10**(digits - 2))*abs(reference)


def test_gauss_float_parameters_preserve_more_than_machine_precision():
    parameters = [sx.Float("0.2", 75), sx.Float("0.3", 80)]
    denominator = sx.Float("1.2", 85)
    actual = sx.to_sympy(sx.hyper(parameters, [denominator], sx.Float(1, 90)).evalf(60))
    reference = sp.hyper([sx.to_sympy(a) for a in parameters],
                         [sx.to_sympy(denominator)], 1).evalf(100)
    assert abs((actual - reference).evalf(100)) < sp.Rational(1, 10**58)*abs(reference)


@pytest.mark.parametrize("ap,bq,expected", [(["-2", "1/3"], ["4/3"], "9/14"),
                                           (["0", "2/3"], ["1/3"], "1"),
                                           (["-3", "-2"], ["-4"], "0"),
                                           (["-2", "5/3", "7/4"], ["2/3", "9/4"], "-10/39")])
def test_terminating_parameters_are_not_classified_as_poles(ap, bq, expected):
    for floating in (False, True):
        expression = (sx.hyper([sx.Float(sx.sympify(a), 70) for a in ap],
                               [sx.Float(sx.sympify(b), 70) for b in bq], 1)
                      if floating else _hyper(sx, ap, bq))
        actual = sx.to_sympy(expression.evalf(40))
        reference = sp.Rational(expected)
        assert abs((actual - reference).evalf(60)) < sp.Rational(1, 10**38)


def test_forbidden_and_symbolic_parameters_keep_existing_behavior():
    for module in (sx, sp):
        for denominator in (-1, module.Float(-1)):
            with pytest.raises(ZeroDivisionError, match="pole in hypergeometric series"):
                module.hyper([1], [denominator], 1).evalf()
        a = module.Symbol("a")
        assert a in module.hyper([a], [], 1).evalf().free_symbols


@pytest.mark.parametrize("point", ["-1", "1 - Rational(1, 10**40)",
                                   "1 + I/10**12", "1 - I/10**12"])
def test_other_boundary_points_and_cut_sides_are_not_unit_poles(point):
    actual_z, reference_z = eval(point, vars(sx)), eval(point, vars(sp))
    actual = sx.to_sympy(_hyper(sx, ["1/3"], [], actual_z).evalf(30))
    reference = _hyper(sp, ["1/3"], [], reference_z).evalf(90)
    assert actual.is_finite is True
    assert abs((actual - reference).evalf(90)) < sp.Rational(1, 10**28)*abs(reference)


def test_oscillating_series_boundary_keeps_its_analytic_continuation():
    # The defining sum at -1 does not converge for a=3, but the continued
    # hypergeometric function is (1-z)**(-3), finite at that point.
    actual = sx.to_sympy(sx.hyper([3], [], -1).evalf(40))
    assert (actual - sp.Rational(1, 8)).is_zero is True


def test_uncertified_hyper_interrupt_stops_without_retry(monkeypatch):
    import symox._mpmath_special as host

    calls = []

    def interrupted(*args):
        calls.append(args)
        raise KeyboardInterrupt

    monkeypatch.setitem(host._GENERIC, "hyper", interrupted)
    with pytest.raises(KeyboardInterrupt):
        sx.hyper([1, 1], [2], sx.Rational(1, 2)).evalf(47)
    assert len(calls) == 1
