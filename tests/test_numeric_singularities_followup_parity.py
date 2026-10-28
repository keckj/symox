"""Numerical parameter tuples at poles and exact E_nu(0) sign facts.

SymPy 1.14 is the differential reference; runtime evaluation remains native
and mpmath-backed. Near a pole the argument must not be rounded onto it.
"""

import pytest
import sympy as sp
import symox as sx


BALANCED = [
    ([(1, 1), (1, 1)], [(2, 1)]),
    ([(1, 2), (1, 2)], [(1, 1)]),
    ([(1, 1), (2, 1)], [(3, 1)]),
    ([(1, 3), (1, 3)], [(2, 3)]),
    ([(1, 3), (1, 3), (1, 3)], [(1, 2), (1, 2)]),
]


def _hyper(module, parameters, argument):
    ap, bq = parameters
    return module.hyper([module.Rational(*a) for a in ap],
                        [module.Rational(*b) for b in bq], argument)


@pytest.mark.parametrize("parameters", BALANCED)
@pytest.mark.parametrize("floating", [False, True])
@pytest.mark.parametrize("digits", [15, 50])
def test_hyper_parameter_tuples_at_logarithmic_singularity(parameters, floating, digits):
    for options in ({}, {"strict": True}, {"maxn": 100}):
        arguments = [module.Float(1, 60) if floating else module.S.One
                     for module in (sx, sp)]
        actual = _hyper(sx, parameters, arguments[0]).evalf(digits, **options)
        reference = _hyper(sp, parameters, arguments[1]).evalf(digits, **options)
        assert sx.to_sympy(actual) == reference == sp.oo


@pytest.mark.parametrize("floating", [False, True])
def test_hyper_singularity_through_substitution(floating):
    actual_x, reference_x = sx.Symbol("x"), sp.Symbol("x")
    actual = sx.hyper([1, 1], [2], actual_x)
    reference = sp.hyper([1, 1], [2], reference_x)
    point = sx.Float(1) if floating else sx.S.One
    reference_point = sp.Float(1) if floating else sp.S.One
    assert sx.to_sympy(actual.evalf(30, subs={actual_x: point})) == \
        reference.evalf(30, subs={reference_x: reference_point}) == sp.oo
    assert sx.to_sympy(actual.subs(actual_x, point).evalf(30)) == sp.oo


@pytest.mark.parametrize("digits", [15, 30, 50])
@pytest.mark.parametrize("offset", ["-Rational(1, 10**8)", "I/10**12", "-I/10**12"])
def test_hyper_near_singularity_keeps_argument_and_branch(digits, offset):
    actual_z = 1 + eval(offset, vars(sx))
    reference_z = 1 + eval(offset, vars(sp))
    actual = sx.to_sympy(sx.hyper([1, 1], [2], actual_z).evalf(digits))
    reference = sp.hyper([1, 1], [2], reference_z).evalf(digits + 20)
    assert actual.is_finite is True
    assert abs((actual - reference).evalf(digits + 10)) < \
        sp.Rational(1, 10**(digits - 2))*abs(reference)


def test_hyper_convergent_boundary_and_symbolic_parameters_are_not_infinite():
    assert sx.to_sympy(sx.hyper([1, 1], [3], sx.Float(1)).evalf(30)) == sp.Float(2, 30)
    a = sx.Symbol("a")
    value = sx.hyper([a, 1], [2], 1).evalf(30)
    assert a in value.free_symbols and value.has(sx.hyper)


@pytest.mark.parametrize("digits", [15, 30])
def test_exact_argument_rounded_onto_pole_is_refined_before_nonfinite_fallback(digits):
    # SymPy evaluates at low working precision and returns oo here. Its
    # higher-precision value verifies that this is a finite nearby point.
    actual = sx.to_sympy(sx.hyper([1, 1], [2],
                                 1 - sx.Rational(1, 10**40)).evalf(digits))
    reference = sp.hyper([1, 1], [2], 1 - sp.Rational(1, 10**40)).evalf(80)
    assert actual.is_finite is True
    assert abs((actual - reference).evalf(80)) < \
        sp.Rational(1, 10**(digits - 2))*abs(reference)


@pytest.mark.parametrize("digits", [15, 30])
def test_exact_parameter_balance_is_not_inferred_from_rounded_tuples(digits):
    parameters = [sx.Rational(1, 3)]*2
    denominator = sx.Rational(2, 3) + sx.Rational(1, 10**40)
    actual = sx.to_sympy(sx.hyper(parameters, [denominator], 1).evalf(digits))
    reference = sp.hyper([sp.Rational(1, 3)]*2,
                         [sp.Rational(2, 3) + sp.Rational(1, 10**40)], 1).evalf(80)
    assert actual.is_finite is True
    assert abs((actual - reference).evalf(80)) < \
        sp.Rational(1, 10**(digits - 2))*abs(reference)


def test_hyper_interrupt_is_not_retried_by_nonfinite_tuple_fallback(monkeypatch):
    import symox._mpmath_special as host

    calls = []

    def interrupted(*args):
        calls.append(args)
        raise KeyboardInterrupt

    monkeypatch.setitem(host._GENERIC, "hyper", interrupted)
    # A different working precision avoids values cached by the pole tests.
    with pytest.raises(KeyboardInterrupt):
        sx.hyper([1, 1], [2], sx.Float('0.5')).evalf(43)
    assert len(calls) == 1


@pytest.mark.parametrize("order", [1, 2, 3, 10])
@pytest.mark.parametrize("floating", [False, True])
def test_uppergamma_negative_integer_at_zero(order, floating):
    zero = sx.Float(0) if floating else sx.S.Zero
    reference_zero = sp.Float(0) if floating else sp.S.Zero
    actual = sx.uppergamma(-order, zero)
    reference = sp.uppergamma(-order, reference_zero)
    assert sx.to_sympy(actual) == reference == sp.zoo
    assert actual.evalf(30, strict=True) == sx.zoo
    x = sx.Symbol("x")
    assert sx.uppergamma(-order, x).subs(x, zero) == sx.zoo


@pytest.mark.parametrize("floating", [False, True])
def test_zero_order_and_undecided_orders_keep_their_reference_forms(floating):
    for module in (sx, sp):
        zero = module.Float(0) if floating else module.S.Zero
        n = module.Symbol("n", integer=True, nonnegative=True)
        assert module.uppergamma(0, zero) == module.expint(1, zero)
        assert module.uppergamma(-n, zero).func == module.uppergamma


def test_expint_positive_sign_at_zero_is_proved_without_changing_its_form():
    for order in (2, 3):
        value = sx.expint(order, sx.Float(0))
        assert value.func == sx.expint
        assert value.is_positive is True
        assert value.is_extended_positive is True
    # All-Float arguments evaluate at construction, as in the reference.
    assert sx.to_sympy(sx.expint(sx.Float(2, 30), sx.Float(0))) == \
        sp.expint(sp.Float(2, 30), sp.Float(0))
    # For E_1(0) the integral diverges; do not use the finite positive rule.
    assert sx.expint(1, 0).is_positive is not True
