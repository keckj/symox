"""Precision must live in expression numbers, rather than display strings."""

import pickle
import subprocess
import sys
from decimal import Decimal
from fractions import Fraction

import pytest
import sympy as sp

import symox as sm


def assert_close(actual, expected, digits):
    converted = sm.to_sympy(actual)
    error = abs(sp.N(converted - expected, digits + 10))
    scale = max(1, abs(expected))
    assert error < scale * sp.Rational(1, 10) ** (digits - 3), (actual, expected, error)


@pytest.mark.parametrize("digits", [20, 50, 100])
@pytest.mark.parametrize("source", [
    "pi", "E", "sqrt(2)", "sin(1)", "exp(1/3)", "log(2)",
    "EulerGamma", "gamma(1/3)", "erf(1)", "sqrt(-2)",
])
def test_evalf_differential(source, digits):
    expression = sm.sympify(source)
    expected = sp.N(sp.sympify(source), digits)
    actual = sm.N(expression, digits)
    assert_close(actual, expected, digits)
    if actual.is_Float:
        assert actual._prec == sp.Float(1, digits)._prec
    assert pickle.loads(pickle.dumps(actual)) == actual


def test_decimal_constructor_and_python_float_keep_their_distinct_values():
    decimal = sm.Float("0.1", 80)
    machine = sm.Float(0.1, 80)
    assert decimal != machine
    assert decimal._prec == sp.Float("0.1", 80)._prec
    assert_close(decimal, sp.Float("0.1", 80), 80)
    assert_close(machine, sp.Float(0.1, 80), 80)
    assert isinstance(decimal, sm.Float)
    assert sm.from_sympy(sm.to_sympy(decimal)) == decimal
    assert sm.sympify(sm.srepr(decimal)) == decimal
    assert pickle.loads(pickle.dumps(decimal)) == decimal
    assert len(str(decimal)) >= 80
    assert sm.Float(decimal) == decimal
    assert sm.Float(Decimal("0.1"), 80) == decimal
    assert sm.Float("inf") == sm.oo
    assert sm.Float(float("-inf")) == -sm.oo
    assert sm.Float("nan") == sm.nan


def test_arithmetic_preserves_precision_and_small_increments():
    value = sm.Float("1.234567890123456789012345678901234567890123456789", 80)
    increment = sm.Float("1e-60", 80)
    changed = value + increment
    assert changed != value
    expected = sp.Float("1.234567890123456789012345678901234567890123456789", 80)
    expected_increment = sp.Float("1e-60", 80)
    assert_close(changed, expected + expected_increment, 80)
    for actual, target in [
        (value * sm.Rational(2, 3), expected * sp.Rational(2, 3)),
        (1 / value, 1 / expected),
        (value**7, expected**7),
        (sm.sqrt(value), sp.sqrt(expected)),
        (sm.sin(value), sp.sin(expected)),
        (sm.exp(value), sp.exp(expected)),
    ]:
        assert actual._prec == value._prec
        assert_close(actual, target, 80)
    assert value > sm.Float("1.234567890123456789012345678901234567890123456788", 80)


def test_partial_symbolic_evaluation_and_substitution():
    x = sm.Symbol("x")
    actual = (sm.pi * x + sm.Rational(1, 3)).evalf(70)
    assert actual.free_symbols == {x}
    assert_close(actual.subs({x: 2}), sp.N(2 * sp.pi + sp.Rational(1, 3), 70), 70)
    direct = sm.N(sm.sin(x) + sm.pi, 70, subs={x: sm.Rational(1, 3)})
    assert_close(direct, sp.N(sp.sin(sp.Rational(1, 3)) + sp.pi, 70), 70)
    assert (x**2).evalf(60).exp == 2
    assert (x**sm.Rational(1, 2)).evalf(60).exp._prec == sp.Float(1, 60)._prec


def test_matrix_evaluation_uses_the_same_precision_engine():
    matrix = sm.Matrix([[sm.pi, sm.sqrt(2)], [sm.Rational(1, 3), 0]])
    evaluated = sm.N(matrix, 60)
    assert evaluated.shape == matrix.shape
    assert_close(evaluated[0, 0], sp.N(sp.pi, 60), 60)
    assert_close(evaluated[1, 0], sp.N(sp.Rational(1, 3), 60), 60)
    assert evaluated[0, 1]._prec == sp.Float(1, 60)._prec
    assert evaluated[1, 1] == 0
    x = sm.Symbol("x")
    substituted = sm.Matrix([sm.sin(x)]).evalf(60, subs={x: sm.Rational(1, 3)})
    assert_close(substituted[0], sp.N(sp.sin(sp.Rational(1, 3)), 60), 60)


def test_bit_precision_and_large_magnitude_are_retained():
    value = sm.Float("1e400", precision=200)
    assert value._prec == 200
    assert sm.from_sympy(sm.to_sympy(value)) == value
    assert value / value == sm.Float("1", precision=200)
    assert int(value) == int(sp.Float("1e400", precision=200))
    assert sm.Float("1e-400", precision=200) > 0


def test_cancellation_does_not_become_spurious_exact_zero():
    expression = sm.sqrt(sm.Integer(10)**100 + 1) - sm.Integer(10)**50
    result = expression.evalf(60)
    assert result > 0
    expected = (sp.sqrt(sp.Integer(10)**100 + 1) - sp.Integer(10)**50).evalf(60)
    assert abs(sp.N(sm.to_sympy(result) / expected - 1, 70)) < sp.Rational(1, 10)**57


@pytest.mark.parametrize("value", [5e-324, -5e-324, sys.float_info.min, sys.float_info.max])
def test_float_conversion_and_pickle_preserve_ieee_boundaries(value):
    for expression in [sm.Float(value), sm.Float(value, 80)]:
        assert float(expression).hex() == value.hex()
        assert float(pickle.loads(pickle.dumps(expression))).hex() == value.hex()
        assert sm.sympify(sm.srepr(expression)) == expression


def test_exact_rational_float_conversion_has_no_double_rounding():
    fine = 2**200
    cases = [
        (1414213562373095048801688724209, 10**30),
        (1414213562373095048801688724210, 10**30),
        (2**53 + 1, 2**53),
        ((2**53 + 1) * fine + 1, 2**53 * fine),
        ((2**53 + 1) * fine - 1, 2**53 * fine),
        (2**53 + 3, 2**53),
        (1, 2**1075),
        (fine + 1, 2**1075 * fine),
        (fine - 1, 2**1075 * fine),
        (3, 2**1075),
        (2**53 - 1, 2**1075),
        (10**400 + 1, 10**400),
        (1, 10**400),
    ]
    for numerator, denominator in cases:
        for sign in [1, -1]:
            expected = float(Fraction(sign * numerator, denominator))
            assert float(sm.Rational(sign * numerator, denominator)).hex() == expected.hex()


def test_unsupported_precision_options_are_explicit():
    for n in [0, -1, 10001]:
        assert sm.N(sm.Rational(1,3), n).to_sympy() == sp.N(sp.Rational(1,3), n)
    with pytest.raises(ValueError):
        sm.Float("1.2", 20, precision=80)
    assert sm.pi.evalf(40, chop=True) == sm.N(sm.pi, 40)
    with pytest.raises(TypeError):
        sm.N(sm.pi, 40, unknown_precision_option=True)
    assert sm.pi.evalf(40, chop=False) == sm.N(sm.pi, 40)


def test_unrepresentable_exponential_and_power_results_remain_symbolic():
    huge = sm.Float("1e100", 60)
    for function in [sm.exp, sm.sinh, sm.cosh, sm.tanh]:
        for argument in [huge, -huge]:
            expression = function(argument)
            assert expression.is_Function
            assert expression.evalf(60) == expression
    for exponent in [huge, -huge]:
        expression = sm.Float("10", 60)**exponent
        assert expression.is_Pow
        assert expression.evalf(60) == expression
    # Binary exponents are bounded by 2**100, not by a 64-bit counter.
    power = sm.Float("10", 60)**(sm.Integer(2)**63 - 1)
    assert power.is_Float
    assert power.to_sympy() == sp.Float("10", 60)**(sp.Integer(2)**63 - 1)
    assert (sm.Float("0", 60)**2).is_zero is True
    assert_close(sm.exp(sm.Float("1000", 60)), sp.exp(sp.Float("1000", 60)), 60)


def test_high_precision_runtime_does_not_import_sympy():
    source = r"""
import builtins
real_import = builtins.__import__
def guarded(name, *args, **kwargs):
    if name == 'sympy' or name.startswith('sympy.'):
        raise RuntimeError('SymPy runtime dependency')
    return real_import(name, *args, **kwargs)
builtins.__import__ = guarded
import symox as s
result = s.N(s.pi, 80)
assert str(result).startswith('3.141592653589793238462643383279502884197169399375105820974944')
assert result._prec > 250
assert s.Float('1.000000000000000000000000000000000001', 80) > 1
"""
    subprocess.run([sys.executable, "-c", source], check=True)
