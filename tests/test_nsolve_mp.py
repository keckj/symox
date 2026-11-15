"""Numerical solving must retain requested precision throughout iteration."""

import subprocess
import sys

import pytest
import sympy as sp

import symox as sm


def assert_digits(actual, expected, digits):
    numerical = sm.to_sympy(actual)
    error = abs(sp.N(numerical - expected, digits + 15))
    assert error < max(1, abs(expected)) * sp.Rational(1, 10)**(digits - 3), (actual, expected, error)


@pytest.mark.parametrize("digits", [50, 100])
@pytest.mark.parametrize("source,start", [("x**2-2", 1), ("cos(x)-x", 1), ("exp(x)-3", 1)])
def test_scalar_newton_retains_requested_precision(source, start, digits):
    x, X = sm.Symbol("x"), sp.Symbol("x")
    result = sm.nsolve(sm.sympify(source), x, start, prec=digits)
    expected = sp.nsolve(sp.sympify(source), X, start, prec=digits)
    assert result._prec == expected._prec
    assert_digits(result, expected, digits)


@pytest.mark.parametrize("digits", [15, 50, 100])
def test_complex_start_and_complex_coefficients(digits):
    x = sm.Symbol("x")
    root = sm.nsolve(x**2 + 1, x, 1 + 1j, prec=digits)
    assert_digits(root, sp.I, digits)
    # A complex equation must use the complex path even with a real start.
    assert_digits(sm.nsolve(x - sm.I, x, 0, prec=digits), sp.I, digits)
    branch = sm.nsolve(sm.exp(x) + 1, x, 3j, prec=digits)
    assert_digits(branch, sp.N(sp.I * sp.pi, digits), digits)


@pytest.mark.parametrize("overdetermined", [False, True])
def test_complex_system_and_overdetermined_hermitian_qr(overdetermined):
    x, y = sm.symbols("x y")
    equations = [x + sm.I*y - 2, x - sm.I*y - 4]
    if overdetermined:
        equations.append(2*x - 6)
    result = sm.nsolve(equations, [x, y], [0, 0], prec=80)
    assert isinstance(result, sm.Matrix)
    assert result.shape == (2, 1)
    assert_digits(result[0], sp.Integer(3), 80)
    assert_digits(result[1], sp.I, 80)


def test_nonlinear_overdetermined_system_and_dict_output():
    x, y = sm.symbols("x y")
    equations = [x*x - 2, y*y - 3, x*x + y*y - 5]
    result = sm.nsolve(equations, [x, y], [1, 1], prec=70, dict=True)
    assert len(result) == 1 and set(result[0]) == {x, y}
    assert_digits(result[0][x], sp.N(sp.sqrt(2), 70), 70)
    assert_digits(result[0][y], sp.N(sp.sqrt(3), 70), 70)
    # Default precision also supports overdetermined systems.
    result = sm.nsolve([x - 2, 2*x - 4], [x], [1])
    assert_digits(result[0], sp.Integer(2), 15)


def test_mp_bisection_respects_tolerance_and_rejects_poles():
    x = sm.Symbol("x")
    result = sm.nsolve(x*x - 2, x, (0, 2), solver="bisect", prec=70, maxsteps=400)
    assert_digits(result, sp.N(sp.sqrt(2), 70), 70)
    with pytest.raises(ValueError):
        sm.nsolve(1/x, x, (-1, 1), solver="bisect", prec=70, maxsteps=400)
    with pytest.raises(ValueError):
        sm.nsolve(x*x + 1, x, (-1, 1), solver="bisect", prec=50)
    with pytest.raises(ValueError):
        sm.nsolve(x*x - 2, x, (0, 2j), solver="bisect", prec=50)


def test_custom_mp_tolerance_iteration_budget_and_verification():
    x = sm.Symbol("x")
    result = sm.nsolve(x*x - 2, x, 1, prec=100, tol=sm.Float("1e-105", 120), maxsteps=30)
    assert_digits(result, sp.N(sp.sqrt(2), 100), 100)
    with pytest.raises(ValueError, match="Could not find root"):
        sm.nsolve(x*x - 2, x, 100, prec=70, maxsteps=1)
    unfinished = sm.nsolve(x*x - 2, x, 100, prec=70, maxsteps=1, verify=False)
    assert unfinished > 40
    assert unfinished._prec == sp.Float(1, 70)._prec
    with pytest.raises(ValueError):
        sm.nsolve([x - 1, x - 2], [x], [0], prec=50)
    for tolerance in [0, -1, float("nan"), sm.I]:
        with pytest.raises(ValueError):
            sm.nsolve(x*x - 2, x, 1, prec=50, tol=tolerance)


def test_requested_low_precision_and_starting_point_precision():
    x = sm.Symbol("x")
    start = sm.Float("1.234567890123456789", 100)
    root = sm.nsolve(x*x - 2, x, start, prec=10)
    assert root._prec == sp.Float(1, 10)._prec
    assert_digits(root, sp.N(sp.sqrt(2), 10), 10)
    with pytest.raises(ValueError):
        sm.nsolve([x - 1], [x, sm.Symbol("y")], [0, 0], prec=50)
    # Finite native MP values outside the IEEE range keep the MP path even
    # when the requested result has the default 15 decimal digits.
    large = sm.Float("1e400", 50)
    result = sm.nsolve(x - large, x, large)
    assert result.is_Float and result._prec == 53
    assert_digits(result, sp.Float("1e400", 15), 15)


@pytest.mark.parametrize("text,start,expected", [
    ("x**2 + 1", "1 + I", "1.0*I"),
    ("x**2 + 4", "1 + I", "2.0*I"),
    ("x**2 + x + 1", "I", "-0.5 + 0.866025403784439*I"),
])
def test_complex_secant_root_like_reference(text, start, expected):
    # mpmath's findroot works at prec + 20 bits: the residual imaginary or
    # real part of the iterate is below what a 53-bit result prints
    x = sm.Symbol("x")
    result = sm.nsolve(sm.sympify(text), x, sm.sympify(start))
    assert str(result) == expected
    sx = sp.Symbol("x")
    assert str(sp.nsolve(sp.sympify(text), sx, sp.sympify(start))) == expected


def test_mp_solver_has_no_sympy_or_mpmath_runtime_dependency():
    source = r"""
import builtins
original = builtins.__import__
def guarded(name, *args, **kwargs):
    if name.split('.')[0] in ('sympy', 'mpmath'):
        raise RuntimeError('external symbolic/numeric runtime')
    return original(name, *args, **kwargs)
builtins.__import__ = guarded
import symox as s
x = s.Symbol('x')
real = s.nsolve(x*x-2, x, 1, prec=80)
assert str(real).startswith('1.4142135623730950488016887242096980785696718753769480')
complex_root = s.nsolve(x*x+1, x, 1j, prec=80)
assert abs(complex(complex_root)-1j) < 1e-15
"""
    subprocess.run([sys.executable, "-c", source], check=True)
