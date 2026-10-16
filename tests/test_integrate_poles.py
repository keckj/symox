"""Ordinary improper integrals must never turn into finite principal values."""

import pytest
import sympy as sp
import symox as sm


@pytest.mark.parametrize("integrand,a,b,answer", [
    ("1/x", "-1", "1", "nan"),
    ("1/x**2", "-1", "1", "oo"),
    ("-1/x**2", "-1", "1", "-oo"),
    ("1/(x-1)**3", "0", "2", "nan"),
    ("1/(x-1)**2", "0", "2", "oo"),
    ("(x+Rational(1,2))/(x**2-1)**2", "-2", "2", "nan"),
    ("1/(x**2-2)", "-2", "2", "nan"),
    ("1/(x**2-2)**2", "-2", "2", "oo"),
    ("-1/(x**2-2)**2", "-2", "2", "-oo"),
    ("1/(x-sqrt(2))", "0", "2", "nan"),
    ("1/(x-sqrt(2))**2", "0", "2", "oo"),
    ("1/(x-sqrt(2))", "sqrt(2)", "2", "oo"),
    ("1/(x-sqrt(2))", "0", "sqrt(2)", "-oo"),
    ("1/(x-sqrt(2))", "sqrt(2)+1", "sqrt(2)", "-oo"),
    ("1/(x**2-2)", "0", "sqrt(2)", "-oo"),
    ("1/(x**2-2)", "sqrt(2)", "2", "oo"),
    ("1/x", "0", "1", "oo"),
    ("1/x", "-1", "0", "-oo"),
    ("1/x", "1", "0", "-oo"),
    ("1/x", "0", "-1", "oo"),
    ("1/x**2", "0", "1", "oo"),
    ("1/x**2", "-1", "0", "oo"),
    ("1/x**2", "1", "0", "-oo"),
    ("1/x**2-1", "0", "oo", "nan"),
    ("1/x", "-oo", "-1", "-oo"),
    ("1/x", "1", "oo", "oo"),
    ("1/x", "-oo", "oo", "nan"),
    ("x", "-oo", "oo", "nan"),
    ("x**2", "-oo", "oo", "oo"),
    ("-x**2", "-oo", "oo", "-oo"),
    ("1/x**2", "-oo", "oo", "oo"),
])
def test_real_one_sided_pole_and_tail_divergences(integrand, a, b, answer):
    x = sm.Symbol("x")
    result = sm.integrate(sm.sympify(integrand), (x, sm.sympify(a), sm.sympify(b)))
    assert result == sm.sympify(answer)
    # SymPy can retain a finite imaginary logarithm branch beside infinity;
    # its real part agrees with the ordinary real improper integral.
    sx = sp.Symbol("x")
    expected = sp.integrate(sp.sympify(integrand), (sx, sp.sympify(a), sp.sympify(b)))
    expected = sp.simplify(sp.re(expected)) if expected is not sp.nan else sp.nan
    assert expected == sp.sympify(answer)


@pytest.mark.parametrize("integrand,a,b", [
    ("(x**2-1)/(x-1)", "-2", "2"),
    ("(x**2-1)/(x-1)", "2", "-2"),
    ("(x-1)**2/(x**2-2*x+1)", "0", "2"),
    ("1/x", "1", "2"),
    ("1/x", "-2", "-1"),
    ("1/(x**2+1)", "-oo", "oo"),
    ("1/x**2", "1", "oo"),
    ("1/x**2", "-oo", "-1"),
    ("1/(x**2+1)", "0", "sqrt(2)"),
    ("1/(x**2+1)", "sqrt(2)", "0"),
])
def test_removable_poles_and_regular_rational_integrals(integrand, a, b):
    x = sm.Symbol("x")
    result = sm.integrate(sm.sympify(integrand), (x, sm.sympify(a), sm.sympify(b)))
    sx = sp.Symbol("x")
    expected = sp.integrate(sp.sympify(integrand), (sx, sp.sympify(a), sp.sympify(b)))
    assert sp.simplify(sp.sympify(str(result))-expected) == 0


def test_unresolved_real_poles_do_not_produce_a_finite_value():
    x, a = sm.symbols("x a")
    for f, lower, upper in [(1/(x**5-x+1), -2, 2), (1/(x-a), -1, 1)]:
        result = sm.integrate(f, (x, lower, upper))
        assert "Integral" in str(result)
    assert sm.integrate(1/x**2, (x, 0, 0)) == 0


def test_pole_ordering_is_certified_beyond_double_precision():
    x = sm.Symbol("x")
    lower = sm.Rational(1414213562373095048801688724209, 10**30)
    upper = sm.Rational(1414213562373095048801688724210, 10**30)
    assert float(lower) == float(upper)
    assert sm.integrate(1/(x*x-2), (x, lower, upper)) == sm.nan
