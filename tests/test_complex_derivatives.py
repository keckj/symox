"""Absolute value derivatives distinguish real and complex directions."""
import pickle

import numpy as np
import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize("flags", [{}, {"real": True}, {"imaginary": True},
                                    {"positive": True}, {"negative": True}])
@pytest.mark.parametrize("source", ["Abs(x)", "Abs(2*x)", "Abs(I*x)",
                                    "Abs(x+I)", "Abs(x+2)", "Abs(exp(x))",
                                    "re(x)", "im(x)", "conjugate(x)",
                                    "re(x**2+I*x)", "im(x**2+I*x)",
                                    "conjugate(x**2+I*x)"])
def test_directional_derivatives_against_sympy(flags, source):
    x, X = sx.Symbol("x", **flags), sp.Symbol("x", **flags)
    expression = sx.sympify(source).xreplace({sx.Symbol("x"): x})
    reference = sp.sympify(source).xreplace({sp.Symbol("x"): X})
    actual, expected = expression.diff(x), reference.diff(X)
    assert sp.simplify(actual.to_sympy()-expected) == 0, (flags, source, actual, expected)
    assert pickle.loads(pickle.dumps(actual)) == actual


def test_unconstrained_abs_retains_complex_derivatives_and_exact_public_tree():
    x = sx.Symbol("x")
    X = sp.Symbol("x")
    actual, expected = sx.Abs(x).diff(x), sp.Abs(X).diff(X)
    assert actual.to_sympy() == expected
    assert actual.has(sx.Derivative(sx.re(x), x, evaluate=False))
    assert actual.has(sx.Derivative(sx.im(x), x, evaluate=False))
    assert actual != sx.sign(x)
    assert sx.sympify(actual.srepr()) == actual


@pytest.mark.parametrize("source", ["Abs(x+I)", "Abs((1+I)*x+2)",
                                    "Abs(exp(I*x)+2)", "Abs(sin(x)+I*cos(2*x))"])
def test_real_parameter_abs_derivative_matches_real_finite_differences(source):
    x = sx.Symbol("x", real=True)
    expression = sx.sympify(source).xreplace({sx.Symbol("x"): x})
    derivative = expression.diff(x)
    assert not derivative.has(sx.Derivative(sx.re(x), x, evaluate=False))
    values = np.asarray([-.7, -.2, .3, .9])
    function_code = compile(sx.numpycode(expression), "<absolute-value>", "eval")
    gradient_code = compile(sx.numpycode(derivative), "<absolute-derivative>", "eval")
    function = lambda values: eval(function_code, {"numpy": np, "x": values})
    gradient = lambda values: eval(gradient_code, {"numpy": np, "x": values})
    epsilon = 1e-6
    numerical = (function(values+epsilon)-function(values-epsilon))/(2*epsilon)
    np.testing.assert_allclose(gradient(values), numerical, rtol=2e-6, atol=1e-7)


def test_imaginary_symbol_projections_and_conjugation():
    x = sx.Symbol("x", imaginary=True)
    assert sx.re(x) == 0
    assert sx.im(x) == -sx.I*x
    assert sx.conjugate(x) == -x
    assert sx.re(x).diff(x) == 0
    assert sx.im(x).diff(x) == -sx.I
    assert sx.conjugate(x).diff(x) == -1
    assert sx.Abs(x).diff(x) == -sx.sign(x)


def test_real_abs_derivative_stays_codegen_ready():
    x = sx.Symbol("x", real=True)
    assert sx.Abs(x).diff(x) == sx.sign(x)
    assert sx.Abs(x*x+1).diff(x) == 2*x
