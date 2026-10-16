"""Expression API contracts needed by Trisplyne's symbolic code generation."""

import pickle

import numpy as np
import pytest
import sympy as sp

import symox as sm


@pytest.mark.parametrize("constructor,value,predicate", [
    (sm.Integer, 3, "is_Integer"),
    (sm.Rational, "2/3", "is_Rational"),
    (sm.Float, 0.5, "is_Float"),
    (sm.Symbol, "x", "is_Symbol"),
])
def test_constructor_classes(constructor, value, predicate):
    result = constructor(value)
    assert isinstance(result, sm.Expr)
    assert isinstance(result, constructor)
    assert getattr(result, predicate)
    assert result.func is constructor


def test_arithmetic_classes_and_reconstruction():
    x, y = sm.symbols("x y")
    for expression, constructor in [(x + y, sm.Add), (x * y, sm.Mul), (x**2, sm.Pow)]:
        assert isinstance(expression, constructor)
        assert expression.func is constructor
        assert expression.func(*expression.args) == expression
    assert sm.Add() == 0
    assert sm.Mul() == 1
    assert sm.Add.make_args(x) == (x,)
    assert sm.Mul.make_args(x) == (x,)
    assert set(sm.Add.make_args(x + y)) == {x, y}
    assert set(sm.Mul.make_args(x * y)) == {x, y}
    assert sm.S.Half == sm.Rational(1, 2)
    assert sm.S(2) == sm.Integer(2)


@pytest.mark.parametrize("source", [
    "x**2*y", "2*x+y+1", "x**2+y**2+x*y", "x*(x+y)",
    "x*sin(y)+x*cos(y)", "3*x**2*y+2*x*z+1",
])
def test_args_order_matches_sympy(source):
    actual, expected = sm.sympify(source), sp.sympify(source)
    assert [str(arg) for arg in actual.args] == [str(arg) for arg in expected.args]


def test_held_expression_retains_horner_coefficient():
    x, a, b = sm.symbols("x a b")
    coefficient = sm.UnevaluatedExpr(a + b)
    assert coefficient.args == (a + b,)
    assert coefficient.func is sm.UnevaluatedExpr
    assert isinstance(coefficient, sm.UnevaluatedExpr)
    expression = 2 * coefficient + x
    assert coefficient in expression.args[1].args or coefficient in expression.args[0].args
    held = sm.Mul(-1, a + b, evaluate=False)
    assert held.args == (sm.Integer(-1), a + b)
    assert isinstance(held, sm.Mul)
    assert held != -a - b


def test_pickle_preserves_held_tree_and_assumptions():
    x = sm.Symbol("x", real=True, nonzero=True)
    quoted = sm.Symbol("a'b\\c", positive=True)
    expression = sm.Mul(-1, sm.UnevaluatedExpr(x + quoted), evaluate=False)
    restored = pickle.loads(pickle.dumps(expression))
    assert restored == expression
    assert restored.args == expression.args
    assert restored.free_symbols == {x, quoted}
    assert all(symbol.is_real is True for symbol in restored.free_symbols)


def test_coefficient_power_and_two_term_accessors():
    x, y = sm.symbols("x y")
    assert (3 * x**2 * y).as_coeff_exponent(x) == (3 * y, sm.Integer(2))
    assert (x**-2).as_coeff_exponent(x) == (sm.Integer(1), sm.Integer(-2))
    assert (x + y).as_coeff_exponent(x) == (x + y, sm.Integer(0))
    assert (x**2).base == x
    assert (x**2).exp == 2
    for expression in [1 + x + y, 3 * x * y]:
        first, rest = expression.as_two_terms()
        assert expression.func(first, rest) == expression
        assert first == expression.args[0]
    assert x.is_commutative is True


@pytest.mark.parametrize("source", [
    "x", "x+y+z", "2*x*y", "x**3", "x/y", "-2*x", "-x-y", "1/2", "x**(-2)",
])
def test_operation_counts(source):
    assert sm.sympify(source).count_ops() == sp.count_ops(sp.sympify(source))


def test_sort_key_and_numpy_object_interop():
    x, y, z = sm.symbols("x y z")
    expressions = [z, x**2, 2 * y, sm.Integer(1), x]
    expected = sorted(map(str, expressions), key=lambda s: sp.default_sort_key(sp.sympify(s)))
    assert list(map(str, sorted(expressions, key=sm.default_sort_key))) == expected
    a = np.array([x, y], dtype=object)
    assert np.dot(a, a) == x**2 + y**2
    assert list(a * 2 + 1) == [2 * x + 1, 2 * y + 1]


def test_lambdify_held_coefficients():
    x, y = sm.symbols("x y")
    expression = x * sm.UnevaluatedExpr(2 * y + 1)
    compiled = sm.lambdify((x, y), expression, "numpy")
    assert compiled(3, 2) == 15


def test_unevaluated_printing_matches_sympy():
    x, y = sm.symbols("x y")
    X, Y = sp.symbols("x y")
    actual, expected = sm.UnevaluatedExpr(x + y), sp.UnevaluatedExpr(X + Y)
    for a, b in [(actual, expected), (2 * actual, 2 * expected), (actual**2, expected**2)]:
        assert str(a) == str(b)
        assert sm.srepr(a) == sp.srepr(b)
        assert a.is_commutative == b.is_commutative
    assert sm.sympify(sm.srepr(actual)) == actual


def test_sympy_interop_preserves_held_arithmetic():
    x, y = sm.symbols("x y", real=True)
    expressions = [
        sm.UnevaluatedExpr(x + y),
        sm.Mul(-1, x + y, evaluate=False),
        sm.Mul(y, x, evaluate=False),
        sm.Pow(x, 1, evaluate=False),
        sm.Mul(2, sm.UnevaluatedExpr(x + y), evaluate=False),
        sm.sin(sm.UnevaluatedExpr(x + y)),
    ]
    for expression in expressions:
        converted = sm.to_sympy(expression)
        assert sm.from_sympy(converted) == expression
        assert expression._sympy_() == converted
    held = sm.to_sympy(expressions[1])
    assert held.func is sp.Mul
    assert held.args == (-1, sp.Symbol("x", real=True) + sp.Symbol("y", real=True))
    assert sm.sympify("Mul(-1, x + y, evaluate=False)").args == (-1, sm.Symbol("x") + sm.Symbol("y"))
