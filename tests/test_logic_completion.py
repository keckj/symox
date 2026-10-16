"""Care-set minimisation and sparse/dense algebraic normal forms."""

from itertools import product
import pickle

import pytest
import sympy as sp
from sympy.logic import boolalg as oracle
import symox as sx


def assignment_rows(variables):
    for values in product((False, True), repeat=len(variables)):
        yield dict(zip(variables, values))


def assert_equivalent(a, b, variables):
    for row in assignment_rows(variables):
        assert a.xreplace(row) == b.xreplace(row), (a, b, row)


def test_dontcare_priority_and_only_dontcare_variables():
    x, y, z = sx.symbols("x y z")
    assert sx.simplify_logic(x | y, dontcare=y) == x
    assert sx.simplify_logic(x | y, dontcare=x | z) == y
    assert sx.simplify_logic(x | y, dontcare=True) == sx.false
    assert sx.simplify_logic(x | y, dontcare=False) == x | y
    assert sx.simplify_logic(x | y, dontcare=z) == x | y
    assert sx.simplify_logic(x | y, dontcare=(x & ~y) | (~x & y)) in (x, y)
    assert sx.simplify_logic(True, dontcare=True) == sx.true
    assert sx.simplify_logic(False, dontcare=False) == sx.false


@pytest.mark.parametrize("form", [None, "dnf", "cnf"])
def test_dontcare_care_agreement_and_sympy_differential(form):
    x, y, z, w = sx.symbols("x y z w")
    expressions = [x | y, sx.ITE(x, y, z), (x | y) & (~x | z),
                   sx.Equivalent(x, y, z), x ^ y ^ z]
    dontcares = [x, y | w, (x & y) | (~x & z), sx.true, sx.false]
    for e, dc in product(expressions, dontcares):
        result = sx.simplify_logic(e, form=form, dontcare=dc)
        expected = oracle.simplify_logic(e.to_sympy(), form=form, dontcare=dc.to_sympy())
        # Several minimum covers can tie; compare the specified care set.
        for row in assignment_rows([x, y, z, w]):
            if dc.xreplace(row) == sx.false:
                assert result.xreplace(row) == e.xreplace(row)
        assert sp.simplify_logic((result.to_sympy() ^ expected) & ~dc.to_sympy()) == sp.false
        if form == "dnf":
            assert sx.is_dnf(result)
        if form == "cnf":
            assert sx.is_cnf(result)


def test_dontcare_relational_complements_and_capture():
    x, internal, y = sx.symbols("x _logic_0 y")
    e = (x > 0) | internal
    result = sx.simplify_logic(e, dontcare=x <= 0)
    assert result == sx.true
    e = ((x > 0) & y) | internal
    dc = (x <= 0) & internal
    r = sx.simplify_logic(e, dontcare=dc)
    assert sp.simplify_logic((r.to_sympy() ^ e.to_sympy()) & ~dc.to_sympy()) == sp.false
    assert r.free_symbols <= {x, internal, y}


def test_dontcare_threshold_and_argument_validation():
    xs = sx.symbols("x:10")
    e = (xs[0] | xs[1]) & (~xs[0] | xs[2])
    # Too many variables introduced by dontcare: preserve the existing
    # bounded default minimisation of the original expression.
    dc = sx.And(*xs[3:])
    assert sx.simplify_logic(e, dontcare=dc) == sx.simplify_logic(e)
    with pytest.raises(ValueError):
        sx.simplify_logic(e, form="anf", dontcare=xs[0])
    with pytest.raises(NotImplementedError):
        sx.simplify_logic(e, dontcare=sx.And(*sx.symbols("q:21")), force=True)


def test_anf_all_three_variable_truth_tables_and_coefficients():
    variables = sx.symbols("x y z")
    for table in product((0, 1), repeat=8):
        coeffs = sx.anf_coeffs(table)
        assert coeffs == oracle.anf_coeffs(table)
        assert sx.anf_coeffs(coeffs) == list(table)
        expression = sx.ANFform(variables, table)
        assert sx.is_anf(expression)
        for index, row in enumerate(assignment_rows(variables)):
            assert expression.xreplace(row) == sx.sympify(bool(table[index]))
        expected = oracle.ANFform([v.to_sympy() for v in variables], table)
        assert expression.to_sympy() == expected


def test_anf_symbolic_nested_connectives():
    x, y, z = sx.symbols("x y z")
    expressions = [~(x | y), sx.Implies(~x, sx.Equivalent(y, z)),
                   sx.Equivalent(~x, y | z), ~x ^ (y & z),
                   (x & y) | ~z, sx.ITE(x, y | z, ~y)]
    for e in expressions:
        a = sx.to_anf(e)
        assert sx.is_anf(a)
        assert_equivalent(a, e, [x, y, z])
        expected = oracle.to_anf(e.to_sympy())
        if oracle.is_anf(expected):
            assert a.to_sympy() == expected
        assert sx.to_anf(a) == a


def test_anf_shallow_conversion_preserves_opaque_arguments():
    x, y, z = sx.symbols("x y z")
    e = sx.Implies(~x, sx.Equivalent(y, z))
    result = sx.to_anf(e, deep=False)
    assert result.to_sympy() == oracle.to_anf(e.to_sympy(), deep=False)
    for e in [x & (y | z), ~x ^ (y & z), sx.ITE(x, y, z)]:
        assert sx.to_anf(e, deep=False) == e
    assert sx.to_anf(~(x | y), deep=False).to_sympy() == oracle.to_anf(~(sp.Symbol("x") | sp.Symbol("y")), deep=False)


def test_anf_constants_relations_and_more_than_twenty_variables():
    x, y = sx.symbols("x y")
    for e in [sx.true, sx.false, x, x > 0, y & (x > 0), y ^ (x > 0)]:
        assert sx.is_anf(e) == oracle.is_anf(e.to_sympy())
        assert sx.to_anf(e).to_sympy() == oracle.to_anf(e.to_sympy())
    variables = sx.symbols("v:40")
    e = sx.Implies(sx.Xor(*variables), x)
    a = sx.to_anf(e)
    assert sx.is_anf(a)
    assert len(a.args) == 81
    assert a.free_symbols == set(variables) | {x}
    assert sx.ANFform([], [0]) == sx.false
    assert sx.ANFform([], [1]) == sx.true
    assert sx.ANFform([x], [1, 0]).to_sympy() == oracle.ANFform([sp.Symbol("x")], [1, 0])


def test_anf_validations_and_exact_integer_xor():
    x = sx.Symbol("x")
    for values in [[], [0, 1, 0]]:
        with pytest.raises(ValueError):
            sx.anf_coeffs(values)
    for values in [[], [0], [0, 1, 0]]:
        with pytest.raises(ValueError):
            sx.ANFform([x], values)
    for values in [[0.0, 1], ["a", 0]]:
        with pytest.raises(TypeError):
            sx.anf_coeffs(values)
    for values in [[True, False], [2, 0], [2**100, 2**80]]:
        assert sx.anf_coeffs(values) == oracle.anf_coeffs(values)
        assert sx.ANFform([x], values).to_sympy() == oracle.ANFform([sp.Symbol("x")], values)
    with pytest.raises(NotImplementedError):
        sx.ANFform(sx.symbols("q:21"), [])


def test_anf_constant_term_survives_pickle_and_sympy_conversion():
    x, y = sx.symbols("x y")
    expression = sx.to_anf(~(x | y))
    restored = pickle.loads(pickle.dumps(expression))
    assert restored == expression
    assert sx.is_anf(restored)
    assert sx.true in restored.args
    assert oracle.is_anf(expression.to_sympy())
    assert expression.to_sympy() == oracle.to_anf(~(sp.Symbol("x") | sp.Symbol("y")))
    assert sx.sympify(expression.to_sympy()) == expression
