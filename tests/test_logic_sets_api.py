"""Integration and regression coverage for the logic/set public API."""

import pickle

import pytest
import sympy
import symox as sx


def test_named_set_registry_roundtrip_and_pickle():
    for name in ("EmptySet", "UniversalSet", "Reals", "Integers", "Naturals",
                 "Naturals0", "Rationals", "Complexes"):
        value = getattr(sx, name)
        assert getattr(sx.S, name) == value
        assert value.to_sympy() == getattr(sympy.S, name)
        assert pickle.loads(pickle.dumps(value)) == value
    assert sx.S.Half == sx.Rational(1, 2)
    assert not sx.EmptySet
    assert sx.FiniteSet(0)


def test_bound_symbols_do_not_leak_from_lambdas_and_conditions():
    x, y = sx.symbols("x y")
    assert sx.Lambda(x, x + y).free_symbols == {y}
    assert sx.ImageSet(sx.Lambda(x, x + y), sx.Integers).free_symbols == {y}
    assert sx.ConditionSet(x, x > y, sx.Reals).free_symbols == {y}
    # The same symbol remains free in the base set, outside the binding scope.
    assert sx.ConditionSet(x, x > 0, sx.FiniteSet(x, y)).free_symbols == {x, y}
    assert (sx.Lambda(x, x + y)._e + x).free_symbols == {x, y}


def test_product_set_preserves_nested_tuples_and_empty_product():
    a, b, c = (sx.FiniteSet(i) for i in (1, 2, 3))
    nested = sx.ProductSet(sx.ProductSet(a, b), c)
    assert nested.to_sympy() == sympy.ProductSet(
        sympy.ProductSet(sympy.FiniteSet(1), sympy.FiniteSet(2)), sympy.FiniteSet(3))
    assert nested != sx.ProductSet(a, b, c)
    assert sx.ProductSet().to_sympy() == sympy.FiniteSet(sympy.Tuple())
    assert (a ** 0) == sx.ProductSet()
    with pytest.raises(ValueError):
        a ** sx.Rational(1, 2)


def test_range_validation_and_wide_endpoint_arithmetic():
    for args in ((0, 4, 0), (0, 4, sx.Rational(1, 2)), (sx.Rational(1, 2), 4)):
        with pytest.raises(ValueError):
            sx.Range(*args)
    for args in ((-2**63, 2**63 - 1, 2), (2**63 - 1, -2**63, -2), (4, -4, -3)):
        r = sx.Range(*args)
        assert r.to_sympy() == sympy.Range(*args)
        assert r.contains(args[0]) == sx.true
    assert not sx.Range(0)
    assert sx.Range(0).is_empty is True
    assert sx.Interval(-sx.oo, 1).is_closed
    assert sx.Interval(0, sx.oo).is_closed
    assert not sx.Interval.open(0, sx.oo).is_closed


def test_periodic_solutions_exclude_denominator_poles():
    x = sx.Symbol("x")
    for pole in (sx.Integer(0), sx.pi):
        result = sx.solveset(sx.sin(x) / (x - pole), x, sx.Reals)
        assert result.contains(pole) == sx.false
        assert result.contains(2 * sx.pi) == sx.true


def test_minimization_internal_symbols_do_not_capture_user_symbols():
    x, internal = sx.symbols("x _logic_0")
    expression = (x > 0) | internal
    result = sx.simplify_logic(expression)
    assert result.free_symbols == {x, internal}
    assert sympy.simplify_logic(sympy.Xor(result.to_sympy(), expression.to_sympy())) == sympy.false


@pytest.mark.parametrize("constructor", [sx.SOPform, sx.POSform])
def test_minterm_validation(constructor):
    x = sx.Symbol("x")
    with pytest.raises(ValueError):
        constructor([x], [2])
    with pytest.raises(ValueError):
        constructor([x], [0], [0])
    with pytest.raises(ValueError):
        constructor([x, x], [0])
    with pytest.raises(ValueError):
        constructor([x], [[0, 1]])
    with pytest.raises(NotImplementedError):
        constructor(sx.symbols("x:32"), [])


def test_unsupported_solver_options_are_explicit():
    x = sx.Symbol("x")
    for kwargs in ({"use_lra_theory": True}, {"minimal": True}, {"algorithm": "z3"}):
        with pytest.raises(NotImplementedError):
            sx.satisfiable((x > 1) & (x < 0), **kwargs)
    with pytest.raises(ValueError):
        sx.satisfiable(x, algorithm="invalid")
    with pytest.raises(NotImplementedError):
        list(sx.truth_table(x, sx.symbols("x:32")))
    for clauses in ([[0]], [[2]], [[-2]]):
        with pytest.raises(ValueError):
            sx._core._sat_dimacs(1, clauses)


def test_set_constructor_validation():
    for operation in (lambda: sx.Union([1, 2]),
                      lambda: sx.Complement(sx.Reals, 1),
                      lambda: sx.ProductSet(sx.Reals, 1),
                      lambda: sx.solveset(0, sx.Symbol("x"), 1)):
        with pytest.raises(TypeError):
            operation()
