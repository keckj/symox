"""Candidate checking is optional and assumption filtering is conservative."""

import pytest
import sympy as sp

import symox as sm


@pytest.mark.parametrize("assumptions,constant", [
    ({"positive": True}, 4), ({"negative": True}, 4),
    ({"real": True}, -1), ({"real": False}, 4),
    ({"integer": True}, 4), ({"even": True}, 4), ({"odd": True}, 4),
    ({"prime": True}, 9), ({"prime": True}, 16),
    ({"imaginary": True}, -1), ({"zero": True}, 4),
])
def test_known_assumption_contradictions_match_sympy(assumptions, constant):
    x, X = sm.Symbol("x", **assumptions), sp.Symbol("x", **assumptions)
    for check in [True, False]:
        actual = sm.solve(x*x - constant, x, check=check)
        expected = sp.solve(X*X - constant, X, check=check)
        assert set(map(sm.to_sympy, actual)) == set(expected)


def test_check_false_preserves_denominator_and_radical_candidates():
    x = sm.Symbol("x")
    assert sm.solve(sm.sin(x)/x, x) == [sm.pi]
    assert sm.solve(sm.sin(x)/x, x, check=False) == [0, sm.pi]
    assert sm.solve(sm.sqrt(x) + 1, x) == []
    assert sm.solve(sm.sqrt(x) + 1, x, check=False) == [1]
    positive = sm.Symbol("positive", positive=True)
    assert sm.solve(positive*(positive - 1), positive) == [1]
    assert sm.solve(positive*(positive - 1), positive, check=False) == [0, 1]


def test_solve_set_infers_variables_before_formatting():
    x, y = sm.symbols("x y")
    variables, rows = sm.solve(x*x - 1, set=True)
    assert variables == [x]
    assert {tuple(map(int, row)) for row in rows} == {(-1,), (1,)}
    variables, rows = sm.solve([x + y - 3, x - y - 1], set=True)
    assert variables == [x, y]
    assert {tuple(map(int, row)) for row in rows} == {(2, 1)}
    assert sm.solve(1, set=True) == ([], set())


def test_systems_propagate_candidate_checking():
    positive = sm.Symbol("positive", positive=True)
    y = sm.Symbol("y")
    assert sm.solve([positive + 1, y - 2], [positive, y]) == []
    assert sm.solve([positive + 1, y - 2], [positive, y], check=False) == {positive: -1, y: 2}
    assert sm.solve([positive*positive - 1, y - positive], [positive, y]) == [(1, 1)]
    assert sm.solve([positive*positive - 1, y - positive], [positive, y], check=False) == [(-1, -1), (1, 1)]


def test_unknown_signs_are_not_decided_by_floating_point():
    positive = sm.Symbol("positive", positive=True)
    # The current structural assumption engine cannot prove this difference
    # positive. Candidate filtering must keep an undecidable result, even
    # when a numerical approximation suggests its sign.
    candidate = sm.sqrt(2) - 2
    assert sm.solve(positive - candidate, positive) == [candidate]
