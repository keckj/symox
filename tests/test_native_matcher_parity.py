"""Differential parity of the native reference matcher (symox::matches).

The corpus (audits/performance/2026-11-06-matcher) holds 2311 unique
``expr.match(pattern)`` calls harvested from the SymPy 1.14 ODE solvers,
SymPy's test_match.py and seeded random variations, with SymPy's answers.
The native port must either reproduce SymPy's answer or decline; it must
never return a different match.  The public ``Expr.match`` must agree with
SymPy on every case it serves natively.
"""
import gzip
import json
import pathlib

import pytest

import symox
from symox import _core, srepr

ROOT = pathlib.Path(__file__).resolve().parents[1]
AUDIT = ROOT / "audits" / "performance" / "2026-11-06-matcher"


def _load():
    cases = json.loads(gzip.open(AUDIT / "corpus.json.gz").read())["cases"]
    ref = {r["id"]: r for r in json.loads(gzip.open(AUDIT / "sympy-results.json.gz").read())}
    return cases, ref


def _render(value):
    if value is None:
        return None
    return sorted([srepr(k), srepr(v)] for k, v in value.items())


def _expected(row):
    # symox's WildFunction is a Wild subclass whose srepr is ``Wild('g')``;
    # the binding itself is the same object (representation only).
    assert row["ok"]
    if row["value"] is None:
        return None
    return sorted([k.replace("WildFunction(", "Wild("), v] for k, v in row["value"])


def test_corpus_native_declines_or_agrees():
    cases, ref = _load()
    ns = vars(symox)
    served = different = 0
    wrong = []
    for case in cases:
        expr, pattern = eval(case["expr"], ns), eval(case["pattern"], ns)
        handled, value = _core._reference_match(expr, pattern, case["old"])
        if not handled:
            continue
        served += 1
        if _render(value) != _expected(ref[case["id"]]):
            different += 1
            wrong.append(case["id"])
    assert different == 0, wrong[:20]
    # The native port serves the overwhelming majority of the corpus.
    assert served >= 2200, served


def test_corpus_public_match_agrees_where_native():
    cases, ref = _load()
    ns = vars(symox)
    wrong = []
    for case in cases:
        expr, pattern = eval(case["expr"], ns), eval(case["pattern"], ns)
        handled, _ = _core._reference_match(expr, pattern, case["old"])
        if handled and _render(expr.match(pattern, old=case["old"])) != _expected(ref[case["id"]]):
            wrong.append(case["id"])
    assert not wrong, wrong[:20]


def test_linear_pattern_coefficients():
    x = symox.Symbol("x")
    a, b = symox.Wild("a", exclude=[x]), symox.Wild("b", exclude=[x])
    assert (2 * x + 3).match(a * x + b) == {a: 2, b: 3}
    handled, value = _core._reference_match(2 * x + 3, a * x + b, False)
    assert handled and value == {a: 2, b: 3}


def test_power_pattern_choice():
    x = symox.Symbol("x")
    a, c = symox.Wild("a4", exclude=[x]), symox.Wild("c4", exclude=[x])
    handled, value = _core._reference_match(x**2, c**2 * x**(2 * a), False)
    assert handled and value == {a: 1, c: 1}


def test_common_integer_power_of_product():
    # Mul.as_base_exp: z**2*cos(y)**2 is (z*cos(y))**2.
    x, y, z = symox.symbols("x y z")
    d = symox.Wild("d", exclude=[x])
    expr = z**2 * symox.cos(y)**2
    handled, value = _core._reference_match(expr, d**2, False)
    assert handled and value == {d: z * symox.cos(y)}
    assert expr.match(d**2) == {d: z * symox.cos(y)}


def test_no_match_is_none():
    x, y = symox.symbols("x y")
    a = symox.Wild("a", exclude=[x])
    handled, value = _core._reference_match(y, a * x, False)
    assert handled and value is None


def test_noncommutative_declines():
    A, B = symox.symbols("A B", commutative=False)
    w = symox.Wild("w")
    handled, reason = _core._reference_match(A * B, w * B, False)
    assert not handled and isinstance(reason, str)


def test_wild_property_error_propagates():
    x = symox.Symbol("x")

    def boom(_):
        raise ValueError("property")

    w = symox.Wild("w", properties=[boom])
    with pytest.raises(ValueError):
        _core._reference_match(x + 1, w + 1, False)
    with pytest.raises(ValueError):
        (x + 1).match(w + 1)
