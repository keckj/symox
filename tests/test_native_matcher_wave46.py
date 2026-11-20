"""Numeric inverse comparisons preserve SymPy's Float normalization protocol."""

import pytest

import symox
from symox import _core

sympy = pytest.importorskip("sympy")


def _case(module, expr, coefficient):
    ns = dict(vars(module), w=module.Wild("w"))
    return eval(expr, ns), eval(f"({coefficient})*w", ns)


def _render(module, value):
    if value is None:
        return None
    return sorted((module.srepr(k), module.srepr(v)) for k, v in value.items())


@pytest.mark.parametrize("old", [False, True])
@pytest.mark.parametrize(
    "expr, coefficient",
    [
        ("Float(2)", "Integer(2)"),
        ("Float(-2)", "Integer(-2)"),
        ("Float(3)", "Integer(2)"),
        ("Integer(2)", "Float(2)"),
        ("Float('0.1')", "Rational(1, 10)"),
        ("Float('-0.1')", "Rational(-1, 10)"),
        ("Rational(1, 10)", "Float('0.1')"),
        ("Float(0)", "Integer(2)"),
        ("Float(0, 50)", "Integer(2)"),
        # Equal values at different declared precisions must stay distinct.
        ("Float(2, 5)", "Integer(2)"),
        ("Float(2, 50)", "Integer(2)"),
        ("Integer(2)", "Float(2, 50)"),
        ("Float(2)", "Float(2, 50)"),
        ("Float(2, 50)", "Float(2)"),
        ("Float(2, precision=1)", "Integer(2)"),
        ("Integer(2)", "Float(2, precision=1)"),
        # Wide exponents must never pass through a machine double.
        ("Float('1e1000')", "Integer(10)**1000"),
        ("Float('1e1000', 15)", "Integer(10)**1000"),
        ("Float('1e-1000')", "Rational(1, 10**1000)"),
        ("Float('-1e-1000')", "Rational(-1, 10**1000)"),
        ("Rational(1, 10**1000)", "Float('1e-1000')"),
        # Reference evalf truncates exact atoms to 57 bits, then rounds
        # nearest-even to 53 bits. Both stages affect the binding.
        ("Float(1)", "1 + Rational(1, 2**53) + Rational(1, 2**57)"),
        ("Float(1) + Float(2)**-52", "1 + Rational(1, 2**53) + Rational(1, 2**57)"),
        ("Float(-1)", "-1 - Rational(1, 2**53) - Rational(1, 2**57)"),
        ("Float(1)", "1 + Rational(1, 2**53) + Rational(3, 2**58)"),
        ("Float(-1)", "-1 - Rational(1, 2**53) - Rational(3, 2**58)"),
        ("Float(2**58)", "Integer(2**58 + 35)"),
        ("Float(1)", "1 + Rational(1, 2**53) + Rational(1, 2**55)"),
    ],
)
def test_numeric_inverse_is_native_and_matches_sympy(expr, coefficient, old):
    target, pattern = _case(symox, expr, coefficient)
    ref_target, ref_pattern = _case(sympy, expr, coefficient)
    expected = _render(sympy, ref_target.match(ref_pattern, old=old))
    handled, value = _core._reference_match(target, pattern, old)
    assert handled, value
    assert _render(symox, value) == expected
    assert _render(symox, target.match(pattern, old=old)) == expected


def test_unmodelled_comparable_expression_keeps_python_fallback():
    target, pattern = _case(symox, "Float(2)", "sqrt(2)")
    ref_target, ref_pattern = _case(sympy, "Float(2)", "sqrt(2)")
    handled, reason = _core._reference_match(target, pattern, False)
    assert not handled
    assert isinstance(reason, str)
    assert _render(symox, target.match(pattern)) == _render(
        sympy, ref_target.match(ref_pattern)
    )


@pytest.mark.parametrize("old", [False, True])
def test_unevaluated_zero_inverse_keeps_reference_fallback(old):
    results = []
    for module in (sympy, symox):
        w = module.Wild("w")
        with module.evaluate(False):
            target = module.Float(0)
            pattern = module.Mul(0, w, evaluate=False)
            if module is symox:
                handled, reason = _core._reference_match(target, pattern, old)
                assert not handled
                assert isinstance(reason, str)
            results.append(_render(module, target.match(pattern, old=old)))
    # In this context, even adding exact zero remains held. The binding
    # contains 0.0/0**1, rather than normalizing the two zeros to exact 1.
    assert results[0] == results[1]


@pytest.mark.parametrize("index", [1673, 1674])
def test_float_ode_patterns_now_run_natively(index):
    # These were the two Float-only declines in the original matcher corpus.
    import gzip
    import json
    from pathlib import Path

    audit = Path(__file__).resolve().parents[1] / "audits/performance/2026-11-06-matcher"
    cases = json.loads(gzip.open(audit / "corpus.json.gz").read())["cases"]
    case = next(row for row in cases if row["id"] == index)
    expr, pattern = eval(case["expr"], vars(symox)), eval(case["pattern"], vars(symox))
    ref_expr, ref_pattern = eval(case["expr"], vars(sympy)), eval(case["pattern"], vars(sympy))
    handled, value = _core._reference_match(expr, pattern, case["old"])
    assert handled, value
    assert _render(symox, value) == _render(
        sympy, ref_expr.match(ref_pattern, old=case["old"])
    )
