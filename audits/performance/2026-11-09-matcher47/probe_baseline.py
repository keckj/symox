#!/usr/bin/env python3
"""Lightweight wave-47 matcher witnesses; run against the frozen baseline."""

import argparse
import importlib.util
import json
from pathlib import Path

import sympy
import symox


ROOT = Path(__file__).resolve().parents[3]
AUDIT = Path(__file__).resolve().parent


def outcome(operation, module, render):
    try:
        return {"ok": True, "value": render(module, operation())}
    except Exception as error:
        return {"ok": False, "type": type(error).__name__, "message": str(error)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", nargs="?", type=Path, default=AUDIT / "baseline-witnesses.json")
    args = parser.parse_args()
    spec = importlib.util.spec_from_file_location(
        "matcher_wave47_cases", ROOT / "tests/test_native_matcher_wave47.py"
    )
    cases_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cases_module)
    cases = next(
        mark.args[1]
        for mark in cases_module.test_noncommutative_match_is_native_and_agrees.pytestmark
        if mark.args[0] == "expr, pattern"
    )
    error_cases = [
        (f"A**({exponent})*B", "w*B")
        for exponent in ("Rational(1, 2)", "Rational(3, 2)", "Float(2)", "Float('1e-1000')")
    ]
    rows = []
    for group, entries in (("ordered_matching", cases), ("numeric_expansion_error", error_cases)):
        for expr, pattern in entries:
            row = {"group": group, "expr": expr, "pattern": pattern}
            for module in (sympy, symox):
                namespace = cases_module._namespace(module)
                target, form = eval(expr, namespace), eval(pattern, namespace)
                row[module.__name__] = outcome(
                    lambda: target.match(form), module, cases_module._render
                )
                if module is symox:
                    try:
                        handled, value = symox._core._reference_match(target, form, False)
                        row["native"] = {"handled": handled, "value": (
                            cases_module._render(symox, value) if handled else value
                        )}
                    except Exception as error:
                        row["native"] = {"type": type(error).__name__, "message": str(error)}
            row["public_agrees"] = row["sympy"] == row["symox"]
            rows.append(row)
    result = {
        "sympy_version": sympy.__version__,
        "symox_module": symox.__file__,
        "old": False,
        "counts": {
            group: {
                "cases": sum(row["group"] == group for row in rows),
                "public_agrees": sum(row["group"] == group and row["public_agrees"] for row in rows),
                "native_handled": sum(row["group"] == group and row["native"].get("handled", False) for row in rows),
            }
            for group in ("ordered_matching", "numeric_expansion_error")
        },
        "rows": rows,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result["counts"], sort_keys=True))


if __name__ == "__main__":
    main()
