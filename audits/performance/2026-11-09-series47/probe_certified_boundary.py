"""Probe exact exponent ordering just above/below a binary64-sized integer.

Run after compiling wave 47:
    PYTHONPATH=python python probe_certified_boundary.py output.json
"""

import argparse
import hashlib
import json
from pathlib import Path

import sympy
import symox
from symox import _core


def snapshot(module, operation):
    try:
        result = operation()
        return {"status": "value", "srepr": module.srepr(result), "text": str(result)}
    except Exception as error:
        return {"status": "error", "type": type(error).__name__, "message": str(error)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    rows = []
    for offset in (-1, 1):
        numerator = 2**53 + offset
        x = symox.Symbol("x")
        alpha = symox.log(numerator)/symox.log(2)
        expr = x**alpha
        ref_x = sympy.Symbol("x")
        ref_alpha = sympy.log(numerator)/sympy.log(2)
        ref_expr = ref_x**ref_alpha
        reference = snapshot(sympy, lambda: ref_expr.series(ref_x, 0, 53))
        # The reference itself can fail above the boundary. The exact sign
        # proof supplies the mathematical oracle independently of series().
        expected_math = snapshot(sympy, lambda: ref_expr if offset < 0 else sympy.O(ref_x**53))
        native = snapshot(symox, lambda: _core._series(expr, x, 0, 53))
        public = snapshot(symox, lambda: expr.series(x, 0, 53))
        def matches(expected):
            return all(value["status"] == expected["status"] == "value"
                       and value["srepr"] == expected["srepr"]
                       for value in (native, public))
        rows.append({
            "offset": offset,
            "exponent": str(alpha),
            "exact_sign_of_exponent_minus_53": offset,
            "proof": "log is increasing: sign(log((2**53 + offset)/2**53)) = sign(offset)",
            "binary64_exponent": float(ref_alpha),
            "expected_math": expected_math,
            "reference": reference,
            "native": native,
            "public": public,
            "matched_math": matches(expected_math),
            "matched_reference": matches(reference),
        })
    report = {
        "sympy_version": sympy.__version__,
        "extension_sha256": hashlib.sha256(Path(_core.__file__).read_bytes()).hexdigest(),
        "cases": rows,
        "all_math_matched": all(row["matched_math"] for row in rows),
        "math_match_count": sum(row["matched_math"] for row in rows),
        "reference_match_count": sum(row["matched_reference"] for row in rows),
        "case_count": len(rows),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({key: report[key] for key in
                      ("all_math_matched", "math_match_count", "reference_match_count", "case_count")}
                     | {"output": str(args.output)}))
    return 0 if report["all_math_matched"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
