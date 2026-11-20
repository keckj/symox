#!/usr/bin/env python3
"""Audit only: model the matcher's exact-atom evalf convention with integers.

No Symox extension is imported. Compare with SymPy 1.14's default evalf.
"""

import argparse
import json
from pathlib import Path
import random

import sympy


def model_evalf_tuple(value):
    """57-bit truncation followed by independent 53-bit nearest-even rounding."""
    numerator, denominator = abs(int(value.p)), int(value.q)
    shift = 57 + denominator.bit_length() - numerator.bit_length()
    if shift >= 0:
        mantissa = (numerator << shift) // denominator
    else:
        mantissa = numerator // (denominator << -shift)
    excess = max(mantissa.bit_length() - 57, 0)
    mantissa >>= excess
    exponent = -shift + excess

    # Final rounding is implemented directly, independently of libmp.
    excess = max(mantissa.bit_length() - 53, 0)
    if excess:
        retained, discarded = divmod(mantissa, 1 << excess)
        half = 1 << (excess - 1)
        if discarded > half or (discarded == half and retained % 2):
            retained += 1
        mantissa = retained
        exponent += excess
    if not mantissa:
        return (0, 0, 0, 0)
    trailing = (mantissa & -mantissa).bit_length() - 1
    mantissa >>= trailing
    exponent += trailing
    return (int(int(value.p) < 0), mantissa, exponent, mantissa.bit_length())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path,
        default=Path(__file__).with_name("rounding-results.json"),
    )
    args = parser.parse_args()
    seed, random_count, max_bits = 46, 2000, 1200
    rng = random.Random(seed)
    values = [
        sympy.Rational(1, 10**1000),
        sympy.Rational(-1, 10**1000),
        sympy.Integer(2**58 + 35),
    ]
    for _ in range(random_count):
        numerator = rng.getrandbits(rng.randint(1, max_bits)) * (-1)**rng.randrange(2)
        denominator = max(1, rng.getrandbits(rng.randint(1, max_bits)))
        values.append(sympy.Rational(numerator, denominator))

    failures = []
    for index, value in enumerate(values):
        actual = model_evalf_tuple(value)
        evaluated = value.evalf()
        expected = evaluated._mpf_ if evaluated else (0, 0, 0, 0)
        if actual != expected:
            failures.append({
                "index": index, "value": str(value),
                "model": actual, "sympy": expected,
            })
    result = {
        "sympy_version": sympy.__version__,
        "seed": seed,
        "random_cases": random_count,
        "random_requested_bit_range": [1, max_bits],
        "fixed_cases": [str(value) for value in values[:3]],
        "total": len(values),
        "equal": len(values) - len(failures),
        "different": len(failures),
        "observed_max_numerator_bits": max(abs(int(v.p)).bit_length() for v in values),
        "observed_max_denominator_bits": max(int(v.q).bit_length() for v in values),
        "intermediate_rounding": "toward zero, 57 bits",
        "final_rounding": "nearest, ties to even, 53 bits",
        "failures": failures[:20],
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"{result['equal']}/{result['total']} equal; {result['different']} different")
    return bool(failures)


if __name__ == "__main__":
    raise SystemExit(main())
