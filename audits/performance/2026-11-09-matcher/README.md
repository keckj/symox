# Native reference matcher: numeric inverse comparisons (wave 46)

The Rust reference matcher now implements the numerical-atom branch of
`Mul._combine_inverse` from SymPy 1.14.0. Previously, it declined whenever a
Float faced any expression with no free symbols. The public API then ran the
existing Python reference protocol. This change moves the comparison of
Integer, Rational and Float atoms into Rust, without changing that fallback
for other symbol-free expressions.

## Algorithm and precision

The test is the reference's
`left.__add__(0) == right.evalf().__add__(0)`, in each applicable direction.
It is not approximate equality. The implementation:

- preserves the declared precision of the left Float;
- truncates an exact numeric atom toward zero at 57 bits, then rounds to
  53 bits (nearest, ties to even), following mpmath's `round_fast` default
  in SymPy's `evalf_integer` / `evalf_rational`;
- rounds a right Float's stored value to 53 bits;
- normalizes floating zeros by adding exact zero on both sides;
- preserves the reference's undecided comparability of one-bit Floats;
- retains arbitrary exponents, with no conversion through `f64`;
- returns the exact integer 1 only when this normalization test succeeds.

When the test is false, the existing power cancellation and division code
continues. A nonnumeric symbol-free operand still causes a native decline,
and the public Python reference matcher remains available. The Python
binding also declines under `evaluate(False)`, because native evaluated
arithmetic cannot reproduce held numerical normalization in that context.
Infinities and the other unsupported node classes retain their previous
admission rules.

## Evidence and regression coverage

Baseline replay on the installed wave-45 extension, with SymPy 1.14.0:

- 2,311 corpus cases: 2,226 native results, 85 declines;
- 83 declines report `unsupported node class`;
- the remaining two report `Float in a multiplicative inverse` (case IDs
  1673 and 1674, ODE patterns with coefficients `0.0020833333333333333` and
  `-9.8000000000000007`).

`tests/test_native_matcher_wave46.py` adds 61 regression cases: 28 numeric
pairs under both `old` flags, a preserved fallback, and the two formerly
declined ODE cases, plus the held zero-inverse fallback under both flags.
The pairs include low and high Float precisions, one-bit
Floats, floating zero, rational coefficients, exponents +/-1000, and values
near a half ulp. Each native result and public result is compared with
SymPy, including Float precision in `srepr`.

The original 22 pairs were constructed on the existing extension and SymPy; both
expression trees agree before matching. The new Rust unit test independently
checks the half-ulp normalization, declared precision, zero normalization,
tiny rationals, and the preserved decline for `pi`. A separate lightweight
integer-division model of the 57-bit truncation / 53-bit rounding was compared
with SymPy 1.14 on 2,003 exact numeric atoms (seed 46, numerator/denominator
sizes up to 1,200 bits), with zero differences.

Reproduce that probe independently of the extension:

```
/home/keckj/symox-ci/venv/bin/python audits/performance/2026-11-09-matcher/probe_rounding.py
```

The audit-only [script](probe_rounding.py) implements both stages with
integer arithmetic, including final rounding without libmp. Its recorded
[JSON result](rounding-results.json) contains the seed, generated bit-size
range, fixed tiny/large cases, observed sizes, and counts (2,003 equal,
zero different). The random cases request up to 1,200 bits; the three fixed
cases additionally include the 3,322-bit denominator `10**1000`.

## Final integrated validation

The integrating agent compiled the final extension and replayed the complete
2,311-case corpus against the recorded SymPy 1.14 reference:

| Implementation | Before: equal / different / declined | After: equal / different / declined |
| --- | --- | --- |
| Native reference matcher | 2226 / 0 / 85 | 2228 / 0 / 83 |
| Public `Expr.match` | 2311 / 0 / 0 | 2311 / 0 / 0 |
| Existing Python protocol alone | 2309 / 2 / 0 | 2309 / 2 / 0 |
| Older `pattern_match` | 2059 / 252 / 0 | 2059 / 252 / 0 |

Only native cases **1673** and **1674** changed: both formerly declined,
both now return the same proven no-match result (`None`) as SymPy. They
belong to the ODE source, whose native counts become 1,546 equal and 17
declined (before: 1,544 and 19). The other 2,309 native rows are unchanged,
including all 83 remaining declines. Every public/Python-protocol/older-matcher
row is unchanged; no binding-order difference was observed for equal results.

Receipts: [before report](../2026-11-09-wave46/raw/matcher-report-before.json),
[after report](../2026-11-09-wave46/raw/matcher-report-after.json),
[before rows](../2026-11-09-wave46/raw/matcher-before.json),
[after rows](../2026-11-09-wave46/raw/matcher-after.json), and
[SymPy rows](../2026-11-09-wave46/raw/matcher-sympy.json).

The final [focused validation](../2026-11-09-wave46/raw/focused.log) passed
318 tests, including this lot's 61 regressions and the existing matcher
corpus tests ([command and extension identity](../2026-11-09-wave46/raw/focused.command.json)).
The integrated [Rust tests](../2026-11-09-wave46/raw/cargo-test.log) passed
879 unit tests (one ignored), plus 1 integration test and 5 doctests; the
new Float inverse unit test passed. The
[correctness clippy run](../2026-11-09-wave46/raw/clippy.log) completed
successfully. The final [complete suite receipt](../2026-11-09-wave46/raw/suite-run3-summary.json)
records 61,384 passed, 14 xfailed and 41 skipped; all 61,439 collected tests
finished, with no missing files, failures or run errors. Compilation and
heavy validation were centralized by the integrating agent.

This establishes the recorded numeric matcher behavior for the tested
domains. It makes no speedup or global behavioral-parity claim.

## Remaining gaps

- Comparisons of a Float with a nonnumeric symbol-free expression, including
  `sqrt(2)` and `pi`, still run the Python reference protocol.
- Noncommutative expressions, nonfinite atoms, binders, matrices, sets,
  logical/relational/Piecewise nodes, custom host classes, symbol-order
  dependent collection and step-budget exhaustion retain their existing
  declines. This lot makes no global parity claim.
- A separate public numerical-evaluation convention difference was observed:
  for `r = 1 + 2**-53 + 2**-57`, SymPy `r.evalf()` returns the 53-bit Float
  `1.0`, while Symox returns `1.0000000000000002`. Direct correct rounding
  to 53 bits gives the latter; the reference's intermediate 57-bit rounding
  (truncation toward zero) gives the former. The same convention affects
  `Rational(1, 10**1000).evalf()`, which is one ulp below the directly rounded
  `Float('1e-1000')` in SymPy. This is a reference-protocol/double-rounding divergence,
  not a conclusion that the directly rounded value is mathematically wrong.
  The matcher follows SymPy's protocol locally. The general evaluator is
  unchanged and needs a separate scope/contract decision and audit.
