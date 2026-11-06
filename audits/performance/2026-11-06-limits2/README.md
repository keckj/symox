# Limits and series against SymPy 1.14, 2026-11-06 (`feat-limits2-20261106`)

Wave 43 lot. Every assert of the SymPy 1.14 test files for limits and
series is replayed as a statement (`scripts/tl.py`, driven by
`scripts/run_tl.sh`). Each statement is evaluated by SymPy and by symox, and
the printed results are compared: `scripts/compare.py`, `scripts/classify.py`
for the categories, `scripts/declines.py` for the Gruntz declines that the
native trace records. The test files are aseries, demidovich, formal,
gruntz, limits, lseries, nseries, order and series.

Raw outputs are in `raw/tl/` (tags `sympy`, `before` = main `5038ef02`,
`after` = this branch). The derived reports are in `results/`.

## Counts (identical / asserts)

| file | before | after | asserts |
|---|---|---|---|
| aseries | 5 | 5 | 9 |
| demidovich | 53 | 54 | 54 |
| formal | 59 | 59 | 85 |
| gruntz | 103 | 103 | 107 |
| limits | 533 | 538 | 554 |
| lseries | 5 | 5 | 5 |
| nseries | 165 | 172 | 186 |
| order | 172 | 205 | 210 |
| series | 44 | 70 | 83 |
| **total** | **1139** | **1211** | **1293** |

- 72 asserts were fixed and none broke (`results/compare_before_after.txt`).
- Gruntz declines went from 87 to 84 (`results/declines.txt`).
- The three removed declines are:
  - two "zero of a function in symbolic leading term" (`re` limits);
  - the `Eq(im(t)/...)` condition (`re((1 - exp(-2*I*pi*z))/(1 - exp(-2*I*pi*z/5)))` at 0).
- The 67 oscillating "series in mrv variable" declines are unchanged. All of them are taken over by the native bounded/AccumBounds walker of wave 41.

## Fixes

| fix | computation | files | evidence (fixed asserts) |
|---|---|---|---|
| Order of a product with noncommutative factors absorbs them (`Mul.flatten`) | Rust | `ops.rs` | order_noncommutative ×4 |
| Conjugate of an Order (`Order._eval_conjugate`) | Rust | `complex_parts/evaluate.rs` | order_conjugate_transpose ×4 |
| `Order.subs` with a new value carrying its own order: order at that order's point | Rust | `subs.rs` | order_subs_limits ×9 |
| Order construction: terms of the expanded Add with their orders, factors independent of x dropped (bound x of `Integral`), product of functions kept on a pole (PoleError branch) | Rust | `series/order.rs` | order_symbols ×2, simple_1, multivar_1 ×2, multivar_3, multivar_mul_1, issue_4279 ×2, issue_9192, order_at_infinity |
| `expand` renormalizes the orders that it rebuilds | Rust | `expand.rs` | order_subs_limits (sin/cos composed series) |
| Held Subs/Derivative/Integral: `_eval_nseries` and the `Expr.series` loop | Rust | `series/held.rs`, `series/mod.rs` | series_of_Subs ×7, issue_11313 ×6, issue_5223 ×5 |
| `Derivative._eval_as_leading_term` (first nonvanishing term of the operand's expansion) | Rust | `series/held.rs`, `series/mod.rs` | issue_11313 ×4, order_leadterm |
| Generic Taylor expansion of undefined functions with a Dummy `xi` (Subs of the derivatives) | Rust | `series/engine.rs` | issue_3978 ×3 |
| acoth expansions (log form at ±1 and on the cut; declined for an infinite argument) | Rust | `series/engine.rs` | hyperbolic, expbug4 ×2, expbug5, logbug4, seriesbug2, seriesbug2c |
| Leading term of a sum of monomials with real symbolic exponents | Rust | `series/mod.rs` (`power_sum_lead`) | demidovich leadterm |
| Limits through a positive dummy (`limitinf` replaces the variable by a positive Dummy); a value proved zero by the assumptions is 0 | Rust | `limits/mod.rs` | issue_9205 ×2, issue_13715, issue_21756 |
| `heuristics` of the reference, with AccumBounds arguments | Rust (`limits/heuristics.rs`); Python wrapper `series/limits.py` (argument passing only) | `limits/heuristics.rs`, `symox-py/src/series.rs` | heuristic |
| `Limit` with a non-symbol variable: construction and `free_symbols` | Rust | `limits/mod.rs`, `traverse.rs` | (Limit.free_symbols ×2, then counted under order/limits) |
| `Order.as_expr_variables` | Python (API envelope: tuple merge of the variables, no computation) | `python/symox/__init__.py` | as_expr_variables ×4 |
| `Expr._eval_as_leading_term`: an Order is its own leading term, otherwise the native leading term | Python (dispatch) | `python/symox/__init__.py` | order_leadterm |
| `lseries` of an expression free of x yields it once | Python (admission) | `python/symox/__init__.py` | issue_5223:2 |
| `as_leading_term`/`leadterm` admission tests free occurrences (`has_free`), not bound ones | Python (admission) | `python/symox/_expr_leading.py` | issue_11313:7 |

New tests are in `tests/test_native_limits2.py`.

## Remaining differences (82)

- **formal (26).** The `fps` domain is Python and has not been ported. These are recorded as gaps.
- **Form only, same value (30).**
  - Canonical forms of rational, trig or sign expressions: aseries_trig, issues 14393, 16708, 17751, 21785, 5740 and 6364, rational ×2, series_AccumBounds 4 and 5.
  - Order-term or expansion forms: abs:6, bug5 ×2, exp2, 3503, 3639, 4441:2, 12791/19534 (floats), 21245, 24266:3, 9173, 6318, aseries simple:5.
  - `log(x*y)` against `log(x) + log(y)`: ln_args:3.
  - Set order: free_symbols 4 and 5.
  - simple_8:5.
  - Dummy `_xi` printed against `_xi_1` (the values are equal): issue_3978:2.
- **Gaps, declined by symox (16).**
  - Symbolic or log exponents: aseries simple 2, 6 and 7; bug2, issue_3463, log_power1, 14384:3, 6350.
  - A log power in the O term: 3508, 5223:19.
  - exp_sqrt_1.
  - 21938 (weaker but true).
  - series1_failing ×2 (XFAIL in SymPy).
  - 3978:5: no generic Taylor expansion for a host `Function` subclass.
  - 5223:12: `sin(x)**y` with `logx`.
- **Deviations to validate (9).**
  - issue_5172:2, in both gruntz and limits: SymPy returns 0, while symox gives the value that SymPy's own test expects.
  - 19770:2: 0 where SymPy holds the limit.
  - 7535 ×3: `AccumBounds(-oo, oo)` where SymPy holds the limit.
  - gruntz_eval_special_slow:2: SymPy times out.
  - gruntz_eval_special_fail: PoleError where SymPy raises NotImplementedError.
  - order_failing_due_to_solveset: SymPy raises TypeError; symox returns `O(exp(-3*x**2))`.
- **Harness artifact.** test_sympyissue_22986: symox computes `sqrt(2)`.

## Reproduce

```
export WORK=<scratch>/work
bash scripts/run_tl.sh after symox <worktree>
python3 scripts/compare.py before after
python3 scripts/classify.py after -v
python3 scripts/declines.py after
```

## Validation

- The full suite is SUITE CLEAN: 60 503 passed, 6 xfail, 45 skipped, 60 554 collected, 694 files. The evidence is in `~/symox-ci/evidence/limits2-20261106`.
- `cargo test --release`: 871 + 1 + 0 + 5.
- `cargo clippy --release -- -A clippy::all -D clippy::correctness` is clean.
- `tests/test_bound_calculus_compat.py` expected NotImplementedError from `Limit(f(x) + y, x, 0).subs(x, 2)`. SymPy 1.14 returns `Limit(y + f(2), 2, 0)`, and symox now does too; the test was updated to match.
