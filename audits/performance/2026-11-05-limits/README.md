# limits: native Gruntz API, floats, symbolic exponents, reference errors (2026-11-05)

Branch `feat-limits-20261105`, base main 15e75d3a. Reference SymPy 1.14.0
(`~/symox-ci/venv`). The parity notes are in `docs/algorithmic-parity.md`,
section "limits (2026-11-05)".

"before" is the base tree (main 15e75d3a, with its decline trace
`_core._limit_trace`). "after" is this branch.

## Reproduction

```
cd audits/performance/2026-11-05-limits
export WORK=<scratch dir>
bash scripts/run_tl.sh after symox <worktree>   # statement harness, 9 files
bash scripts/run_pt.sh after <worktree>         # pytest on symox copies
python scripts/compare.py before after          # identical asserts, changes
python scripts/classify.py after                # differing asserts by cause
python scripts/declines.py after                # native Gruntz declines
PYTHONPATH=<worktree>/python python scripts/probe.py symox raw/probe/after.json
python scripts/probe.py sympy raw/probe/sympy.json  # symbolic-exponent probe
```

`raw/tl/sympy_*.json` are the SymPy 1.14 values of every evaluated
assertion side; `raw/tl/{before,after}_*` the symox values and decline
traces; `raw/pt/` the pytest runs of the copies of the SymPy test files.

## Fixes and where the computation runs

| Fix | Computation |
| --- | --- |
| `compare`, `sign`, `mrv`, `rewrite`, `mrv_leadterm`, `limitinf`, `gruntz` of `sympy.series.gruntz` | Rust `limits/gruntz/api.rs` (Python: `SubsSet` container and argument conversion only) |
| `Limit.doit` replaces floats by rationals (`nsimplify`) | Rust `limits/floats.rs`, called from `functions.rs` and `_limit_nsimplify` |
| `(x - x0)**p`, symbolic `p`: the `Limit.doit` leading-term step (`0`, `coeff`, `oo*sign(c)`, `-oo*sign(c)*(-1)**(1 + ex)`, `Not sure of sign of ex`) | Rust `limits/symbolic_power.rs` |
| `b**p` at `±oo`, symbolic `p`: `(±oo)**p`, held product, `Not sure of sign of -n` | Rust `limits/symbolic_power.rs` |
| A held or error outcome of these two steps stands only when the native engine does not decide | Rust `limits/mod.rs` (`limit_one_sided_raw`) |
| `Limit.pow_heuristics` for a constant base, `exp(+-oo*t)` as `exp.eval` (`(I/2)**x` is 0, `(-2)**x` is `zoo`) | Rust `limits/symbolic_power.rs` (`constant_base_power`, `exp_of_limit`) |
| `powsimp` of the expression before `pow_heuristics` and a second engine pass (`exp(b*x)*exp(-x)` raises the reference sign error) | Rust `limits/mod.rs` (`engine_or_pow_heuristics`) |
| `Result depends on the sign of ...` from the rewrite step, variable positive in the message | Rust `limits/gruntz.rs` |
| Sign of an expanded leading coefficient read after `factor_terms` (`limitinf((b*x - x)/log(x))` is `oo*sign(b - 1)`) | Rust `limits/gruntz.rs` (`infinity_with_sign_of`) |
| `limitinf`: a coefficient `I*a` gives `c0*oo` | Rust `limits/gruntz.rs` (`infinity_with_sign_of`) |
| `exp(I*pi*k)` coefficients written `(-1)**k` in a symbolic sign (issue 8634) | Rust `limits/gruntz.rs`, `limits/mod.rs` (`unit_powers`) |
| Bessel J/Y/I/K at `±I*oo` | Rust `limits/mod.rs` (`bessel_imaginary_value`), `limits/gruntz.rs` |
| `sinh`/`cosh`/`tanh` of an argument with a finite real limit kept (issue 26991) | Rust `limits/gruntz.rs` (`prepare_node`) |
| `zoo` limit point: `NotImplementedError` | Rust `limits/mod.rs` |
| The reference `NotImplementedError` raised instead of a held `Limit` | Rust message test `reference_not_implemented`; PyO3 `_limit(..., raise_nie)` raises it |
| `Limit(..., dir=0)`: `TypeError` | Python (argument admission) |
| `Limit.doit` of a held limit: oscillating route, then the native reference `NotImplementedError` | Python (routing) |

No Python algorithm was added or extended.

## Results

Statement harness (assertion sides equal to SymPy 1.14, `compare.py`):

| file | total | before | after |
| --- | --- | --- | --- |
| limits | 554 | 520 | 533 |
| gruntz | 107 | 95 | 103 |
| series | 83 | 44 | 44 |
| nseries | 186 | 165 | 165 |
| order | 210 | 172 | 172 |
| formal | 85 | 59 | 59 |
| demidovich | 54 | 53 | 53 |
| aseries | 9 | 5 | 5 |
| lseries | 5 | 5 | 5 |
| total | 1293 | 1118 | 1139 |

21 sides fixed, none broken.

pytest on the symox copies of the SymPy files:

| file | before | after |
| --- | --- | --- |
| test_limits | 32 failed, 141 passed | 15 failed, 158 passed |
| test_gruntz | 19 failed, 19 passed | 1 failed, 37 passed |
| test_series | 14 failed, 26 passed | 14 failed, 26 passed |
| test_nseries | 17 failed, 52 passed | 17 failed, 52 passed |
| test_order | 20 failed, 37 passed | 20 failed, 37 passed |

Native Gruntz declines (trace, statement harness): 90 before, 88 after. Of
the 88, 10 now raise the reference `NotImplementedError` (`Result depends
on the sign of ...` 7, `Limits at complex infinity are not implemented` 2,
`Not sure of sign of 3 - w` 1), where the base returned a held `Limit`. 68
are `series in mrv variable`: oscillating `sin`/`cos`/`tan`/`frac` where
the reference gruntz also fails and its `heuristics` answer; the Python
`AccumBounds` walker answers most of them (see the remaining gaps).

Symbolic-exponent probe (`scripts/probe.py`, 54 limits with symbolic
exponents or constant bases, `raw/probe/`): 15 equal to SymPy 1.14 before,
48 after. The 6 others stay held (`I**x`, `sin(x)**w` at `0+`, `x**n` at
`-oo`, `x**a/(x**a + 1)`, `(3/x)**w` at `0+`, `Mod(x, a)` at `a`).

## Validation

- `run_suite.py` (sharded suite): SUITE CLEAN, 60 227 passed, 45 skipped,
  6 xfailed, 60 278 collected and finished (main: 60 211 / 60 262; the 16
  new tests are this lot's), no failure.
- `cargo test --release --lib --bins --tests`: 860 passed, 1 ignored; 1
  passed; 0 failed. Doctests (`--test-threads=1`): 5 passed.
- `cargo clippy --release -- -A clippy::all -D clippy::correctness`: exit 0.
- `tests/test_native_limits_20261105.py` (16 tests) passes unchanged on
  SymPy 1.14 (imports switched to `sympy`).

## Remaining gaps

- 68 oscillating declines are answered by the existing Python walker
  (`calculus/_accum_limit.py`); a Rust port needs `AccumBounds` arithmetic
  in Rust.
- `heuristics()` is not exposed (it would be a Python algorithm).
- Value forms of the meromorphic leading-term path (`test_rational`,
  issues 5740, 6364, 14393, 16708, 17751), `series_AccumBounds`, 21785.
- `Limit` with a non-symbol variable (9205), `re`/`im` limits (21756).
- `Mod(x, a)` at `a`: held, SymPy raises `Not sure of sign of 1/a` (25230).
- `x**n` (integer `n`) at `-oo` and `2*(x + 1)**w` at `-oo` stay held;
  SymPy gives `(-oo)**n` and `2*(-oo)**w`. Also held: `I**x` at `oo`
  (SymPy `nan`), `(3/x)**w` at `0+` (SymPy `oo**w`), `sin(x)**w` at `0+`
  and `x**a/(x**a + 1)` at `oo` (SymPy raises).
- The API `mrv_leadterm` returns the expanded coefficient (`b/x - 1/x`
  where SymPy gives `(b - 1)/x`); equal values, different form.
- `test_aseries_trig` of test_gruntz (a `log(pi/2)` cancel gap).
- Order, nseries, fps and undefined-function series classes unchanged.

## Deviations for validation

- `gruntz()` falls back to the native limit (with its heuristic steps)
  after the pure Gruntz algorithm fails, and raises `PoleError` only if
  that is undecided too; SymPy's `gruntz` raises directly;
- a float with no small rational and more than 12 significant digits, a
  magnitude outside `[1e-12, 1e15)` or an arbitrary precision float is
  kept as a float in `Limit.doit` (the expression is then computed with
  floats, as before), where SymPy's `identify` may find a closed form.

Unchanged by this lot (they differed before): issue 19770 (`cos(m*x)/x`
at `oo` gives 0, SymPy holds), issue 7535 (`AccumBounds(-oo, oo)`, SymPy
holds), issue 13715 (`p` where SymPy gives 0).
