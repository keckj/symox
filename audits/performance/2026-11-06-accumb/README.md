# accumb: native AccumulationBounds arithmetic and oscillating-limit walker (2026-11-06)

Branch `feat-accumb-20261106`, base main 5ebd60c0. Reference SymPy 1.14.0
(`~/symox-ci/venv`). The parity notes are in `docs/algorithmic-parity.md`,
section "accumb (wave 41)".

"before" is the base tree (main 5ebd60c0: Python `AccumulationBounds`
methods, Python walker `_accum_limit._walk`, Python function evaluation
`_accum_eval._sin` ...). "after" is this branch.

## Reproduction

```
cd audits/performance/2026-11-06-accumb
export WORK=<scratch dir>
bash scripts/run_tl.sh after symox <worktree>   # statement harness, 9 files
bash scripts/run_pt.sh after <worktree>         # pytest on symox copies
python scripts/compare.py before after          # identical asserts, changes
python scripts/classify.py after                # differing asserts by cause
python scripts/declines.py after                # native Gruntz declines
PYTHONPATH=<tree>/python python scripts/bench.py raw/bench/<tag>.json
```

The scripts are those of `2026-11-05-limits`. `tl.py` now counts, besides the
statement values, how each limit the native `_limit` declines ends in
`oscillating_limit`: the native bounded decomposition
(`_core._limit_bounded_affine`), the walker (Python before, native after),
or held (`raw/tl/*.stats.json`). `bench.py` is new.

`raw/tl/sympy_*.json` are the SymPy 1.14 values of every evaluated
assertion side (copied from `2026-11-05-limits`); `raw/tl/{before,after}_*`
the symox values, decline traces and walker counts; `raw/pt/after_*` the
pytest runs of the copies of the SymPy test files; `raw/bench/` the timings.

## Fixes and where the computation runs

| Fix | Computation |
| --- | --- |
| `AccumulationBounds.__new__` (real/order checks, `AB(a, a)` is `a`, the node) | Rust `limits/accum.rs` (`ab_new`); Python: `sympify` of the arguments |
| `__add__`, `__sub__`, `__rsub__`, `__mul__`, `__truediv__`, `__rtruediv__`, `__pow__`, `__rpow__`, `__neg__`, `__abs__` | Rust (`ab_add` ... `ab_rpow`, `ab_abs`, dispatched by `op_add` ... `op_pow` as Python operators dispatch); Python: `_sympifyit` and the `isinstance(other, Expr)` admission |
| `_eval_is_ge` (the three dispatches), `Lt`/`Gt`/... of bounds | Rust (`is_ge`, `rel_value`, `decide`) |
| `intersection`, `union`, `__contains__` | Rust (`ab_intersection`, `ab_members`, `ab_union`, `ab_contains`); Python: the reference `TypeError` admission, `FiniteSet`/`EmptySet` and `S.true`/`S.false` wrapping of the native result |
| `Add.flatten`/`Mul.flatten` folding of Numbers and bounds into the coefficient | Rust (`fold`); Python: the scan that decides whether an argument holds bounds |
| Function values at bounds and at real infinities: `sin`, `cos`, `tan`, `cot`, `sec`, `csc`, `exp`, `log`, `frac`, hyperbolic functions at `+-I*oo`, `atan`, `atanh`, `asech` at `zoo` (the `AccumBounds` branches of the reference `eval`s) | Rust (`evaluate`, `eval_sin` ... `eval_hyperbolic`); Python: `_accum_eval.evaluate` routes by head name |
| `log(AccumBounds(a, b))` with `a` negative or of unknown sign is `nan` (reference `log.eval`; before: left unevaluated) | Rust (`eval_log`) |
| The oscillating-limit walker (`_accum_limit._walk`, `_combine`, `_bounded_value`, `_held_oscillation`, `oo*sign(s)` arguments, the two-sided case, the fallback to the bounded decomposition) | Rust (`oscillating_limit`, `walk`, `bounded_affine_value`) |
| `__contains__` returns `S.true`/`S.false` for a finite element (the reference `And` value) and a Python bool for `+-oo` | Python (return-type wrapping of the native boolean) |

No Python algorithm was added or extended; 873 lines of Python bounds
arithmetic, function evaluation and walker were removed.

## Results

Statement harness (assertion sides equal to SymPy 1.14, `compare.py`):

| file | total | before | after |
| --- | --- | --- | --- |
| limits | 554 | 533 | 533 |
| gruntz | 107 | 103 | 103 |
| series | 83 | 44 | 44 |
| nseries | 186 | 165 | 165 |
| order | 210 | 172 | 172 |
| formal | 85 | 59 | 59 |
| aseries | 9 | 5 | 5 |
| lseries | 5 | 5 | 5 |
| demidovich | 54 | 53 | 53 |
| total | 1293 | 1139 | 1139 |

No assertion side changed (154 still differ, see `classify.py after`); no
series, order or nseries regression. Six of the 154 involve bounds, all
unchanged: issue 7535 (three asserts, `AccumBounds(-oo, oo)` where SymPy
holds), `test_series_AccumBounds` (two asserts, the same interval in another
form: symox `simplify` returns another form of the sum, which the native
bounded decomposition then takes) and `test_heuristic` (`heuristics()` is
not exposed). The pytest
runs of the SymPy test copies fail the same tests as `2026-11-05-limits`
(limits 15, gruntz 1, nseries 17, order 20, series 14 failures).

Gruntz declines (`declines.py`): 88 before (68 "series in mrv variable",
oscillating), 87 after (67). The missing one is `limit(sin(log(cos(x))),
x, oo)`: `log(cos(oo))` is now `nan` as in SymPy, so the walk stops before
the native Gruntz call that recorded the decline; the value is the held
`Limit` before, after and in SymPy (`test_limits.py:881`). The decline
itself is the reference Gruntz `NotImplementedError`; what changed is where
the declined limits go next:

| route of the declined limits (`raw/tl/*.stats.json`) | before | after |
| --- | --- | --- |
| native bounded decomposition | 19 | 19 |
| walker | 38 (Python) | 38 (native) |
| held | 25 | 25 |

So every oscillating limit of the SymPy test files is now computed in
Rust. The held ones match SymPy when it holds or raises, except the
recorded deviations below.

Replay of the 290 distinct `oscillating_limit` calls recorded while running
the repository tests and the SymPy test files with the Python walker (base):
288 identical (value and exception class), 2 not replayable (unparsable
`srepr`), 0 different. Two probes of 188 and 158 bounds operations and
limits (arithmetic, powers, comparisons, sets, function values, limits)
give identical output before and after, except `log` of bounds reaching
below zero, which now matches SymPy.

Timings (`bench.py`, best of 5, a fresh symbol each run). Two interleaved runs of each tree on a loaded machine (other builds
running); the best of the two runs, in ms (`raw/bench/`):

| case | before | after | ratio |
| --- | --- | --- | --- |
| `limit(sin(x), x, oo)` | 0.77 | 0.46 | 1.7x |
| `limit(sin(1/x), x, 0)` | 0.81 | 0.47 | 1.7x |
| `limit(cos(1/x), x, 0, '+-')` | 1.05 | 0.28 | 3.7x |
| `limit(exp(sin(x)), x, oo)` | 1.37 | 0.80 | 1.7x |
| `limit(sin(x)**2, x, oo)` | 0.88 | 0.54 | 1.6x |
| `limit(1/(2 + sin(x)), x, oo)` | 1.50 | 0.64 | 2.3x |
| `limit(sin(x) + cos(x), x, oo)` | 1.02 | 0.48 | 2.1x |
| `limit(frac(x), x, oo)` | 0.78 | 0.42 | 1.9x |
| `limit(tan(x), x, oo)` | 0.71 | 0.38 | 1.9x |
| `limit(log(2 + sin(x)), x, oo)` | 1.65 | 0.78 | 2.1x |
| `limit(x*sin(x), x, oo)` | 1.08 | 0.51 | 2.1x |
| `limit(sin(x)**x, x, oo)` (held) | 1.74 | 0.96 | 1.8x |
| bounds arithmetic, 29 x 14 operations | 94.18 | 40.63 | 2.3x |

## Deviations from SymPy (all present before this lot; see the parity notes)

- `limit(x*sin(x), x, oo)`, `limit(exp(x)*sin(x), x, oo)`,
  `limit(x + x*sin(x), x, oo)`, `limit(sin(1/x)/x, x, 0)`,
  `limit(x*log(2 + sin(x)), x, oo)`: `AccumBounds(-oo, oo)` (or
  `AccumBounds(0, oo)`); SymPy gives `oo*sign(AccumBounds(...))`.
  Recorded intentional divergence.
- `limit(sin(x)*log(x), x, oo)` is `AccumBounds(-oo, oo)` and
  `limit(Max(sin(x), 0), x, oo)`, `limit(x**sin(x), x, oo)` are held; SymPy
  raises an internal `TypeError` (`_eval_nseries() missing ... 'cdir'`).
- `AccumBounds(1, 2)*nan` prints `nan*AccumBounds(1, 2)` (argument order of
  `Mul`); `latex(AccumBounds(1, 2))` is the generic function form, not
  `\left\langle 1, 2\right\rangle`.
