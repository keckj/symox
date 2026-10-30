# odesolve (wave 35): native characteristic-equation solutions

Base d78f7e5c, branch feat-odesolve-20261030j. Reference SymPy 1.14.0.

## What moved to Rust

`crates/symox-core/src/ode/constcoeff.rs` (binding
`crates/symox-py/src/ode_solve.rs`, admission
`python/symox/solvers/ode/_odesolve_native.py`, hooks at the start of the
reference helpers in `python/symox/solvers/ode/nonhomogeneous.py`):

| Reference helper | Native part | Used by |
| --- | --- | --- |
| `_get_const_characteristic_eq_sols` | characteristic polynomial, `roots(multiple=True)` order (`ordered`), `rootof(chareq, k)` for irreducible polynomials without radical roots, `exp`/`sin`/`cos` fundamental solutions, `collectterms` | constant-coefficient homogeneous, undetermined coefficients, variation of parameters |
| `_get_simplified_sol` | `collect` on the sorted `collectterms` patterns (`exact=False`), `powsimp` | the same three hints |
| `_get_euler_characteristic_eq_sols` | indicial polynomial, explicit `rootof` roots, general solution with the caller's numbered constants, fundamental solutions | Euler homogeneous, Euler undetermined coefficients, Euler variation of parameters |

`poly/rootof.rs` gains `all_roots_irreducible_in`, the existing
irreducible `RootOf` enumeration with the defining polynomial written in a
given generator (the reference `PurePoly` keeps the dummy).

Python keeps the admission checks, `get_numbered_constants` (so `C1..Cn`
are numbered by the reference code), and the particular solutions (trial
set and linear solve of undetermined coefficients, Wronskian integrals of
variation of parameters). Any decline (Float, `I`, `sqrt(2)`, `pi` or
`x`-dependent coefficients, reducible polynomials with a `RootOf` factor,
Euler roots needing `RootOf`, non-integer keys) runs the unchanged
reference body.

## Differential test (gen_eqs.py, runner.py, compare.py)

338 equations (160 random constant-coefficient, 28 symbolic, 80 random
Euler, 10 symbolic Euler, 60 undetermined coefficients, 30 variation of
parameters); fields: hinted solution, `simplify=False` solution, default
`dsolve`. `stock` = this build with the three kernels disabled.

| Family | field | after == SymPy | before == SymPy | after == before |
| --- | --- | ---: | ---: | ---: |
| const | hint | 152/160 | 152/160 | 160/160 |
| const | nosimp | 158/160 | 118/160 | 118/160 |
| const | default | 152/160 | 152/160 | 160/160 |
| csym | hint | 24/28 | 24/28 | 28/28 |
| csym | nosimp | 24/28 | 21/28 | 21/28 |
| csym | default | 23/28 | 23/28 | 28/28 |
| euler | all three | 240/240 | 240/240 | 240/240 |
| esym | hint/nosimp/default | 9/9/8 of 10 | 9/9/8 of 10 | 30/30 |
| undet | hint | 55/60 | 55/60 | 60/60 |
| undet | nosimp | 59/60 | 45/60 | 45/60 |
| undet | default | 54/60 | 54/60 | 60/60 |
| total (non-vop) | | 967/1014 | 905/1014 | 950/1014 |
| vop | hint/nosimp/default | 25/19/26 of 29 | 25/18/26 of 29 | 82/87 |

before == stock on 1014/1014 and 87/87 fields. Every after/before
difference is a `simplify=False` solution whose fundamental solutions are
now listed in the reference order (the stock `roots` sorted by
`default_sort_key` only); each one now equals SymPy. 0 regressions.
`vop_016` is skipped (`ODE_SKIP`): symox exceeds 6 GB inside the native
integration before and after; SymPy times out at 30 s.

Remaining differences, all identical before: cubic/quartic radical forms,
the `collect` partial exponent match (`const_025`, `const_036`,
`undet_007`), the `re`/`im` sign order of `esym_00`, the `RecursionError`
message of `esym_09`, `undet_014`/`undet_017` where SymPy raises
`NotImplementedError` (its `re`/`im` of casus irreducibilis radicals) and
symox solves, and the vop integration forms (`Ei` with `exp_polar`,
`Integral` forms, SymPy timeouts).

## SymPy test_single examples (examples.py)

The six ported hints, `hint` and `simplify=False` (124 examples,
`var_of_parameters_15` skipped: symox exceeds 6 GB before and after):

| Hints | after == SymPy | before == SymPy |
| --- | ---: | ---: |
| all (248 fields) | 223 | 216 |

0 regressions. Remaining: CRootOf printed with `x` instead of `_x` by
SymPy (a SymPy cache artifact of the module-level expected solutions),
`lin_hom_38` nosimp (SymPy repeats C4/C5), `lin_hom_41` (an `I`
coefficient, declined), `undet_38`/`undet_39` symbolic quadratic form, the
variation-of-parameters particular-solution forms, `euler_undet_03/06`,
`euler_var_05`.

## Validation

- tests/test_native_odesolve.py: 54 passed (29 SymPy strings, default hint,
  native helpers equal to the reference bodies on 15 polynomials, declines).
- Caller tests (tests/test_native_ode*.py, tests/test_ode_*.py,
  tests/test_linear_ode_systems.py): 738 passed, 2 skipped.
- Copied SymPy test_single.py and test_ode.py: 22 failed, 51 passed, the
  same 22 failures as the base.
- Full suite on demeter, 24 shards plus the local compiler files
  (raw_suite_shards.txt): 58598 passed, 39 skipped, 4 xfailed, 5 failed in
  shard 6 (`test_polymatrix_stats_protocol.py`, 5 tests). They pass alone
  and in the shard locally; the shard rerun on demeter fails
  `test_native_ratint_logpart.py::test_pending_interrupt_is_never_a_refusal`
  (an asynchronous `_thread.interrupt_main` whose late KeyboardInterrupt
  lands in later tests) 3 times out of 3 on the base build too.
- cargo test -p symox --release: 806 passed, 1 ignored, 1 + 5 doctests
  (doctests rerun with one thread after a "Couldn't compile").
- clippy correctness: clean.
- ODE corpus: {'cls': 303, 'hint': 279, 'def': 274, 'n': 315}, baseline,
  0 per-case differences in classify/hinted/default.
- Integration corpus: {"COND": 20, "DIFF": 10, "EQUIV": 8, "MISSING": 1,
  "SAME": 1171}, baseline.

## Timings (bench.py, ms, median of 3 fresh processes; repeat = median of 5)

| Case/op | SymPy first | before first | after first | SymPy repeat | before repeat | after repeat |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| lin_const_coeff_hom_01/hint | 47.6 | 13.9 | 12.4 | 30.1 | 12.3 | 10.9 |
| lin_const_coeff_hom_01/dsolve | 157.5 | 30.9 | 28.8 | 54.5 | 29.0 | 26.9 |
| lin_const_coeff_hom_10/hint | 103.4 | 38.3 | 35.3 | 55.2 | 35.9 | 32.9 |
| lin_const_coeff_hom_10/dsolve | 111.5 | 39.4 | 36.6 | 61.5 | 37.2 | 34.4 |
| lin_const_coeff_hom_20/hint | 86.8 | 32.6 | 30.6 | 49.9 | 30.4 | 28.3 |
| lin_const_coeff_hom_20/dsolve | 265.3 | 73.3 | 71.0 | 128.8 | 71.2 | 67.9 |
| lin_const_coeff_hom_30/hint | 73.5 | 31.3 | 29.1 | 47.0 | 29.7 | 27.1 |
| lin_const_coeff_hom_30/dsolve | 174.2 | 35.6 | 34.3 | 60.9 | 33.5 | 31.4 |
| lin_const_coeff_hom_40/hint | 39.0 | 23.1 | 22.8 | 23.5 | 21.6 | 20.8 |
| lin_const_coeff_hom_40/dsolve | 132.4 | 26.3 | 25.7 | 30.7 | 24.7 | 23.5 |
| lin_const_coeff_hom_47/hint | 62.4 | 27.0 | 26.1 | 32.8 | 24.4 | 23.1 |
| lin_const_coeff_hom_47/dsolve | 186.8 | 56.7 | 53.1 | 68.5 | 53.0 | 49.4 |
| euler_hom_01/hint | 35.9 | 13.7 | 12.4 | 21.9 | 11.7 | 10.8 |
| euler_hom_01/dsolve | 145.0 | 29.6 | 27.6 | 47.0 | 27.3 | 25.6 |
| euler_hom_05/hint | 63.7 | 48.3 | 28.2 | 31.8 | 30.0 | 26.3 |
| euler_hom_05/dsolve | 265.5 | 81.0 | 60.0 | 142.6 | 61.8 | 57.3 |
| euler_hom_09/hint | 72.4 | 35.3 | 15.2 | 32.1 | 16.7 | 13.3 |
| euler_hom_09/dsolve | 80.8 | 37.7 | 16.6 | 39.0 | 17.9 | 14.5 |
| undet_01/hint | 203.5 | 61.8 | 59.0 | 107.0 | 58.6 | 56.3 |
| undet_01/dsolve | 353.3 | 68.3 | 66.0 | 192.8 | 65.6 | 62.8 |
| undet_10/hint | 145.7 | 29.8 | 28.1 | 85.1 | 27.5 | 26.2 |
| undet_10/dsolve | 296.7 | 66.5 | 63.1 | 145.9 | 63.8 | 60.3 |
| var_of_parameters_01/hint | 588.9 | 68.3 | 63.3 | 358.5 | 64.8 | 60.4 |
| var_of_parameters_01/dsolve | 342.5 | 67.9 | 64.8 | 181.7 | 65.7 | 61.9 |
| euler_undet_01/hint | 205.5 | 27.5 | 26.1 | 72.8 | 25.6 | 23.9 |
| euler_undet_01/dsolve | 237.5 | 53.1 | 53.7 | 105.3 | 50.5 | 50.4 |
| euler_var_01/hint | 416.9 | 26.4 | 26.5 | 243.6 | 24.3 | 24.6 |
| euler_var_01/dsolve | 461.2 | 56.6 | 55.1 | 330.5 | 53.7 | 52.0 |

Geometric means: first call SymPy 146.8, before 38.7, after 33.9 ms
(after/before 0.88); repeated call SymPy 72.0, before 33.7, after 31.3 ms
(after/before 0.93). The homogeneous part was already a small share of
`dsolve` (classification and the particular solution dominate); the
largest gains are the Euler cases (euler_hom_09 first call 35.3 to 15.2 ms).

## Files

- gen_eqs.py, runner.py (`ODE_SKIP`), compare.py: differential test.
- examples.py: SymPy test_single examples.
- bench.py: timings.
- raw_ref_sympy.jsonl, raw_before_symox.jsonl, raw_after_symox.jsonl,
  raw_stock_symox.jsonl and the `_vop` variants: differential outputs.
- raw_ex_sympy.jsonl, raw_ex_before.jsonl, raw_ex_after.jsonl: examples.
- raw_bench_{sympy,before,after}.json: timing samples.
- raw_suite_shards.txt, raw_callers.txt: test summaries.
