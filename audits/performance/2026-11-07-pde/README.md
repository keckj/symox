# pde corpus, 2026-11-07 (`feat-pde-20261107`, wave 45)

The lot ports SymPy 1.14's `sympy/solvers/pde.py` to Rust
(`crates/symox-core/src/ode/pde.rs`, binding `crates/symox-py/src/pde.rs`):

- `classify_pde` (tuple and `dict=True`, power reduction, `Eq` with a
  right-hand side, the reference exceptions);
- `pdsolve` with the `_desolve` logic (`default`, named hints, `all`,
  `all_Integral`, `best`, the reference `KeyError`/`ValueError`/
  `NotImplementedError`);
- the four solvers `1st_linear_constant_coeff_homogeneous`,
  `1st_linear_constant_coeff`, `1st_linear_constant_coeff_Integral`
  (`_handle_Integral`) and `1st_linear_variable_coeff`;
- `checkpdesol` (sequences, `solve_for_func`, the final
  `NotImplementedError`);
- `pde_separate`, `pde_separate_add`, `pde_separate_mul` and `_separate`.

The `lie_group` ODE solver (wave 44) called `pdsolve` through a host
callback to Python; it now calls the native `pdsolve` directly
(`ode/lie_group.rs`, `native_pdsolve_rhs`). The host callbacks are
`simplify`, `expand`, `integrate`, `solve`, `doit` and `dsolve`.

## Corpus

`scripts/corpus.py` evaluates 791 cases with live SymPy 1.14 and with
symox, on the equations of SymPy's `test_pde.py` and 30 in-house equations
(symbolic coefficients, power reduction, `Eq` with a right-hand side,
isotropic complex characteristics, `xi` in the equation, an arbitrary
function `F` in the equation, nonlinear and second-order equations):

| Kind | Cases |
|---|---:|
| `classify_pde` tuple / dict | 53 / 53 |
| `pdsolve` (`default`, each hint, `all`, `all_Integral`, `best`, wrong hints) | 477 |
| `pdsolve` with `func` | 53 |
| `checkpdesol` of `pdsolve` results / of 13 hand-written candidates (with and without `func`) | 54 / 26 |
| `pde_<hint>` called directly | 36 |
| `pde_separate`, `_add`, `_mul` | 23 |
| error paths of `_desolve`, `checkpdesol`, `pde_separate` | 16 |

Reproduce:

```sh
python scripts/corpus.py sympy sympy.json                  # live SymPy 1.14
PYTHONPATH=<worktree>/python python scripts/corpus.py symox symox.json
python scripts/corpus.py report sympy.json symox.json report.json
```

The classes are those of the 2026-11-06 liegroup audit: SAME (equal
`srepr` with dict order ignored, or the same error class), DECLINED, WRONG
(SymPy's `checkpdesol` rejects a symox `pdsolve` answer), DIFF
(`validated_by_sympy` is true when SymPy's `checkpdesol` accepts the
differing `pdsolve` answer; false means it was not checked or did not
conclude).

## Results

"Before" is main `f7888825` (the Python body of `pde.py`). "After" is this
branch. Files: `sympy.json.gz`, `symox-before.json.gz`,
`symox-after.json.gz`, `report-before.json`, `report-after.json`.

| | SAME | DIFF | DECLINED | WRONG |
|---|---:|---:|---:|---:|
| before | 714 | 77 | 0 | 0 |
| after | **716** | 75 | 0 | 0 |

| Kind | SAME after (before) | DIFF after (before) |
|---|---:|---:|
| classify | 106 (106) | 0 |
| pdsolve | 282 (282) | 36 (36) |
| all / all_Integral / best | 180 (180) | 32 (32) |
| check | 74 (74) | 5 (5) |
| method | 34 (34) | 2 (2) |
| extra (error paths) | 17 (15) | 0 (2) |
| separate | 23 (23) | 0 |

- 714 stayed SAME, 2 went from DIFF to SAME, none got worse.
- The two fixed cases are `err:check_nosolve` and `err:check_sol4`. The
  Python body of `checkpdesol` had its final
  `raise NotImplementedError('Unable to test if ...')` as dead code after
  the `return` of `_combine_integrals`, so it returned `None`. The native
  port raises as SymPy does.
- 759 of the 791 symox results are identical before and after. The other
  30 are `all`/`best` dicts whose key order changed: the native port keeps
  the `classify_pde` order. SymPy iterates a set of hint names, so its
  order follows the string hash seed; the report ignores dict order.
- The 75 remaining DIFFs give the same symox values as before; none comes
  from the port (below).

Native calls during the after run (no call declined to the Python body):

| Call | Native |
|---|---:|
| `pde_classify` | 144 |
| `pde_pdsolve` | 372 |
| `pde_checkpdesol` | 77 |
| `pde_solve_hint` | 12 |
| `pde_separate` | 18 |

Time (sum over the 791 cases): SymPy 573.7 s, symox before 208.5 s,
after **68.2 s**.

### Remaining differences (after)

| Cases | Cause |
|---|---|
| `k1`, `k2`, `k4`, `all1`, `m3`, `m4`, `m17`, `m18`, `m28` (`pdsolve`, named hints, `all`/`best`) and `method:k1:*` | Form of `simplify`/`Mul` outside the lot: `exp(x/2 - y)` for `exp(x/2)*exp(-y)`, the `_Integral` integrand with merged exponentials. SymPy's `checkpdesol` accepts every single-solution answer (33 validated) |
| `k5` (`pdsolve`, `all`, `best`, `check`) | SymPy's `Integral.doit` raises TypeError on the one-sided limit `(xi, x + y)` of the held log integral; symox keeps it held and answers (validated) |
| `m15`, `m16` | Isotropic characteristics `u_x + I*u_y`: SymPy returns `False` (and its `checkpdesol` raises AttributeError) or raises NotImplementedError; symox keeps its `_native_pde` characteristic fallback (validated) |
| `m18` | `xi` appears in the equation: SymPy's integration variable captures it and its `checkpdesol` gives `(False, x - xi + y)`; symox renames the coordinates (`(True, 0)`) |
| `m19` | `G(x, y)` in the equation: SymPy integrates `G(x, y)/x` in `x` with `y` held although the characteristic variable is `x*y`, and its own answer fails its `checkpdesol`; symox integrates along the characteristic (`checkpdesol` gives `(True, 0)`; SymPy's check of it did not conclude) |
| `m23` | `1st_linear_variable_coeff` of `sin(x)*u_x + cos(y)*u_y`: another form of the `solve`/`dsolve` result (validated) |

These are the voluntary deviations of symox's `pde.py` kept by the port
(free name `F_` for the arbitrary function when `F` is in the equation,
`Dummy` characteristic coordinates when `xi`/`eta` are in the equation,
renamed captured integration variables, `_combine_integrals` in
`checkpdesol`, the `_native_pde` fallbacks), and forms of the shared
`simplify`, `integrate` and `doit`.

## lie_group corpus and `dsolve` recursion

The lie_group corpus of wave 44 (`../2026-11-06-liegroup`, 428 cases)
was rerun on this branch (outputs in `liegroup/`):

- **413 SAME, 13 DIFF, 2 DECLINED, 0 WRONG**, the same verdicts case by
  case as wave 44 (base 412-413 / 13-14 / 2 / 0). The symox values are
  identical except for one `Dummy` index (`dsolve:sim2`).
- Native calls: 108 `infinitesimals`, 15 `checkinfsol`, 204 splits
  (176 in wave 44: the native `pdsolve` of `lg01`/`lg03` reaches more
  recursion levels, below) and 19 solves, with no decline.
- The `dsolve` cases, the only ones that go through `pdsolve`, take
  25.3 s on the SAME cases (25.5 s in wave 44). The 400 SAME cases
  below the deadline take 218.6 s against 161.8 s in wave 44 and 772.6 s
  in SymPy; the increase is in `infinitesimals` cases that do not call
  `pdsolve` (`fsum`, `lin1`) and comes from the machine load during the
  run (other builds and corpora in parallel). Run alone, `infinitesimals`
  of `fsum` with `hint='all'` makes the same 1 622 774 Python calls on
  main (`integ44` build) and on this branch, in 25.8 s and 25.7 s.

`dsolve(hint='lie_group')` of lie_group_01 and lie_group_03 raises
RecursionError, as in SymPy (6 min 22 s and 8 min 10 s):

| Case | wave 44 corpus | this branch, corpus | this branch, alone |
|---|---:|---:|---:|
| lg01 | 79.1 s | 72.9 s | 90.7 s |
| lg03 | 50.3 s | 56.2 s | 61.4 s |

The recursion goes through `pdsolve` -> `1st_linear_variable_coeff` ->
`dsolve` -> `lie_group`. The `pde.py` frames are now native, so each
level uses fewer Python frames: the traceback has 831 frames over 81
levels, and the native split count rises from 61 to 75 per case. No Rust
stack overflow occurs.
