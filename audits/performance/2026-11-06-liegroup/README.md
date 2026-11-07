# lie_group corpus, 2026-11-06 (`feat-liegroup-20261106`, wave 44)

The lot ports SymPy 1.14's `lie_group` ODE solver to Rust
(`crates/symox-core/src/ode/lie_group.rs`, binding
`crates/symox-py/src/ode_lie.rs`). The ported code covers:

- `infinitesimals` and the nine `lie_heuristic_*` heuristics, in SymPy's order;
- `checkinfsol`;
- `_lie_group_remove`;
- the `LieGroup` single-solver steps: the `_matches` split and
  `_get_general_solution` with the canonical coordinates.

## Corpus

`scripts/corpus.py` evaluates 428 cases with live SymPy 1.14 and with symox:

- `infinitesimals` for every hint (`default`, `all` and the nine heuristics);
- `checkinfsol` of the `all` generators;
- `'lie_group' in classify_ode`;
- `dsolve(..., hint='lie_group')`.

The equations come from three sources:

- SymPy's `test_lie_group.py`;
- the 20 `lie_group_*` rows of `test_single.py`;
- `audits/performance/2026-11-05-ode`.

Each symox case runs in a forked child that is killed at the 120 s
deadline, because SIGALRM cannot interrupt native code.

The report classes every case:

| Class | Meaning |
|---|---|
| SAME | Equal `srepr`, ignoring dict order and with `Dummy` indices numbered by first appearance, or the same error class. |
| DECLINED | symox raises or times out where SymPy answers. |
| WRONG | symox answers and SymPy's `checkinfsol` or `checkodesol` rejects the answer. |
| DIFF | Any other difference. `validated_by_sympy` records whether SymPy accepted symox's differing answer. |

Reproduce:

```sh
python scripts/corpus.py sympy sympy.json                  # live SymPy 1.14
PYTHONPATH=<worktree>/python python scripts/corpus.py symox symox.json
python scripts/corpus.py report sympy.json symox.json report.json
```

## Results

"Before" is main `583315a8`, whose `lie_group` was symox's own symmetry
search in Python and not SymPy's algorithm. "After" is this branch.
Files: `sympy.json.gz`, `symox-before.json.gz`, `symox-after.json.gz`,
`report-before.json` and `report-after.json`.

| | SAME | DIFF | DECLINED | WRONG |
|---|---:|---:|---:|---:|
| before | 319 | 96 | 11 | 2 |
| after | **413** | 13 | 2 | **0** |

By kind, after (before in brackets):

| Kind | Cases | SAME | DIFF | DECLINED | WRONG |
|---|---:|---:|---:|---:|---:|
| `infinitesimals` | 341 | 329 (249) | 10 (80) | 2 (10) | 0 (2) |
| `checkinfsol` | 31 | 31 (16) | 0 (14) | 0 (1) | 0 |
| `classify_ode` | 31 | 31 (31) | 0 | 0 | 0 |
| `dsolve` | 25 | 22 (23) | 3 (2) | 0 | 0 |

How the verdicts moved from before to after:

- 314 stayed SAME; 91 went from DIFF to SAME, 6 from DECLINED to SAME and
  2 from WRONG to SAME.
- 5 went from SAME to DIFF, and all of them are explained below.
  `h1_eq7` `default` and `abaco1_simple` now carry the `integrate` log
  form. `lg11` `all` and `bivariate` give the hash-seed generator. The
  old Python search matched SymPy's printed form on these four by chance;
  the port follows SymPy's algorithm, and SymPy's `checkinfsol` accepts
  every new generator. `dsolve` of `lg03` now raises SymPy's
  RecursionError (see below).

Native calls during the after run:

| Call | Native | Declined |
|---|---:|---:|
| `ode_lie_infinitesimals` | 108 | 0 |
| `ode_lie_checkinfsol` | 15 | 0 |
| `ode_lie_group_split` | 176 | 0 |
| `ode_lie_group_solve` | 19 | 0 |

No call declined to the Python reference body.

Remaining differences (after):

| Case | Class | Cause |
|---|---|---|
| `inf:h1_eq4:abaco1_simple`, `inf:h1_eq7:{default,all,abaco1_simple}` | DIFF, validated by SymPy | `integrate` log form (`log(y - 1) + log(y + 1)` for `log(y**2 - 1)`), outside the lot |
| `inf:h1_eq7:bivariate`, `inf:lg11:{all,bivariate}`, `inf:lg03:{default,all,bivariate}` | DIFF, validated by SymPy | SymPy's bivariate/chi generator follows the hash seed (sign, scale, duplicate removal); `lg03` was DECLINED before |
| `inf:uu3:default` | DECLINED (deadline) | SymPy's `default` stops at `abaco2_unique_unknown` (4.7 s). symox raises ValueError there (next row) and goes on to `abaco2_unique_general`, which runs past the deadline on both sides. The old search raised NotImplementedError |
| `inf:uu3:abaco2_unique_unknown` | DECLINED (ValueError) | `solve` gap outside the lot: `solve(eq, df)` gives `sqrt(x**2*(a + x**2 + x*f))/x**2` for SymPy's `sqrt(a + x**2 + x*f)/x`, so the `Pow` atom differs and `cancel(p.diff(y)/p.diff(x))` does not separate |
| `dsolve:h1_eq3` | DIFF | `log(2 + exp(-C1)/(2*x + 1))` written as `log((4*x*e**C1 + 2*e**C1 + 1)*e**-C1/(2*x + 1))`, SymPy's root with `C1 → exp(-C1)`; a direct `checkodesol` gives `(True, 0)`, the corpus check did not conclude; `solve` form gap |
| `dsolve:lg01`, `dsolve:lg03` | DIFF (SymPy past the deadline) | the same RecursionError as SymPy (below) |

Equal only because both sides ran past the 120 s deadline (13 cases, no
parity shown):

- `all`, `abaco2_unique_general` and `checkinfsol` of `sim2`, `uu1`,
  `uu2` and `uu3`;
- `dsolve` of `fsum`.


### lie_group_01 and lie_group_03

`dsolve(hint='lie_group')` of these equations recurses in SymPy:

    pdsolve → pde_1st_linear_variable_coeff → dsolve → lie_group → pdsolve …

It ends in RecursionError after 6 min 22 s (lg01) and 8 min 10 s (lg03).
symox follows the same recursion through the `pdsolve` host call and
raises the same RecursionError after 79 s (lg01) and 50 s (lg03). The
frame counts and timings are in `longrun/`.

### Time

On the 400 SAME cases that finished on both sides:

| | Total |
|---|---:|
| SymPy 1.14 | 772.6 s |
| symox before | 790.9 s |
| symox after | 161.8 s |

Examples, SymPy against symox after:

| Case | SymPy | symox after |
|---|---:|---:|
| `inf:uu1:default` | 3.6 s | 0.095 s |
| `inf:uu2:default` | 13.5 s | 0.22 s (130 s deadline before) |
| `inf:lin1:all` | 61.0 s | 5.3 s |
| `check:lin1` | 66.9 s | 7.2 s |

## Port

| Step | Where |
|---|---|
| `infinitesimals`, nine heuristics, `lie_heuristic_*` | Rust `ode/lie_group.rs` |
| `checkinfsol` (`(True, 0)` with an int `0`) | Rust `ode/lie_group.rs`, `ode_lie.rs` |
| `_lie_group_remove` | Rust |
| `LieGroup._matches` split, `classify_ode` | Rust |
| `_get_general_solution`, canonical coordinates, user `xi`/`eta` | Rust, `pdsolve` as a host call |
| second derivative `factor_terms(signsimp(...))` (`diff2`) | Rust |
| admission, dispatch, hook table, reference body on decline | Python (`lie_group.py`, `single.py`, `ode.py`) |
| `pdsolve`, `pde_1st_linear_variable_coeff` | Python, not ported |
