# checksol (public, native), 2026-10-30

Base 2ea88aaf ("before") against the feat-checksol-20261030j build
("after", 597991e2) and SymPy 1.14. Run on the local machine under a 4G
memory cap.

## Differential tests

`cmp.py cases.json sympy.json symox.json` compares the result or
exception (type and message), the warnings and the solution mapping
after the call (SymPy simplifies it in place).

| set | cases | before agrees | after agrees | after disagrees |
|---|---|---|---|---|
| generated (`gen_cases.py`: 156 bases, 34 flag combinations, every flag) | 5304 | 4166 | 5304 | 0 |
| logged calls of SymPy's test_solvers.py and test_solveset.py (`cslog_public.py`, 1897 distinct) | 1897 | 267 | 1889 | 1 decision, 6 mappings (below) |
| solvegaps pairs (`replay_pairs.py`, internal and public check) | 1414 tests + 1363 ODE | - | 1408 + 1217 | 1 decision (the same as above) |

The pair replays also include 5 + 146 calls that cannot be rebuilt (locally
defined `diffx` in the ODE tests, Dummy-only names, `ExprCondPair` before
the namespace fix): they are not counted.

Before, the generated set had 1074 different decisions and 64 different
warnings. Among them, 201 cases gave False where SymPy gives None, and 24
gave True where SymPy gives None. Every Relational or Boolean case raised
TypeError instead of SymPy's AttributeError. In the logged set, 1627 calls
raised TypeError ("Unknown checksol options", because SymPy ignores the
solve flags that solve passes through). 45 `minimal=True` calls gave False
where SymPy gives None. `checksol([1e-10, x - 2], x, 2)` gave False
(SymPy True), and `checksol([], ...)` gave None (SymPy ValueError).

Remaining differences on the logged set, after:

- `test_issue_17799`: for `-erf(x**(1/3))**pi + I` at
  `x = erfinv(I**(1/pi))**3` (simplify=False), SymPy gives False and symox
  gives None. SymPy's `Function._eval_evalf` for `erf` converts its argument
  inside `try/except ValueError`, so mpmath's "erfinv(x) is defined only
  for -1 <= x <= 1" leaves `erf(...)` unevaluated. The numerical test
  `abs(...) < 1e-9` is then not `S.true`, which gives False. In symox the
  nested erfinv error propagates out of `.n()` and checksol returns None.
  Before this branch it was None too. This is an evalf gap, not a checksol
  one.
- `test_unrad1` (4 calls, simplify=False): the decision is the same and the
  mapping is not changed. Only the printed term order of the input differs
  (the str printer's ordering).
- `test_exclude` (2 calls): the decision is the same. The simplified
  mapping value differs: SymPy's `simplify(Vplus*(C**2*R**2*s**2 + 3*C*R*s +
  1)/(C*R*s))` is `C*R*Vplus*s + 3*Vplus + Vplus/(C*R*s)`, while symox's
  simplify gives `Vplus*(C*R*s + 3 + 1/(C*R*s))`. This is a simplify gap.

## Ten calls

`bench.py` makes each call once (cold, in a process that has already
imported the library), then 5 more times (median, warm). 3 processes per
label. `summarize.py raw` prints the medians below. The raw JSON is in `raw/`.

| call | SymPy cold / warm (ms) | before cold / warm (ms) | after cold / warm (ms) | before r | after r | SymPy r |
|---|---|---|---|---|---|---|
| x**2-4 at 2 | 109.25 / 0.05 | 0.45 / 0.11 | 0.34 / 0.04 | True | True | True |
| Eq(x**2,4) at -2 | 2.77 / 0.09 | 0.19 / 0.13 | 0.12 / 0.05 | True | True | True |
| system 2x2 | 0.84 / 0.16 | 0.49 / 0.34 | 0.15 / 0.08 | True | True | True |
| exp(x)-y at log(y) | 0.93 / 0.05 | 0.20 / 0.12 | 0.11 / 0.05 | True | True | True |
| x*exp(x)-1 at LambertW(1) | 10.83 / 6.62 | 212.70 / 146.76 | 0.55 / 0.36 | True | True | True |
| trig identity sol | 25.85 / 12.69 | 1.50 / 0.92 | 2.91 / 2.18 | True | True | True |
| sqrt(x)+1 at 1 | 0.23 / 0.05 | 0.18 / 0.11 | 0.09 / 0.05 | False | False | False |
| z**x-y minimal | 0.85 / 0.05 | 2.52 / 2.34 | 0.09 / 0.06 | None | True | True |
| Derivative f(x) | 5.68 / 3.36 | 1.61 / 1.34 | 1.54 / 1.34 | True | True | True |
| unverified warn | 22.85 / 7.70 | 3.56 / 3.53 | 2.45 / 2.33 | None | None | None |

"trig identity sol" (`x - 2` at `{x: 2*sin(y)**2 + 2*cos(y)**2}`) is
slower after. The mapping is now simplified with the public `simplify`
(the native check calls it back), as SymPy does. Before, a weaker
simplification was used.

Whole sets, run one after another (sum of the per-call times, seconds):

| set | SymPy | before | after |
|---|---|---|---|
| generated, 5304 calls | 13.25 | 17.75 | 7.01 |
| logged, 1897 calls | 126.23 | 19.29 (1627 of them raised TypeError early) | 4.58 |

## Files

- `gen_cases.py`, `cslog_public.py` (logs the public calls of SymPy's test
  files), `run_cases.py`, `replay_pairs.py`, `cmp.py`, `bench.py`,
  `summarize.py`.
- `raw/bench_<label>_<i>.json`.
- `raw-evidence.tar.gz`: the case sets, every per-label result file, the
  pair replays and the SymPy probe values.

## Validation (597991e2)

- Full suite (demeter, 24 shards): 58434 passed and 5 failed. None of
  the failures is caused by this branch:
  - `tests/test_native_ratintgen.py::test_logpart_gaussian_three_poles`
    fails in the same shard order on the base 2ea88aaf too. It is the
    dag_id traversal pollution fixed on main by 75c27d37.
  - The 4 `tests/test_polymatrix_vector_protocol.py` reload-guard tests
    depend on the hash seed. They pass locally on both trees and are fixed
    on main by a850d52c ("hash-seed independent code key").
- `cargo test -p symox --release`: 822 passed. The doctests ("Couldn't
  compile" under load) were rerun on one thread: 5 passed.
- `cargo clippy -- -D clippy::correctness`: clean.
- ODE corpus `{'cls': 303, 'hint': 279, 'def': 274, 'n': 315}`; integration
  corpus `{"COND": 20, "DIFF": 10, "EQUIV": 8, "MISSING": 1, "SAME": 1171}`
  (both at baseline).
- Caller tests (`tests/test_*solve*.py`, `tests/test_*checksol*.py`,
  `tests/test_solver_checks.py`, `tests/test_equation_solver_helpers.py`):
  1755 passed after the 3 SymPy corrections in `tests/test_solver_checks.py`.
