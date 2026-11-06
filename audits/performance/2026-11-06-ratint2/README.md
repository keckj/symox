# ratint2: RootSum constructor, doit and evalf in Rust

- Harness: `diffcorpus.py sympy|symox OUT.jsonl [cases...]` (75 cases:
  RootSum constructions, errors, `args`, `doit`, `evalf`, `diff`, `subs`;
  `integrate`/`ratint` with algebraic and parametric coefficients;
  `ratint_logpart` over `QQ<sqrt(2)>`, `QQ<sqrt(2)>[p]` and
  `QQ<sqrt(2)>(t)`; `apart(full=True)`; EX elements). Each case is run
  once for warm-up and then timed as the best of three runs. A case that gives no result within 60 s
  is recorded as `CaseTimeout`. `compare.py ref.jsonl X.jsonl ...`
  compares the printed results; for errors it compares the exception class.
- `raw/sympy.jsonl`: SymPy 1.14 (`~/symox-ci/venv`). `raw/before.jsonl`:
  main 5038ef02. `raw/after.jsonl`: this branch (release build). Table:
  `raw/compare.md`.
- `rscount.py` wraps the five native functions and counts the calls that
  return a value (`native`) and those that decline (`None`, after which the Python body runs).
  The counts are in `raw/declines.txt`:
  - corpus: 237 transforms, 228 builds, 4 `new`, 16 `doit` and 8 `evalf`,
    all native, with 0 declines;
  - the 34 test files that use RootSum, ratint or apart (3535 passed): 0 declines, apart from the
    two deliberate ones in `test_native_functions_decline_without_poly`.

Result: 62/75 identical to SymPy (before: 61/75).

| group | cases | SymPy ms | before ms | after ms |
|---|---:|---:|---:|---:|
| RootSum constructor | 28 | 46.7 | 49.9 | 13.3 |
| RootSum doit / evalf | 5 | 14.0 | 10.5 | 9.1 |
| integrate / ratint | 9 | 1656.8 | 146.0 | 144.4 |
| apart(full=True) | 8 | 46.6 | 19.9 | 17.4 |

Each group covers the cases that succeed in all three runs.

Changed: `RootSum(x**3 + x + 1, Function('f'))` raised ValueError; it now
gives SymPy's `RootSum(x**3 + x + 1, Lambda(x, f(x)))`.

The 13 differences left:

- 10 `ratint_logpart` cases over `QQ<sqrt(2)>[p]` and
  `QQ<sqrt(2)>(t)`: NotImplementedError at Poly construction, both before and
  after. SymPy returns pairs over `ZZ[p,_t]`/EX. This is a gap, documented in
  `docs/algorithmic-parity.md` (ratint2).
- `EX(sqrt(2)*(1+sqrt(2)))`: symox gives `EX(sqrt(2) + 2)`, while SymPy keeps the
  argument. This is a gap.
- `RootSum x gen keyword`: an artefact of SymPy's cache. A fresh SymPy prints
  `RootSum(y**3 + y + 1, Lambda(x, log(x)))`, as symox does.
- `integrate 1/(x**3+sqrt(2)*x+1)`: SymPy times out (more than 60 s); symox
  gives a result in 61 ms. Not compared.

Every other case, errors included, gives the same printed result as
SymPy, before and after. Timings are noisy (other agents were building).
