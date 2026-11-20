# Audit index

[README](../README.md) · [SymPy parity audit](sympy-parity-audit.md) · [Algorithm performance](algorithm-performance.md) · [Algorithmic parity](algorithmic-parity.md)

`audits/` keeps the recorded evidence of correctness replays, performance
measurements and source inventories. This page lists the tracked
directories (`git ls-files audits`); the untracked `*.tar.gz` evidence
archives are not described. Each report states its own snapshot, machine and
protocol. Its figures apply to that snapshot only, and none of these reports
certifies global parity with SymPy 1.14. Figures below are quoted from the
named report, with its date.

Unless stated, the main report of a directory is its `README.md`, and the
harness that produced the measurements is stored in the same directory
(`bench.py`, `bench_case.py`, `harness.py`, `run.sh`, `summarize.py`,
`worker.py`, ...), not under `tools/` or `benches/`.

## Source inventories

| Directory | Content |
|---|---|
| `architecture/python-runtime-2026-10-27` | Static review of the Python computations still active behind public entry points (`REPORT.md`, in French); `count_sources.py` counts source lines |
| `architecture/python-runtime-2026-10-28-r9` | Line counts of the accepted R9 snapshot with the same script: 1,113 production files, 136,799 Python and 159,397 Rust lexical SLOC (28 October 2026) |

## SymPy 1.14.0 reference audits

Files directly in `audits/sympy-1.14.0/` (`catalog.json.gz`,
`*-examples.json.gz`, `*-keywords.json.gz`, `*-review.json`,
`source-certificates.json`, `reference-*.json`, `triage-wave*.{json,md}`,
`WORKLOG.md`) are the catalog, documentation-example evidence and reviews
of the parity audit. They are produced by `tools/sympy_parity.py`,
`tools/sympy_stdlib_certificates.py` and the `tools/sympy_*_review.py`
scripts, and are described in [SymPy parity audit](sympy-parity-audit.md).
The `validation-2026-10-2x-*.json` files are the full-suite gate records of
the accepted releases; the directories below link them.

| Directory | Content |
|---|---|
| `documentation-2026-10-26` | Replay of all 2,916 documented objects: 2,719 matching in ordinary-example mode, 2,636 in keyword-call mode, 278 distinct objects with differences to review (26 October 2026) |
| `documentation-triage-2026-10-27` | Display triage of that replay, without a rerun (`set-order-review.json`, `triage-set-order.py`) |
| `corpus-2026-10-26` | ODE and integration corpus replay: 315 ODE inputs from SymPy's `test_single.py`; 1,210 cases of `tools/integration_corpus.py`, 1,164 exact output matches and 36 differences, triaged in the report (26 October 2026) |
| `corpus-2026-10-27`, `corpus-r2-2026-10-27`, `corpus-r3-2026-10-27`, `corpus-r5-2026-10-27`, `corpus-r6-2026-10-27` | Integration replays of the successive accepted releases of 27 October |
| `corpus-r9-2026-10-28`, `corpus-r10-2026-10-28`, `corpus-r11-2026-10-28` | Integration replays of the R9, R10 and R11 releases of 28 October |
| `native-derivative-2026-10-27`, `native-inexact-2026-10-27`, `native-literal-derivatives-2026-10-27`, `native-lrt-limits-2026-10-27`, `native-matrix-rational-2026-10-27`, `native-polymatrix-real-log-2026-10-27`, `native-lll-core-2026-10-27` | Integration reports of native ports accepted on 27 October (derivatives, inexact linear reduction, log_to_real, rational driver, Berkowitz, LLL core) |
| `native-fflu-2026-10-28`, `native-lll-rounding-protocols-2026-10-28`, `native-lu-risch-lll-2026-10-28`, `native-math-ports-2026-10-28`, `native-normalforms-risch-lll-2026-10-28`, `native-sdm-irref-2026-10-28`, `native-solve-risch-lll-2026-10-28` | Integration reports of native ports accepted on 28 October (FFLU, LLL, LU, Risch, HNF, sparse Gauss-Jordan, QQ triangular solve) |
| `native-sdm-irref-private-2026-10-28` | Sparse QQ Gauss-Jordan candidate, stated as not integrated in the accepted release |
| `fflu-gap-2026-10-28` | FFLU counterexamples shared by SymPy 1.14 and symox, with the proved correction (`PUBLIC-REPORT.md`, in French) |
| `hnf-bootstrap-gap-2026-10-28` | HNF bootstrap gap and its R13 resolution (`REPORT.json`, `RESOLUTION-R13.json`) |
| `ratint-growth-2026-10-26` | Native and reference `ratint` growth probe (`probe.py`, JSON results) |

## Downstream projects

`audits/downstream/` holds the featured-project audit produced by
`tools/downstream/audit.py` (see `tools/downstream/README.md` and
[downstream integration](downstream-integration.md)): `summary.json`, one
result file per project under `results/`, the candidate wheel identity
(`candidate.json`), `validation.json`, `source-integrity.json` and packed
evidence under `evidence/`. `summary.json` lists 26 entries: 12
`native_suite_passes`, 7 `native_failure`, 3 `reference_not_clean`, 2
publications that are not testable distributions, 1 with no SymPy code
exercised and 1 excluded by its reference requirement;
`featured_drop_in_ready` is `false`. The files carry no date; they were
added to the repository on 19 October 2026 and last changed on 20 October
2026.

## Performance

Most directories compare a base build with the changed build on public
calls, often also against SymPy 1.14, timing the first call in fresh
processes and the median of repeated calls. Later reports name their wave
and branch, and several state that the machine was shared, so their
timings are indicative.

### 20 October 2026

| Directory | Content |
|---|---|
| `2026-10-20` | 137-workload screening before and after the first patch round; `benches/bench_algorithms.py` and `benches/performance_report.py` (see [algorithm performance](algorithm-performance.md) and `benches/README.md`) |
| `2026-10-20-weak-cases` | Follow-up fixes for the weak cases of that screening |
| `2026-10-20-binary-arithmetic` | Exact binary arithmetic reuse in incremental expression construction |

### 26 October 2026

| Directory | Content |
|---|---|
| `2026-10-26-native-dispatch` | Native dispatch and exact denominator preparation |
| `2026-10-26-native-kernels` | Probability, floating inversion and gamma kernels |
| `2026-10-26-native-algebra` | Operation counts and rational simplex |
| `2026-10-26-native-storage` | QQ preparation and correlated Gaussian moments |
| `2026-10-26-native-products` | Dense matrix multiplication and QQ Gauss-Jordan (DDM) |
| `2026-10-26-native-rational` | Rational reduction (Horowitz-Ostrogradsky selection) |
| `2026-10-26-native-calculus` | Sparse calculus and fraction-free reduction |

### 27 October 2026

| Directory | Content |
|---|---|
| `2026-10-27-native-derivative` | Derivative and algebra ports |
| `2026-10-27-native-lrt` | `log_to_real` |
| `2026-10-27-native-berkowitz` | Dense `DDM.charpoly()` |
| `2026-10-27-rational-selector` | Rational selector controls (`initial/`, `optimized/`, `integrated-r5/`, `profile-r5-linear/`) |
| `2026-10-27-native-auto` | AUTO method selection (main report `integrated-r6r1/PERFORMANCE-R6r1.md`) |

### 28 October 2026

| Directory | Content |
|---|---|
| `2026-10-28-native-lll`, `-native-lll-r9`, `-native-lll-r10` | LLL between accepted releases |
| `2026-10-28-native-lu-r9`, `-native-lu-r11` | Exact LU between accepted releases |
| `2026-10-28-native-risch-r9`, `-native-risch-r10` | Risch polynomial reduction between accepted releases |
| `2026-10-28-native-solve-r10`, `-native-solve-r11` | QQ triangular linear solve (`DDM`/`DomainMatrix`), not `solve`; R11 against R10: 18 fixtures, 1,152 exact certificates, R10/R11 first-call ratios 0.593 to 2.012 (above 1 means R11 is faster) |
| `2026-10-28-native-hnf-r11` | Hermite normal form |
| `2026-10-28-native-sdm-irref` | Sparse QQ Gauss-Jordan |
| `2026-10-28-fflu-dom`, `2026-10-28-lu-dom` | Fraction-free and LU kernels over other domains |
| `2026-10-28-ratint-dom` | `ratint` over Gaussian grounds and EX growth |
| `2026-10-28-smith` | Smith normal form and number-field HNF (`summary.md`) |

The R-series directories (`-r9`, `-r10`, `-r11` and the unsuffixed LLL,
sparse Gauss-Jordan and HNF ones) contain `summary.json`, a regression list
and `raw-files.json`; their harness is in the evidence archive, not in the
tracked tree.

### 29 October 2026

| Directory | Content |
|---|---|
| `2026-10-29-solvesuite` | Public `solve` calls on the paths changed by the solvesuite wave, plus unchanged controls |
| `2026-10-29-tsolverust` | Native `_tsolve` steps and `bivariate` (`summary.md`) |
| `2026-10-29-radicalrust` | Native `unrad`, `_solve_radical` and negative-base steps (`summary.md`) |
| `2026-10-29-solvesetrust` | `solveset` native result checks and modular step |
| `2026-10-29-slowcases` | Nested exponential `solve` and a slow `ratint` case (harness in `harness/`) |
| `2026-10-29-polycap1`, `-polycap2`, `-polycap3` | `count_roots`, modular resultant, cyclotomic products, subresultant PRS, factor and solver limits |
| `2026-10-29-bigexact`, `2026-10-29-biggcd` | Large exact combinatorial values; big-integer gcd and Rational normalization |
| `2026-10-29-domgaps`, `2026-10-29-smith-fx`, `2026-10-29-rrnorm` | Domain element comparisons, Smith form over QQ, GF(p) and QQ[x], RR[y] gcd |
| `2026-10-29-integsuite`, `-prderust`, `-prderust2`, `-ratintgrow`, `-lrt3` | Integration: suite paths, parametric Risch DE, `ratint` growth, `log_to_real` of degree >= 3 |
| `2026-10-29-limheld`, `-limrest3`, `-seriesrest` | Limits and series previously left unevaluated |
| `2026-10-29-hyperx`, `-holo` | `hyperexpand` shift operators; holonomic functions |
| `2026-10-29-arrayrust`, `-arrexprrust`, `-arraysuite`, `-tracearr`, `-idxmvn`, `-spec7` | Array and matrix-expression derivatives, traces, conjugates and scalar functions |
| `2026-10-29-derivsym`, `-w23rust`, `-slotfold` | Symbolic derivative orders, `piecewise_fold`, slot caches |
| `2026-10-29-lambdaprint`, `-lambdarr` | `lambdify` generated source and array evaluation |

### 30 October 2026

| Directory | Content |
|---|---|
| `2026-10-30-solvegaps` | `solve` gaps and native `checksol` (subdirectory `checksol/`) |
| `2026-10-30-checksol` | Public native `checksol` |
| `2026-10-30-smallgaps` | Float Bareiss determinant, `checksol` corrections, `Pow` construction |
| `2026-10-30-exroots` | EX roots and the simplify denominator step |
| `2026-10-30-limgaps`, `-limitsrest` | Limit, Order and series gaps against SymPy 1.14 (`scripts/` in `limgaps`) |
| `2026-10-30-odeclass`, `-odesolve`, `-odepart` | `classify_ode` matches, characteristic equations, particular solutions |
| `2026-10-30-prdefrac`, `-prdeconst`, `-prdelog`, `-ratintgen` | `risch_integrate` through the parametric Risch DE; `ratint` outside QQ/ZZ |
| `2026-10-30-arr2mat`, `-arrindex` | Array-to-matrix and indexed-to-array conversions |
| `2026-10-30-markovrust` | Finite Markov chain computations |
