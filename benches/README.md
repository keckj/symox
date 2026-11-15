# Algorithm performance and regression checks

`bench_algorithms.py` compares 137 deterministic public API workloads against
SymPy 1.14.0. Its purpose is to identify and fix weak native paths. A result
must pass the independent mathematical and output-quality checks before its
timings enter the comparable ranking. Selected checks do not certify a whole
function, its domains, or asymptotic complexity.

Use an optimized, immutable package extracted from a wheel. This is especially
important while other worktrees are rebuilding the editable installation.
Keep the package directory unchanged throughout a run.

```sh
python benches/bench_algorithms.py \
  --symox-package /path/to/extracted/symox \
  --out /path/to/screening.json.gz \
  --work-dir /path/to/screening-workers \
  --rounds 3
python benches/performance_report.py /path/to/screening.json.gz \
  --out /path/to/screening.md --summary-json /path/to/summary.json
```

The default report highlights cases with a first-call or repeated-call speedup
at most 1.5. Speedup is SymPy time divided by Symox time: below one means Symox
is slower. First-call and repeated-call results remain separate, since the
reference caches many calculations. A fast first call does not excuse an
avoidable repeated-call slowdown.

Each trial/backend runs in a fresh process, serially, on the same CPU. Backend
order alternates; Python GC stays enabled for both. Input construction, imports,
export and result checks are outside the measured operations. Numerical answers
must already be evaluated before export. Native measurement forbids importing
SymPy. The exporter and independent validator may explicitly import it after
timing. Both the first and final answers are retained and checked.

The evidence records backend/dependency fingerprints, timing-source hashes,
individual samples, setup and measurement phases, memory high-water marks,
environment data, actual answers and their hashes. Timeouts remain explicit;
they are never replaced with attractive timing estimates. The total worker
timeout includes startup/setup/export, so inspect the recorded phase before
attributing a timeout to the operation. Concurrent activity on the shared host
can still affect timings.

Extend marginal cases to seven fresh-process trials without changing inputs,
packages or the measurement protocol:

```sh
python benches/bench_algorithms.py \
  --symox-package /path/to/extracted/symox \
  --out /path/to/screening.json.gz \
  --work-dir /path/to/screening-workers \
  --resume --rounds 7 --refine-below 2
```

Use a new output for a patched package. `--case` accepts exact IDs or prefixes
and can be repeated. Resume refuses changed timing sources or package identity.
For a validator correction, recheck the retained answers without measuring
again, keeping the original report and appending the old/new validation history:

```sh
python benches/performance_report.py /path/to/screening.json.gz \
  --revalidate --case solvers.solveset_real \
  --out /path/to/rechecked.json.gz
```

The validator distinguishes mathematical differences, reduced precision,
inexact replacements of exact answers, missing roots, held computations,
inferior factorization/simplification output, alternate formatting, and
inconclusive checks. An inconclusive reference or validation timeout does not
prove a native algorithm incorrect.

After a weak path is patched, run its behavioral regressions and remeasure it.
Preserve domains, exact arithmetic, requested precision, options and explicit
failures. See the [parity contract](../docs/algorithmic-parity.md).

The [performance fix report](../docs/algorithm-performance.md) describes the
2026-10-20 patches, retained evidence and remaining weak paths. Later targeted
measurements of native ports are kept under `audits/performance/`; most use
their own harness, documented in each directory's README, rather than this one.

`bench_vs_sympy.py` (with SymEngine when installed) and the per-area
`bench_*.py` scripts are quick timing comparisons; they do not perform the
validation described above. `audit_solver_parity.py` is a behavioral audit,
not a benchmark (see [testing](../docs/testing.md)).
