# Public algorithm performance evidence

Reference: SymPy 1.14.0, Python 3.10. The patch branch starts at `0d405b7`.
See the [fix report](../../../docs/algorithm-performance.md) and
[runner instructions](../../../benches/README.md).

## Artifacts

- `before.json.gz`: the original 137-workload screening, retaining the original
  timings and an explicit validation correction for bound-variable names in
  `solveset`. Each selected workload has three fresh processes, except an early
  setup failure, which has one.
- `before-package.json`: hashes of all files in the reference native wheel and
  its production/build source commits. This is the pre-patch Symox package,
  not the SymPy reference.
- `before-add.json.gz` and `before-add-package.json`: an intermediate screening
  after the other fixes, before the last exact `Add` reuse patch. This separately
  exposes the effect and eligibility of that final change.
- `after.json.gz`, `after.md`, `after-summary.json` and `after-package.json`:
  the final screening, readable report, machine-readable summary and immutable
  package manifest. All package Python sources were checked byte-for-byte
  against the recorded source snapshot when extracting the wheel.
- `comparison.json`: per-case before/after summaries and checks that inputs,
  timing sources and timing configuration are identical. The final run contains
  501 paired trials; 24 marginal cases have seven fresh-process trials each.
- `nroots-permutations-seven-processes.json.gz`: an earlier seven-process
  confirmation of the precision fix and permutation changes, with its own
  immutable package fingerprint retained in the report.
- `mersenne-seven-processes.json.gz`: an earlier seven-process confirmation of
  the exact perfect-power shortcut, also retaining its own fingerprint.
- `qualified-core.log.gz`, `qualified-regressions.log.gz`, `rust-mersenne.log.gz`
  and `qualified-build.log.gz`: final Python suite, focused regressions, Rust
  core suite and release-build output. The final Add and NumPy binding changes
  do not change the Rust core sources covered by that Rust log. Its one ignored
  test is the pre-existing `series::tests::timing_series` timing experiment.
- `qualified-screening.log.gz` and `qualified-refinement.log.gz`: the final
  three-trial screening and additional trials for cases initially at most 1.5×
  faster than SymPy in either timing mode.

The JSON gzip reports contain actual answers (interned by content hash), first
and repeated timings, samples, status/phase information, validator outcomes,
source hashes and dependency fingerprints. A report's package fingerprint
belongs to that report's snapshot; historical confirmations are not silently
relabeled as measurements of the final package.

## Reproduction

Extract the optimized wheel matching the desired manifest into an immutable
directory. For example, with the final package directory as `$SYMOX_PACKAGE`:

```sh
python benches/bench_algorithms.py \
  --symox-package "$SYMOX_PACKAGE" --out after.json.gz \
  --work-dir after-workers --rounds 3 --cpu 0
python benches/performance_report.py after.json.gz \
  --out after.md --summary-json after-summary.json
```

Use `PYTHONPATH` pointing to the extracted wheel's parent when running pytest.
The full core command is `python -m pytest -q -ra tests --ignore=tests/extra`.
The recorded run enables the local Trisplyne integration with `TRISPLYNE_PATH`
and `TRISPLYNE_PRINTERS`, and sets `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1`.
Its tests run on CPU 1, separately from serial benchmark workers on CPU 0.
The shared machine is not an isolated performance host.

Do not combine these results with the older downstream prerelease qualification
or claim full per-function mathematical parity. Excluded validations and worker
timeouts remain excluded from timing rankings.
