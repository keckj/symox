# Exact binary arithmetic reuse

This round targets repeated incremental expression construction. Its predecessor
is [the weak-case follow-up](../2026-10-20-weak-cases/README.md), branch source
`3da436a`, package source `f301412`. The final production snapshot is `0242a28`.
The reference remains SymPy 1.14.0.

## Preservation argument

Native expression nodes carry a size/eligibility weight for ordinary canonical
exact scalar algebra. This metadata does not enter hashing, mathematical
equality, ordering, printing or serialization. Accepted forms are integers and
rationals, ordinary commuting symbols, sums, products and nonnegative integral
powers. Floats, held forms, dummy/wild symbols, binders, functions and matrices
are ineligible. Names beyond 256 bytes, integer/rational components beyond 256
bits and expressions beyond 1,024 size units retain the original arithmetic.
These limits govern reuse, not accepted mathematical domains or precision.

Completed nontrivial binary additions, multiplications and raw canonical powers
can be reused. Each thread retains at most 2,048 entries and 524,288 size units. These are per-thread
bounds, not a global memory allowance; extra retained graphs are the tradeoff
for avoiding repeated reconstruction.
Numerical operands use exact values as keys; non-numerical operands use physical
node identities, with retained owners. Equal symbols with different original
declarations therefore keep their own presentations. Numerical refinement policy
is included in the key. Only eligible completed outputs are stored, and cache
borrows end before recursive arithmetic begins. Identity results project onto
the current caller's operand, retaining the existing Python identity behavior.

The Python priority/evaluation wrapper additionally performs the indexed tensor
adapter's predicate and invokes the same selected handler directly. It keeps
that predicate even for scalar operands and preserves reflected priority,
conversion hooks, held evaluation, indeterminate-division controls and operand
validation. The outer wrapper does not advertise a bypass of these controls.

Repeated identical work can improve without changing the underlying first-call
algorithm or its worst-case complexity. First and repeated timings remain
separate, with both libraries' caches enabled. Large or unfamiliar inputs still
execute the original arithmetic and may not benefit from reuse.

## Complete screening

The final rerun retains 132 verified workloads out of 137 and the same five
exclusions. It contains 585 paired trials; 45 workloads have seven fresh-process
trials. Sixteen cases have a first-call or repeated-call speedup at most 1.5,
compared with seventeen in the predecessor (this includes modest native wins
and near ties). The complete comparison checks identical inputs, timing sources,
configuration and reference-package fingerprints. Its native timing increases
remain visible alongside improvements.

## Focused results

Seven fresh-process trials per workload confirm that incremental construction
of 100 terms takes a median 203 µs on repeated calls, compared with the
predecessor's 373 µs (seven trials). The median paired SymPy/Symox ratio is
1.80. At 500 terms, repeated calls take 3.77 ms, compared with 4.50 ms previously
(three prior trials); the paired ratio is 141.

First calls expose the cost of tracking eligibility and populating the cache:
100 terms increase from 431 to 643 µs, and 500 terms from 4.71 to 5.98 ms.
Both first calls still exceed SymPy speed by about 100× in this run.
These are workload-specific cache tradeoffs, not a uniformly faster fresh-input
algorithm. The final focus results have all checks verified and an unchanged
package/dependency environment. The initial trial remains separate: its
smaller budget (512 entries, 65,536 size units per thread) regressed repeated
500-term calls to 5.19 ms.

## Costs beyond construction

The historical full comparison suggests broader timing increases. The
[consecutive seven-process cost controls](cost-controls.md) remeasure five
workloads after the Python suite finishes. Every answer verifies, both package
environments stay unchanged, and both wheel manifests match. Generic 6×6
`charpoly` repeated calls increase about 4.5%, harmonic ODE calls about 17%.
Reference medians change less than 1%. Native versions are not interleaved and
the shared host remains unisolated. The original full-screening regressions
remain visible; these controls do not overwrite or replace them.

## Regression qualification

The final optimized build passes 30,058 core Python tests (five documented
skips, no failures or xfails, 15 reference deprecation warnings), including
local Trisplyne integrations. Focused protocol, precision, tensor and arithmetic
regressions pass 1,385 tests. The changed Rust core passes 529 tests with one
pre-existing series timing experiment ignored. No new skip or expected failure
was introduced for the discovered noncommutative gap. Optional extra and
downstream suites were not rerun.

## Artifacts

- `after-package.json`: final immutable wheel/build identity and all package hashes.
- `initial-package.json` and `initial-focus.json.gz`: the separately fingerprinted
  initial cache trial (`ca3b741`), retained to expose first-call overhead and
  medium-size eviction rather than relabeling those measurements as final.
- `final-focus.json.gz`, `final-focus.md`, `final-focus-summary.json`: seven-process
  focused confirmation with the final build. `binary-focus.log.gz` retains progress.
- `after.json.gz`, `after.md`, `after-summary.json`: retained answers, timing
  samples, validations, readable report and machine-readable summary.
- `comparison.json`: predecessor/final summaries with identical input and
  timing-protocol checks; the predecessor's report is referenced rather than
  silently relabeled or duplicated.
- `binary-build.log.gz`, `binary-rust.log.gz`, `binary-regressions.log.gz` and
  `binary-core.log.gz`: optimized build, Rust suite, focused regressions and
  full core Python suite. The full Python run includes Trisplyne and excludes
  `tests/extra`; optional extras and downstream qualification are separate.
- `binary-screening.log.gz`, `binary-refinement.log.gz`: measurement progress.
- `cost-before.json.gz`, `cost-after.json.gz`, their readable reports, summary
  JSON and compressed progress logs: seven-process checks on five workloads.
  `cost-before-package.json` retains the predecessor's manifest;
  `cost-comparison.json` checks the protocol and records both native/reference
  changes. `cost-controls.md` explains the limits of this consecutive comparison.
- `checksums.json`: hashes of every retained artifact except the checksum file.
- `predecessor-domain-probe.json`: the predecessor already rejects
  `Symbol(..., commutative=False)` with `ValueError`. This is now recorded as
  an explicit blocking gap in the [parity contract](../../../docs/algorithmic-parity.md).
  This round does not certify general noncommutative algebra.

The same [runner](../../../benches/README.md) pins serial timing to CPU 0;
tests run on CPU 1, Rust compilation on CPUs 2–3. The shared host is not
isolated. Retained inconclusive checks and timeouts remain excluded from
rankings. The full API audit and downstream prerelease qualification retain
their own snapshots and are not renewed by these measurements.
