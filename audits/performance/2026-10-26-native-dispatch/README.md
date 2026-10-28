# Native dispatch and exact denominator preparation, 2026-10-26

This compares the previously validated `8323edda…` release with the final
`34525336…` release. The [validation record](../../sympy-1.14.0/validation-2026-10-26-native.json)
identifies the exact working sources and full Demeter test execution.
This batch moves denominator preparation into Rust and removes unnecessary
Python domain selection before existing native matrix algorithms. It does
not complete the port of the remaining Python mathematical algorithms.

## Confirmed public-API timings

Three fresh processes per backend and workload, alternating backend order,
on the same pinned CPU (Intel i7-9700K). Inputs and independent mathematical
result validation are outside timing. Immutable package copies retain their
source and native fingerprints. First calls and repeated identical calls
are separate; the latter retain each library's caches. Times are medians
in milliseconds; gains compare the previous Symox release with this one.

| Workload | First before | First after | Gain | Repeated before | Repeated after | Gain |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `simplification.radsimp.3` | 26.442 | 0.626 | 42.2× | 12.373 | 0.388 | 31.9× |
| `simplification.radsimp.5` | 25.060 | 0.186 | 134.8× | 10.860 | 0.0486 | 223.3× |
| `matrices.inverse_rational.4` | 14.098 | 0.0471 | 299.1× | 0.365 | 0.00645 | 56.7× |
| `matrices.inverse_rational.10` | 16.395 | 0.324 | 50.6× | 2.335 | 0.181 | 12.9× |
| `matrices.charpoly_generic.3` | 15.823 | 1.315 | 12.0× | 1.012 | 0.420 | 2.4× |
| `matrices.charpoly_generic.4` | 16.641 | 1.971 | 8.4× | 1.918 | 0.960 | 2.0× |

The [exact comparison](comparison.json) and retained
[before](performance-targeted-before.json.gz) / [after](performance-targeted-after.json.gz)
records contain all trials, outputs and validation results. This is a
single-host measurement with a powersave governor, not hardware isolation
or a universal speed/complexity guarantee. The benchmark and validation
sources are unchanged between runs.

`radsimp` previously performed Python traversal and numerical recognition
through mpmath identify/PSLQ even for provably irrational square-root sums.
Rust now proves their square classes exactly, preserving recognition for
unproved forms. The five-radical default case still respects `max_terms=4`:
its exact result remains unchanged. Matrix dispatch proves rational or
rational-polynomial entries natively; unknown and overridden domains retain
their general routes. The cold profiles are retained as
`radsimp-{before,after}.txt` and `charpoly-{before,after}.txt`; profile timings
are diagnostic and are not the timing samples used in the table.

## Full existing screening

The unchanged 137-workload catalog was run once per backend before and
after. Both runs have 132 verified, comparable workloads. These one-process
screenings are provisional; they are included to expose the remaining
slow and unproved cases, not to claim that every operation improved.
The final screen has 33 workloads whose first or repeated SymPy/Symox
ratio is at most 1.5, including ratios below one where Symox is slower.

Five exclusions are identical in both runs: undecided `assumptions.trig_zero`,
validation timeout for `differentiation.composite.8`, and non-comparable
`integration.rational.8` / `series.transcendental.20` / `.40`. None is used
as a speed claim. Raw results and statuses are preserved in
`performance-{before,after}.json.gz` and `performance-{before,after}-summary.json`.

## Hypergeometric correctness and regression controls

`hyper-unit-benchmark.py` uses 60 distinct preconstructed inputs per family,
in five batches, with `evalf(30)` outside construction. The final native
unit-pole evaluation returns all 60 positive 1F0 infinities, all 60 positive
balanced infinities and all 60 negative signed-balanced infinities. The old
release retained two families and returned false finite values on 31 of the
60 positive balanced inputs. These are capability/correctness gains; no
speed ratio against an uncomputed or incorrect answer is published.

An intermediate native convergent-Gauss shortcut was rejected: three
Spouge gamma evaluations took about 11.56 ms instead of about 0.285 ms.
`hyper-unit-rejected-gauss.json` retains that experiment's `f59ce49a…`
fingerprint. That shortcut is absent from the final release; the existing
convergent numerical path and precision are preserved.

The isolated regression control alternates three fresh processes per
snapshot, with 60 distinct inputs per family and identical value hashes:

| Comparable family | Before | Final |
| --- | ---: | ---: |
| Gamma | 25.973 ms | 26.087 ms |
| Convergent Gauss | 289.54 µs | 294.61 µs |
| Terminating hypergeometric | 312.90 µs | 311.16 µs |

The small differences do not establish a general improvement or distinguish
proof overhead from measurement noise. The rejected 40-fold regression is
absent. `hyper-unit-controls.py`, the six `hyper-unit-controls-isolated-*.json`
process records and their summary retain the controls.

## Reproduction

Run the existing harness against each immutable package directory:

```bash
.venv/bin/python benches/bench_algorithms.py \
  --symox-package /path/to/snapshot/symox \
  --out targeted.json.gz --work-dir targeted-workers --rounds 3 \
  --case matrices.inverse_rational.4 --case matrices.inverse_rational.10 \
  --case matrices.charpoly_generic.3 --case matrices.charpoly_generic.4 \
  --case simplification.radsimp
```

Omit the case selectors and use `--rounds 1` for the complete screening.
`comparison.json` records the package environments and final source identity;
`checksums.json` identifies every retained artifact.
