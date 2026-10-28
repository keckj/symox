# Native probability, floating inversion and gamma kernels

Baseline: `34525336…`. Final release: `c461eb3d…`.
The [full validation record](../../sympy-1.14.0/validation-2026-10-26-kernels.json)
identifies frozen build sources and the complete Demeter run: 45,345 passed,
55 skipped, 6 expected failures, zero unexpected failures on 45,406 collected
nodes. All 171 new regression cases pass. No expected failure was added.

## Probability kernels

The existing benchmark definitions and independent validators are unchanged.
Three fresh processes per backend, same pinned CPU and immutable packages;
construction and result verification are outside timing. Repeated calls retain
caches. These are single-host measurements, not a universal speed guarantee.

| Existing workload | Repeated before, ms | Repeated after, ms | Gain |
| --- | ---: | ---: | ---: |
| Binomial variance | 12.546 | 1.000 | 12.6× |
| Normal fourth moment | 2.176 | 0.595 | 3.7× |

Raw binomial and normal moment formulas now run in Rust. Native polynomial
expectations/variances combine independent law moments without enumerating
all finite worlds in Python. Unsupported expressions and optimization
budgets retain the original general route; floating and conditional queries
retain their established paths. Custom moments, changed tables and lazy or
overridden table properties are not bypassed. Finite zeroth moments retain
the actual table mass. `stats-{before,after}.json.gz` retain timings, outputs,
checks and environment fingerprints.

A later regression probe found an uncovered held-input error in this release:
`X=Die('X',6); E(Mul(X,X,evaluate=False))` returned `7/2` instead of `91/6`.
Polynomial dictionary extraction overwrote repeated generator exponents.
The following algebra batch fixes accumulation and adds a regression. The
benchmark inputs above are canonical, independently checked expressions;
their timings do not establish correctness for every held construction.

## Floating RR/CC inversion

Gauss-Jordan and partial pivoting moved from Python/mpmath to the Rust
floating kernels. Three fresh processes per case; seven repeats of one
matrix and seven changed matrices, with all eight matrices constructed
before timing. Each distinct result is certified outside timing using
exact rationalized binary inputs/outputs and Gaussian-rational DomainMatrix
multiplication. Every entry of the identity residual has magnitude below
`2**(-(precision-10))`. There are 288 matrix certificates per release across
96 distinct inputs, including 53/140 bits and sizes 3/10/20. Both runs use the identical
`inverse-benchmark.py` source hash.

| Size | Domain/bits | Repeated before, ms | Repeated after, ms | Gain |
| --- | --- | ---: | ---: | ---: |
| 3 | RR/53 | 0.232 | 0.0111 | 20.8× |
| 3 | CC/53 | 0.659 | 0.0396 | 16.7× |
| 3 | RR/140 | 0.249 | 0.0144 | 17.3× |
| 3 | CC/140 | 0.649 | 0.0591 | 11.0× |
| 10 | RR/53 | 3.827 | 0.336 | 11.4× |
| 10 | CC/53 | 10.861 | 1.061 | 10.2× |
| 10 | RR/140 | 4.259 | 0.484 | 8.8× |
| 10 | CC/140 | 10.825 | 1.484 | 7.3× |
| 20 | RR/53 | 26.168 | 2.619 | 10.0× |
| 20 | CC/53 | 71.653 | 7.981 | 9.0× |
| 20 | RR/140 | 27.987 | 3.807 | 7.4× |
| 20 | CC/140 | 73.198 | 11.270 | 6.5× |

Changed-input gains range from 6.5× to 18.1×; they are recorded separately from
repeated inputs in `inverse-{before,after}.json`. Cold first calls are also
retained separately. Reference rounding, presentation precision, classes,
singularities and interrupts have dedicated tests. Nonfinite matrix entries
retain their earlier explicit unsupported/error paths; this is not a claim
of general matrix-domain parity.

## Gamma coefficient reuse

`gamma-spouge-benchmark.py` constructs 24 independent held expressions per
family before timing and evaluates them in four batches of six at 30/80
decimal digits. Before/after alternate across three fresh processes each.
Every value is independently checked against mpmath at 80 extra digits;
all value-string hashes also match between the releases and trials.
The cache stores coefficients by full working precision and `a`, never gamma
values or rounded precision keys. Distinct arguments exercise the same
coefficient reuse as repeated arguments.

| Digits | Argument family | Reuse before, ms | Reuse after, ms | Gain |
| --- | --- | ---: | ---: | ---: |
| 30 | Midrange, repeated | 3.873 | 0.128 | 30.3× |
| 30 | Midrange, distinct | 3.857 | 0.140 | 27.6× |
| 30 | Reflection, repeated | 26.245 | 0.410 | 64.0× |
| 30 | Reflection, distinct | 26.015 | 0.416 | 62.5× |
| 80 | Midrange, repeated | 19.163 | 0.289 | 66.2× |
| 80 | Midrange, distinct | 18.984 | 0.289 | 65.7× |
| 80 | Reflection, repeated | 109.670 | 0.959 | 114.3× |
| 80 | Reflection, distinct | 108.976 | 0.930 | 117.1× |

A new precision key still pays coefficient construction. First calls on
fresh keys are approximately 3–6% slower in this measurement; reuse gains
must not be presented as cold-miss gains. Later families in one process can
reuse keys from earlier families, so their first expression is not a fresh
process/key measurement. Full per-family first timings are retained.
The cache is limited to 16 keys and 2 MiB; very high precisions keep the former
streaming computation. No convergent-Gauss shortcut was reintroduced.

## Exact Cartesian integration preparation

The full numerical-coefficient traversal and principal square-root recursion
used by Heurisch moved into Rust. Tests preserve unknown coefficients,
parameters, original object returns, exact complex branches and global
primitive checks. This is a structural port; no speed claim is made here.
The broader Python Heurisch algorithm remains work to port.

`comparison.json` contains the source identity and all summaries. Raw timing
records and runners are retained alongside it; `checksums.json` identifies
every artifact. Absolute cache paths are provenance, not portable inputs.

The unchanged 137-workload screening was also replayed once on this release:
132 cases are verified and comparable, with the same five exclusions as the
previous release. The [complete screen](performance-screen.md) exposes all
results, including slow cases. A one-process screen is provisional and does
not establish a universal improvement or mathematical parity.
