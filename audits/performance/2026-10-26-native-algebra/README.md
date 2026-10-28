# Native operation counts and rational simplex

Baseline release: `c461eb3d…`. Measured release: `9405dccb…`.
The [validation record](../../sympy-1.14.0/validation-2026-10-26-algebra.json)
links the frozen sources and complete Demeter collection: 45,749 passed,
55 skipped, six existing expected failures, zero unexpected failures on
45,810 collected and completed nodes. Rust on Demeter passed 618 unit tests
and five doctests, with one ignored test. Clippy 1.84 passes; Clippy 1.99
finds seven preexisting approximate-constant literals, tracked for the next
portable constant-definition correction. That newer check is not reported
as passing for this release.

Three fresh processes per release run on CPU 24 on the same Demeter host.
Inputs and exact certificates are outside timing. Repeated inputs retain
both libraries' caches; changed inputs are measured separately. After the
initial baseline process, before/after processes alternate. Raw nanosecond
samples, script hashes, binary hashes and package paths are retained.
These measured workloads do not certify universal performance or parity.

## Operation counting

| Input | Changed before, ms | Changed after, ms | Gain |
| --- | ---: | ---: | ---: |
| Scalar expression | 3.4523 | 0.1728 | 20.0× |
| Shared function DAG, depth 8 | 0.3378 | 0.0159 | 21.2× |
| Shared function DAG, depth 12 | 0.4646 | 0.0178 | 26.2× |

All returned counts match independently constructed SymPy expressions. Each
process checks 117 results over 24 distinct expressions. Unchanged inputs
already use a public count cache: repeated calls remain approximately
three microseconds, with no material gain claimed. Native DAG traversal
keeps arbitrary-size operation multiplicities. Host, container and matrix
protocols retain their existing fallback; this is a scalar-kernel port.

## Rational simplex

| Variables | Repeated before, ms | Repeated after, ms | Gain |
| --- | ---: | ---: | ---: |
| 10 | 26.079 | 1.829 | 14.3× |
| 30 | 252.517 | 11.395 | 22.2× |
| 60 | 1384.493 | 40.518 | 34.2× |

Changed-input gains range from 14.8× to 33.1×. Every result is checked with
an exact primal/dual certificate: nonnegative variables satisfy every box
and sum constraint, and the primal objective equals the independent dual
bound. Each process certifies 33 results over 12 distinct programs. No
floating tolerance, relaxed feasible set or known answer bypass is used.

The Rust kernel preserves the two-phase method and Bland pivot ordering,
including equality constraints, redundant rows, infeasibility and
unboundedness. Separate regressions enumerate all vertices of small bounded
rational programs and check a classical cycling example. Only ordinary exact
rational coefficients select this kernel. Other scalar domains retain the
original algorithm; bounds transformation and symbolic linear extraction
remain Python work to port.

## Poisson moments

| Workload | Before | After | Gain |
| --- | ---: | ---: | ---: |
| Raw order 12 | 123.47 µs | 33.94 µs | 3.64× |
| Raw order 36 | 469.78 µs | 86.54 µs | 5.43× |
| Raw order 72 | 1579.31 µs | 183.26 µs | 8.62× |
| Public E(P**12) | 8.064 ms | 2.367 ms | 3.41× |
| Independent mixed fourth moment | 2.804 ms | 0.709 ms | 3.95× |

Three fresh processes per release alternate on CPU 24, with no concurrent
builds during measurements. Each case uses 18 distinct preconstructed inputs
in three batches of six; construction and independent exact verification
remain outside timing. All 648 certificates pass and result digests match.
The normal-moment control changes from 48.96 to 45.67 microseconds (1.07×),
which is not evidence of a general speedup. The raw Poisson algorithm now
builds one Stirling row; no cubic-complexity claim is made about the former
explicit Stirling-number formula. Existing precision, host, callback and
conditional fallbacks are preserved.
