# Native QQ preparation and correlated Gaussian moments

Measured release: `570f1832…`, sources `80f63ba1…`.
[Full validation](../../sympy-1.14.0/validation-2026-10-26-storage.json):
46,132 passed,55 skipped,five existing expected failures,zero unexpected
failures on46,192 collected and completed tests. Rust638+five doctests and
clippy1.99 pass;one timing test remains ignored. These counts are scoped
regression evidence,not global mathematical parity certification.

Three fresh processes per release alternate on Demeter CPU24. The runner
waits for cargo/pytest to finish before timing. Construction and independent
verification remain outside timing. All raw samples, binary identities,
package paths and exact commands are retained. These observations are not
universal speed guarantees.

## Full rational linprog preparation

Baseline `d9b14020…` already has native rational pivots. This release also
prepares stock QQ constraints, variable bounds, free variables and objective
in Rust. Other scalar domains and unsafe protocols keep their complete route.

| Variables | Repeated before,ms | Repeated after,ms | Gain |
| --- | ---: | ---: | ---: |
|10|1.738|0.0638|27.2×|
|30|10.689|0.3103|34.4×|
|60|39.424|1.0460|37.7×|

Every result has the exact primal/dual certificate from the identical runner:
nonnegative variables satisfy all constraints and the primal objective equals
the independent dual bound. There are99 result certificates per release over
12 distinct programs, across three processes. First and changed-input times
are separate in the JSON. The earlier entirely Python pivot implementation
measured about1.38s for the60-variable workload; that historical measurement
is retained in the preceding algebra audit,not substituted for this paired
baseline. No input normalization or feasible-set restriction is skipped.

## Correlated Gaussian recurrence

Baseline `9405dccb…` uses the prior Python recurrence; the new release uses
the native iterative exponent-state graph for guarded ordinary scalar inputs.

| Workload | Before,ms | After,ms | Gain |
| --- | ---: | ---: | ---: |
|Raw pair powers(4,4)|0.2492|0.0482|5.17×|
|Raw triple powers(6,4,2)|1.3290|0.2667|4.98×|
|Raw80-digit Float powers(6,4)|0.3753|0.0650|5.78×|
|Public mixed powers(4,2)|0.3690|0.1952|1.89×|
|Function-coefficient fallback control|19.4559|19.6489|0.99×|

Each process preconstructs18 distinct inputs per case, measured in three
batches of six. All540 independent MGF coefficient checks pass; all output
string digests match across releases and trials,including Float values.
Float accuracy is independently checked against the exact binary-input MGF
at130 digits. Native measured workers forbid importing SymPy; reference
verification is enabled only after all timings complete.

The function-coefficient fallback control has no claimed gain. Public-query
gains apply to the whole release,which also contains other optimizations;
they are not attributed exclusively to this recurrence. Stock Vector/Tensor
dispatchers, custom callbacks, indexing protocols, extended domains and
unproved cases have explicit guards and regression coverage.
