# Same-session cost checks

The full historical comparison suggested broader regressions. Two consecutive
seven-process screenings check those workloads after the Python suite finishes.
Both use the same CPU, inputs, timing sources and configuration, with alternating
SymPy/native backend order inside each screening. Native versions themselves
are not interleaved, and the shared host remains unisolated.

All five workloads pass both-answer checks in every trial, with unchanged
package/dependency environments. Both packages match their immutable wheel
manifests. Reference medians change by less than 1% between these controls.

Native medians in microseconds; increases remain explicit.

| Workload | Before first | After first | Before repeated | After repeated | Repeated change |
|---|---:|---:|---:|---:|---:|
| `matrices.charpoly_generic.6` | 25471.739 | 26315.980 | 22269.772 | 23280.059 | +4.5% |
| `matrices.jacobian.16` | 159.180 | 171.511 | 109.732 | 115.594 | +5.3% |
| `ode.harmonic` | 239.356 | 259.267 | 88.448 | 103.096 | +16.6% |
| `printing.str.8` | 612.134 | 644.348 | 575.355 | 593.496 | +3.2% |
| `solvers.nonlinear_system` | 275.358 | 299.495 | 102.165 | 100.066 | -2.1% |

The generic 6×6 characteristic polynomial costs about 4.5% more on repeated
calls here, versus about 21% in the historical comparison. Harmonic ODE
repeated calls still cost about 17% more. This is a measured cache/metadata
tradeoff, not a uniformly faster fresh-input algorithm. The full screening
and its original regressions remain in `after.json.gz` and `comparison.json`.
