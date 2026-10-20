# Public algorithm performance screening

Reference: SymPy 1.14.0. Symox distribution: 0.1.1.
Package fingerprint: `6f23f3f8e68e201204afa14900aeedbb0b6a4ca554d61b9839c748763040550b`.

5 workloads; 5 verified and comparable; 0 excluded.
0 comparable workloads have a first-call or repeated-call speedup ≤ 1.5.

Speedup means reference time / native time; values below 1 mean Symox is slower.
Fresh processes alternate backend order on one CPU. Inputs and result verification are outside timing.
Repeated calls retain both libraries' caches. These measurements are not a universal parity or complexity certificate.
Three-process screenings are provisional; bootstrap intervals describe observed variation, not host isolation.

## Weak or marginal cases

| Workload | First speedup | Repeated speedup | Native first, ms | Native repeated, ms | Processes |
|---|---:|---:|---:|---:|---:|

## Excluded results

These statuses require diagnosis; a validation timeout or alternate formatting is not evidence of an incorrect native algorithm.

| Workload | Status |
|---|---|

## All comparable cases

| Workload | First speedup | Repeated speedup |
|---|---:|---:|
| `matrices.charpoly_generic.6` | 2.526 | 2.708 |
| `matrices.jacobian.16` | 619.296 | 73.713 |
| `ode.harmonic` | 123.298 | 201.478 |
| `printing.str.8` | 25.495 | 21.457 |
| `solvers.nonlinear_system` | 396.250 | 35.997 |
