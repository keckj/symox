# Public algorithm performance screening

Reference: SymPy 1.14.0. Symox distribution: 0.1.1.
Package fingerprint: `3320a259cd9b79511eb133750086d0ca9925fd658fb3c9275617db257927b1ab`.

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
| `matrices.charpoly_generic.6` | 2.438 | 2.582 |
| `matrices.jacobian.16` | 578.559 | 69.995 |
| `ode.harmonic` | 110.280 | 172.613 |
| `printing.str.8` | 24.289 | 21.036 |
| `solvers.nonlinear_system` | 368.837 | 37.290 |
