# Public algorithm performance screening

Reference: SymPy 1.14.0. Symox distribution: 0.1.1.
Package fingerprint: `3320a259cd9b79511eb133750086d0ca9925fd658fb3c9275617db257927b1ab`.

4 workloads; 4 verified and comparable; 0 excluded.
1 comparable workloads have a first-call or repeated-call speedup ≤ 1.5.

Speedup means reference time / native time; values below 1 mean Symox is slower.
Fresh processes alternate backend order on one CPU. Inputs and result verification are outside timing.
Repeated calls retain both libraries' caches. These measurements are not a universal parity or complexity certificate.
Three-process screenings are provisional; bootstrap intervals describe observed variation, not host isolation.

## Weak or marginal cases

| Workload | First speedup | Repeated speedup | Native first, ms | Native repeated, ms | Processes |
|---|---:|---:|---:|---:|---:|
| `parsing.sympify.50` | 10.707 | 1.498 | 4.97357 | 2.96814 | 7 |

## Excluded results

These statuses require diagnosis; a validation timeout or alternate formatting is not evidence of an incorrect native algorithm.

| Workload | Status |
|---|---|

## All comparable cases

| Workload | First speedup | Repeated speedup |
|---|---:|---:|
| `construction.incremental_sum.100` | 103.131 | 1.804 |
| `construction.incremental_sum.500` | 104.638 | 141.302 |
| `expansion.univariate.30` | 73.660 | 2.396 |
| `parsing.sympify.50` | 10.707 | 1.498 |
