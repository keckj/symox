# Fraction-free kernels over non-ZZ/QQ domains (wave 18, fflu-dom)

Before: release binary of base commit 34cb18b5 (copied with the same Python
tree into a scratch PYTHONPATH). After: release binary of this branch, built
from the same tree with the same command (`cargo build --release -p symox-py`).

`harness.py --drive BEFORE AFTER results.json` runs every case in fresh
processes: "first" is the median of 5 fresh processes timing the first public
call (import and matrix construction excluded), "repeat" is the median of
7 repeated calls in one further process (after one warm-up call).
Raw per-run values are in `results.json`, the console output in
`run-output.txt`. The machine (8 cores) was shared with other builds: load
average about 5.7 during the run, so absolute times are noisy; the ratios
come from before/after runs interleaved case by case.

Cases: ZZ_I 7x7, GF(10007) 14x14, QQ<sqrt(2)> 7x7, ZZ[x] 5x5 (degree <= 2
entries); public `DomainMatrix.fflu()`, `rref_den()`, `det()`, `charpoly()`
and `Matrix.det(method='bareiss')`.

| case | first before | first after | ratio | repeat before | repeat after | ratio |
|---|---|---|---|---|---|---|
| zzi_fflu | 1.579 ms | 0.855 ms | x1.85 | 1.409 ms | 0.663 ms | x2.12 |
| zzi_det | 1.099 ms | 0.469 ms | x2.34 | 1.009 ms | 0.367 ms | x2.75 |
| zzi_charpoly | 3.809 ms | 0.644 ms | x5.92 | 3.623 ms | 0.463 ms | x7.82 |
| zzi_rref_den | 1.167 ms | 0.652 ms | x1.79 | 0.993 ms | 0.485 ms | x2.05 |
| gf_rref_den | 6.062 ms | 1.661 ms | x3.65 | 5.806 ms | 1.567 ms | x3.70 |
| gf_det | 7.438 ms | 2.136 ms | x3.48 | 7.329 ms | 1.982 ms | x3.70 |
| gf_fflu | 8.379 ms | 2.333 ms | x3.59 | 8.071 ms | 2.121 ms | x3.80 |
| qqa_det | 1.057 ms | 0.408 ms | x2.59 | 0.966 ms | 0.322 ms | x3.00 |
| qqa_rref_den | 1.168 ms | 0.632 ms | x1.85 | 0.984 ms | 0.458 ms | x2.15 |
| zzx_det | 0.572 ms | 0.330 ms | x1.73 | 0.479 ms | 0.228 ms | x2.10 |
| zzx_rref_den | 0.714 ms | 0.527 ms | x1.36 | 0.552 ms | 0.362 ms | x1.52 |
| zzx_charpoly | 3.385 ms | 2.150 ms | x1.57 | 3.208 ms | 1.982 ms | x1.62 |
| zzx_fflu | 1.053 ms | 0.726 ms | x1.45 | 0.886 ms | 0.574 ms | x1.54 |
| matrix_bareiss_zzi | 1.227 ms | 1.245 ms | x0.99 | 1.069 ms | 1.062 ms | x1.01 |
| matrix_bareiss_zzx | 0.541 ms | 0.580 ms | x0.93 | 0.311 ms | 0.311 ms | x1.00 |

`Matrix.det(method='bareiss')` does not reach the DomainMatrix kernels (hit
counters show zero calls); it is unchanged, and the x0.93 first-call value is
within the noise of the shared machine (the repeat is identical).

The remaining ZZ[x] cost is the element arithmetic itself (the kernels call
the same `DomainElement` operations as the Python bodies, so results are
identical by construction); the gain there is the removed Python dispatch.
