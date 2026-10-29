# integsuite (wave 25): before/after timings

Harness: `bench.py` (one case per fresh process: first call, then the median
of 5 repeated calls after `clear_cache()`), driven by `run.sh`, which runs
3 fresh processes per case under `MemoryMax=6G`. Raw numbers: `raw.jsonl`.

- base: release `.so` built from an unmodified 3d7da725, with the same
  Python tree.
- after: release `.so` from this branch (feat-integsuite-20261029g), same
  machine, same session.

Times in ms. "first" and "median" are the medians over the 3 processes.

| case | base first | after first | base median | after median | ratio | result |
|---|---|---|---|---|---|---|
| cos3_def | 1730.4 | 1739.5 | 1690.83 | 1695.93 | 1.00 | fixed: 67*pi/24 -> pi/12 |
| cos2_def | 1660.6 | 1619.5 | 1618.07 | 1581.67 | 0.98 | fixed: 19*pi/12 -> pi/6 |
| cos2_indef | 1659.3 | 1643.1 | 1592.33 | 1588.96 | 1.00 | fixed: doubled floor terms removed |
| cos1_def | 86.7 | 84.2 | 62.85 | 62.20 | 0.99 | same (pi/3) |
| inv2cos_final_false | 55.9 | 54.5 | 33.19 | 32.38 | 0.98 | fixed: no floor term with final=False |
| x_pow_ny | 0.4 | 0.4 | 0.12 | 0.14 | 1.17 | fixed guard: Ne(n*log(x), 0) |
| two_pow_ny | 0.4 | 0.4 | 0.11 | 0.16 | 1.44 | fixed: Piecewise added (wrong at n = 0 before) |
| three_x_pow_ny | 0.4 | 0.4 | 0.15 | 0.17 | 1.13 | fixed guard |
| x_pow_ny_plus3 | 0.4 | 0.4 | 0.14 | 0.16 | 1.14 | fixed guard, x**3*y branch (correction) |
| x_pow_2y | 0.4 | 0.4 | 0.14 | 0.15 | 1.02 | same |

The trigonometric cases are unchanged within noise. The exponential cases
take 0.02-0.05 ms more because they now build and return the guarded
Piecewise that SymPy 1.14 returns (before, the guard was missing or wrong).
That cost is the price of the correct result.

SymPy 1.14 test files run against symox (failed+errors / tests). The
laplace and meijerint copies have their private-helper imports made
optional (the base and after columns use the same copies):

| file | SymPy 1.14 | base | after |
|---|---|---|---|
| test_integrals | 0/205 | 52/205 | 51/205 |
| test_risch | 0/36 | 8/36 | 8/36 |
| test_heurisch | 1 xfail/33 | 3/33 (+1 xfail) | 3/33 (+1 xfail) |
| test_manual | 0/42 | 0/42 | 0/42 |
| test_meijerint | 0/33 | 13/33 | 13/33 |
| test_transforms | 1 xfail/16 | 6/16 | 6/16 |
| test_laplace | 0/4 | 2/4 | 2/4 |
| test_rationaltools | 0/12 | 0/12 | 0/12 |
