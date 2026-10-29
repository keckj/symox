# Wave 20 polycap3: dense subresultant PRS over ZZ, factor and solver limits

Date 2026-10-29. Branch `feat-polycap3-20261029b`, base `9e2340a1`.
Machine: 8 CPUs shared with other wave-20 agents (load about 6). Every
process ran from `/tmp` under `timeout` and
`systemd-run --user --scope -p MemoryMax=6G -p MemorySwapMax=0`.

## Harness

- `bench.py <case> <repeats>` prints `case first median digest`.
  `first` is the first call in a fresh process; `median` is the median of
  the repeated calls in the same process. `digest` is the first 16 hex
  digits of sha256(str(result)), or the error.
- `run.sh <python dir> <label>` runs every case in its own process.
  `base` is a snapshot of `python/symox` with the release `.so` built from
  `9e2340a1` in the same target directory. `new` is the worktree build.
- The SymPy 1.14 rows come from the same `bench.py` with
  `from sympy import`.
- `probes/` contains the witness, factor and solver probe scripts.
  `raw/` contains their outputs (`gen_a.*`: small and digest references from
  SymPy, base and new).

Raw files: `raw-base.txt`, `raw-new.txt`, `raw-sympy.txt`.

## resultant(..., includePRS=True) over ZZ (seconds)

| case | base first | base median | new first | new median | SymPy first | SymPy median | digest |
|---|---|---|---|---|---|---|---|
| prs_200_150 (1 rep) | 15.44 error | 15.35 error | 5.62 | 3.96 | 26.50 | 18.67 | 413949b5 all three agree, base "PRS storage" |
| prs_150_200 (1 rep) | 17.05 error | 19.47 error | 5.55 | 3.96 | 25.65 | 17.29 | same |
| prs_120_90 (3) | 4.42 | 4.26 | 0.47 | 0.29 | 2.73 | 1.35 | 1c83f684 |
| prs_60_45 (5) | 0.098 | 0.087 | 0.020 | 0.0093 | 0.45 | 0.12 | 7e4aa2a4 |
| prs_40_30 (5) | 0.015 | 0.012 | 0.0048 | 0.0021 | | | e9ce0b44 |
| subres_60_45 (5) | 0.098 | 0.087 | 0.020 | 0.0094 | | | 0b629bbc (= SymPy digest in tests) |
| poly_prs_60_45 (5) | 0.096 | 0.085 | 0.018 | 0.0082 | 0.090 | 0.023 | fa55955d |
| prs_small (20) | 0.0006 | 0.0002 | 0.0005 | 0.0002 | | | 508e1708 |
| prs_abnormal (20) | 0.0007 | 0.0002 | 0.0006 | 0.0002 | | | 4752325e |
| prs_dense_30 (5) | 0.0022 | 0.0015 | 0.0014 | 0.0007 | 0.077 | 0.0013 | f6a61c3a |

Base and new give identical digests wherever base succeeds. On an unloaded
machine the witness took 3.94 s in Symox and 16.65 s in SymPy
(`probes/a_probe.py`).

## factor_list(x**n + x + 1)

| n | SymPy 1.14 | Symox base | Symox new |
|---|---|---|---|
| 2000 | timeout at 600 s | refusal after 2.68 s (`factorization resource limit: work`) | refusal after 2.49 s, same error |
| 4000 | timeout at 600 s (619.7 s CPU, 173 MB RSS, `raw/b_sympy_4000.txt`) | refusal | refusal |

SymPy does not finish, so the limit is kept.

## Solver and root limits (`probes/c_*.py`, timeout 180-600 s)

| probe | Symox | SymPy 1.14 | note |
|---|---|---|---|
| nested exp, depth 10 | 6.80 s | 0.40 s | same result, performance gap |
| nested log, depth 10 | 0.01 s | no result within 300 s | |
| nested sqrt, depth 10 | 0.00 s | 0.13 s | same |
| exp(e)+1 nested, depth 10 | 426.08 s | 77.70 s | same up to log(4) vs 2*log(2) (pre-existing); about 3.6x per level |
| exp(e)+1 nested, depth 8 | 19 s | 1.85 s | |
| rational root 2**41 (cubic, quintic, QQ leading) and roots() | 0.00 s | 0.03-0.04 s | same roots |
| ODE from those polynomials (order 3, 5) | 0.03-0.05 s | no result within 300 s | |
| solve(x**n*exp(x) - 2), n = 9, 64, 256 | 0.28, 1.95, 7.97 s | 1.88, 13.0, 58.9 s | n+1 solutions for even n |
| same, n = 4096, 4097 | > 180 s | > 180 s (extrapolated > 1000 s) | |
| solve(x**4097 - exp(x)) | 0.05 s | | |

No limit refuses an input that SymPy computes under the same caps.

## Validation

- `cargo test -p symox --release`: 781 passed, 1 ignored. The one failure
  is the doctest of `limits/mod.rs` ("Couldn't compile the test"). It is
  pre-existing (seen in the integ18 validation and other branches), and the
  file is unchanged since 34cb18b.
- `cargo clippy -p symox -p symox-py --release -- -A clippy::all -D clippy::correctness`: clean.
- Python suite over 24 shards: 56721 passed. The three known flaky files
  failed in their shards (`test_polynomial_source_storage.py` TimeoutExpired,
  `test_native_risch_polynomial_reduce.py:312`,
  `test_polymatrix_presentation_protocol.py`). They pass locally: 8, 190 and
  44 passed.
- ODE corpus `{'cls': 303, 'hint': 279, 'def': 274, 'n': 315}` and
  integration corpus
  `{"COND": 20, "DIFF": 10, "EQUIV": 8, "MISSING": 1, "SAME": 1171}` are
  equal to the baseline. The integration corpus is also identical per case
  (1210 cases).
