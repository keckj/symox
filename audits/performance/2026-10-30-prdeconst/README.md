# prdeconst (wave 34): constant_system over rational function fields in Rust

Harness (template: ../2026-10-30-prdefrac): `bench.py` (one fresh process
per call: first call, then the median of 3 repeated calls; prints the sha1
prefixes of `srepr(result)` and `str(result)`), `run.sh <label> <module>
<PYTHONPATH|-> <reps>` (every process under `systemd-run -p MemoryMax=6G -p
MemorySwapMax=0` and `timeout 120`, run from `/tmp`), `summarize.py raw.txt`
(medians, `tables.md`). Labels: `before` = base tree e9e858e5 with its
release `_core` build, `after` = this branch (release build), `sympy` =
SymPy 1.14 in the same venv. The three configurations were interleaved, 5
rounds. `raw-run1.txt` is a first run made while `cargo test` and both
corpora were running on the same machine (noisier, kept for completeness:
first-call outliers up to 1.48x on #5 and #7 that did not reproduce).

`kernel.py` times `constant_system` alone (median of 7 runs of 20 calls) on
the 16 non-`QQ` generated cases of tests/test_native_prdeconst.py: `native`
= this kernel, `stock` = the stock Python body in the same process (the
admission made to decline), `sympy` = SymPy 1.14. Raw: `raw-kernel.txt`.

## Change (computation location)

- Rust, `crates/symox-core/src/integrate/prde_const.rs`: the whole of
  `constant_system` (fraction-field rref via `ddm::field_irref` on
  `ratfunc::RF`, the derivation over the active tower levels, the system
  construction with the stock item expressions, the PolyMatrix domains via
  a `Host` trait).
- Rust, `crates/symox-py/src/prde_const.rs`: the host (domain inference
  through the Rust Poly constructor, `unify_native`, membership), the
  argument conversion and the admission of the values.
- Python, admission only: `_prde_native._constant_system_field`.

## Results

Public `risch_integrate` calls: within noise (0.94x to 1.05x after/before,
first and repeated). `constant_system` is a small part of these calls
(1 or 2 fraction-field systems per integrand, about 1 ms of 13-470 ms).
All results identical to before (srepr); `str` identical to SymPy for all
ten (the srepr difference on #1 and #3 is the known `NonElementaryIntegral`
srepr gap, see prdefrac).

`constant_system` alone: 3.0x to 10.0x faster than the stock body, 2.1x
to 12.4x faster than SymPy (see the table below).

| case | native us | stock us | SymPy us | stock/native | SymPy/native |
|---|---|---|---|---|---|
| 0 | 267 | 980 | 570 | 3.7 | 2.1 |
| 1 | 260 | 993 | 532 | 3.8 | 2.1 |
| 2 | 212 | 745 | 554 | 3.5 | 2.6 |
| 4 | 332 | 1410 | 1424 | 4.3 | 4.3 |
| 5 | 3834 | 11399 | 20922 | 3.0 | 5.5 |
| 6 | 500 | 1931 | 1323 | 3.9 | 2.6 |
| 7 | 373 | 3099 | 3141 | 8.3 | 8.4 |
| 8 | 362 | 1696 | 2024 | 4.7 | 5.6 |
| 9 | 1214 | 7043 | 15097 | 5.8 | 12.4 |
| 10 | 550 | 1991 | 2301 | 3.6 | 4.2 |
| 11 | 407 | 4053 | 4572 | 10.0 | 11.2 |
| 12 | 227 | 960 | 638 | 4.2 | 2.8 |
| 14 | 304 | 1370 | 1709 | 4.5 | 5.6 |
| 16 | 340 | 1380 | 1003 | 4.1 | 3.0 |
| 17 | 301 | 1204 | 776 | 4.0 | 2.6 |
| 19 | 383 | 1702 | 2021 | 4.4 | 5.3 |

Public calls (`raw.txt`):

first call in a fresh process

| # | integrand | before ms | after ms | after/before | SymPy ms | after/SymPy | same srepr as before | same srepr as SymPy | same str as SymPy | n |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | `x*log(x)` | 35.4 | 33.6 | 0.95 | 75.3 | 0.45 | yes | yes | yes | 5 |
| 1 | `log(log(x))` | 143.9 | 144.5 | 1.00 | 149.8 | 0.96 | yes | NO | yes | 5 |
| 2 | `x*exp(x)*log(x)` | 96.3 | 94.4 | 0.98 | 120.7 | 0.78 | yes | yes | yes | 5 |
| 3 | `log(x)*log(x+1)` | 451.0 | 463.3 | 1.03 | 787.6 | 0.59 | yes | NO | yes | 5 |
| 4 | `exp(x)*log(x)` | 92.4 | 93.5 | 1.01 | 104.2 | 0.90 | yes | yes | yes | 5 |
| 5 | `x*exp(x**2)*log(x)` | 97.0 | 94.5 | 0.97 | 102.1 | 0.93 | yes | yes | yes | 5 |
| 6 | `log(x)**2*exp(x)` | 90.9 | 92.5 | 1.02 | 123.2 | 0.75 | yes | yes | yes | 5 |
| 7 | `exp(x)*log(x)/x` | 79.3 | 78.6 | 0.99 | 116.6 | 0.67 | yes | yes | yes | 5 |
| 8 | `log(x)*exp(x)*x**2` | 95.5 | 93.1 | 0.98 | 128.0 | 0.73 | yes | yes | yes | 5 |
| 9 | `log(x)/(x*(1+log(x)))` | 34.9 | 34.4 | 0.98 | 55.9 | 0.61 | yes | yes | yes | 5 |

repeated call (median of 3 per process)

| # | integrand | before ms | after ms | after/before | SymPy ms | after/SymPy | same srepr as before | same srepr as SymPy | same str as SymPy | n |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | `x*log(x)` | 13.1 | 13.0 | 0.99 | 24.1 | 0.54 | yes | yes | yes | 5 |
| 1 | `log(log(x))` | 114.7 | 113.7 | 0.99 | 134.3 | 0.85 | yes | NO | yes | 5 |
| 2 | `x*exp(x)*log(x)` | 75.9 | 72.5 | 0.96 | 66.6 | 1.09 | yes | yes | yes | 5 |
| 3 | `log(x)*log(x+1)` | 429.7 | 437.1 | 1.02 | 770.4 | 0.57 | yes | NO | yes | 5 |
| 4 | `exp(x)*log(x)` | 73.0 | 70.4 | 0.96 | 52.6 | 1.34 | yes | yes | yes | 5 |
| 5 | `x*exp(x**2)*log(x)` | 71.9 | 72.6 | 1.01 | 71.5 | 1.02 | yes | yes | yes | 5 |
| 6 | `log(x)**2*exp(x)` | 70.3 | 69.6 | 0.99 | 67.6 | 1.03 | yes | yes | yes | 5 |
| 7 | `exp(x)*log(x)/x` | 56.1 | 58.8 | 1.05 | 63.0 | 0.93 | yes | yes | yes | 5 |
| 8 | `log(x)*exp(x)*x**2` | 73.5 | 72.7 | 0.99 | 70.6 | 1.03 | yes | yes | yes | 5 |
| 9 | `log(x)/(x*(1+log(x)))` | 15.0 | 14.2 | 0.94 | 27.5 | 0.52 | yes | yes | yes | 5 |
