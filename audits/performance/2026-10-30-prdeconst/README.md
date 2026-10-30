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

## Item 2: derivation(p, DE) with K(x)[t] Poly coefficients

Computation location: Rust, `crates/symox-py/src/prde_const.rs`
(`derivation_poly`): the stock loop over the active `(D, T)` pairs with the
Rust Poly constructor `Poly._from_normalized` called directly for the two
Poly constructions, and the Rust Poly/Expr `diff`, `as_expr`, `*`, `+`.
Python, admission only: `_prde_native.derivation` (hook at the top of
`risch.derivation`, `basic=False`).

`derivation.py <symox|sympy>` times `derivation` alone (best of 5 runs of
2000 calls). `head` = this branch at item 1 (stock Python body), `new` =
item 2, `sympy` = SymPy 1.14. Raw: `raw-derivation.txt`.

| case | coefficientD | stock us | native us | SymPy us | stock/native | SymPy/native |
|---|---|---|---|---|---|---|
| tower ZZ(x,t0) | False | 967 | 531 | 1880 | 1.82 | 3.54 |
| tower ZZ(x,t0) | True | 630 | 329 | 1210 | 1.92 | 3.68 |
| exp ZZ[x] | False | 535 | 240 | 661 | 2.23 | 2.76 |
| exp ZZ[x] | True | 219 | 72 | 404 | 3.04 | 5.60 |
| prim ZZ(x) | False | 656 | 341 | 1003 | 1.92 | 2.94 |
| prim ZZ(x) | True | 327 | 161 | 582 | 2.04 | 3.62 |

Public calls, `raw-item2.txt` (same harness, `before` = base e9e858e5,
`after` = items 1 and 2, 5 interleaved rounds; run while the ODE corpus was
running on the same machine): after/before 0.90-0.95 first call and
0.89-0.93 repeated call, every result identical to before (srepr) and to
SymPy (str). Repeated call after/SymPy: `x*exp(x)*log(x)` 0.98,
`exp(x)*log(x)` 1.26, all others 0.50-0.96.

## Profile of the two repeated calls slower than SymPy

`profile_risch.py <module> <integrand> <function prefixes>` (cProfile over
5 repeated calls after a warm-up; cProfile inflates Python-heavy code, the
ratios are indicative). Before item 2, `exp(x)*log(x)`: 0.594 s profiled,
of which `risch.splitfactor` 0.451 s (60 calls, 12 per integral) and inside
it `_risch_splitfactor_native.reduce` 0.406 s: `_risch_native._prove_fields`
0.222 s (`all(...)` over the `PolyMatrix` class snapshots,
`_polymatrix_native._class_unchanged` 0.143 s, and
`rationaltools._stock_ratint_protocol_with_metadata` 0.042 s),
`_states_match` 0.082 s, `_str_keyed` (51660 calls) 0.078 s; the native
kernel `_core._risch_splitfactor` itself is below the profiler resolution.
SymPy: `splitfactor` 0.082 s of 0.551 s. `derivation`: symox 0.081 s, SymPy
0.120 s. `x*exp(x)*log(x)` has the same profile (`splitfactor` 0.453 s of
0.610 s; SymPy 0.083 s of 0.687 s). So the time goes into the per-call
verification of the splitfactor protocol snapshot (Python admission), not
into the computation: about 5.7x SymPy's whole `splitfactor` on these
inputs. Not changed here (other wave's admission module); recorded as a gap
in docs/algorithmic-parity.md.

## Item 3: is_deriv_k(fa, fd, DE) in Rust

Computation location: Rust, `crates/symox-py/src/prde_deriv.rs` (new): the
whole stock body. `Df/f` uses the native `derivation` loop (item 2) and the
Rust Poly `*`, `-`, `cancel(include=True)`. The `E_part + L_part` row
(`D[i].quo(Poly(T[i], T[i]))`, `D[i]`) is built with the Rust Poly methods,
and the PolyMatrix domains are inferred over `dum`. `constant_system` is the
item 1 fraction-field kernel (refactored into `field_system`). The checks
`derivation(u_i, basic=True) == 0` and `u_i.is_Rational` run exactly on
`RF` values over the active levels. The answer `(ans, result, const)` uses
the native `Add`/`Mul`/`Pow` constructors, `as_numer_denom`, the Rust
`sqf_list` and `cancel`. The empty row (a tower without exp/log extension)
returns None, as the stock `constant_system` (`not A`). Python, admission
only: `_prde_native.is_deriv_k` (hook at the top of `prde.is_deriv_k`). It
creates the Dummy once, and the stock body reuses it when the kernel
declines. The admission also accepts the stock AccumBounds coefficient fold
(`_compat._accum_fold`, installed once an AccumBounds exists) when it is the
captured function with its code, globals and builtins. The kernel only
builds host-free expressions, for which that fold returns None.

`deriv_k.py <symox|sympy> [stock]` times `is_deriv_k` alone (median of 7
runs of 200 calls). `stock` is the stock Python body in the same process.
Raw: `raw-deriv-k.txt` (run while both corpora ran on the machine).

| case | native us | stock us | SymPy us | stock/native | SymPy/native |
|---|---|---|---|---|---|
| two logs (SymPy sym1) | 1985 | 3909 | 16793 | 1.97 | 8.46 |
| log+exp (SymPy sym2) | 1688 | 4109 | 11529 | 2.43 | 6.83 |
| two logs, rational | 2504 | 4593 | 18100 | 1.83 | 7.23 |
| log, None | 529 | 1817 | 4220 | 3.43 | 7.98 |
| log, x**3 | 954 | 1758 | 2801 | 1.84 | 2.94 |
| no extension | 277 | 359 | 929 | 1.30 | 3.36 |

Public calls, `raw-item3.txt` (`before` = item 2, i.e. HEAD 67c742e6 Python
with the same `_core`: the new kernel is only reached through the new
admission; `after` = item 3; 3 interleaved rounds, while both corpora ran):
after/before 0.95-1.01 first call and 0.97-1.02 repeated call (noise:
`is_deriv_k` runs a few times per integrand, inside the tower construction).
Every result is identical to before (srepr) and to SymPy (str). The known
srepr gap on #1 and #3 is unchanged.

`is_log_deriv_k_t_radical_in_field` is not ported (see
docs/algorithmic-parity.md): it needs `splitfactor`, `residue_reduce` with
its `RootSum`/`real_roots` residue terms, `residue_reduce_derivation`,
`parametric_log_deriv` and a recursion on the tower, with `ValueError` and
`NotImplementedError` paths for each case. It remains the stock Python
body.

The splitfactor admission caching was not done. No registration generation
or change notification exists for the class dictionaries and function
records that `_prove_fields` and `_states_match` compare. Python 3.10 has no
dict or type watchers, and the abi3 build cannot read `tp_version_tag`. A
cache keyed on anything else would miss a class or method mutation, so it
would weaken the protocol checks.
