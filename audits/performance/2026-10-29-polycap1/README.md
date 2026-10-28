# Wave 19 polycap1: `count_roots` Sturm chain and factorization work estimate

Base: commit 58e50206 (wave-18 integration) and its release `_core.abi3.so`.
New: the same tree plus this branch, release build. Reference: SymPy 1.14
from `/home/keckj/Documents/symox/.venv`, run from `/tmp`.

Harness: `bench.py <case> [repeats]` runs one case in a fresh process (the
"first" time is the first call in that process, "rep" is the median of the
following calls; 1 repeat for the 5-10 s cases, 2 otherwise). `run.sh`
runs every case under `timeout` and `systemd-run -p MemoryMax=6G
-p MemorySwapMax=0`; `table.py` renders the raw JSON lines
(`raw-base.jsonl`, `raw-new.jsonl`, `raw-sympy.jsonl`). The machine was
shared with other agents (8 cores); base runs overlapped a `cargo test`.
Ratios are first-call ratios. A result is only a speed-up when it equals
SymPy's (last column; for base the factor cases raised the resource-limit
ValueError, recorded as `error`).

## Changes measured

* `count_roots` over ZZ/QQ with real bounds: the QQ Sturm sequence is
  replaced by positive integer multiples of the same elements (primitive
  pseudo-remainders with lazily scaled coefficients, homogeneous integer
  evaluation). Same signs, same variation counts, same result.
* Checked factorization: dense divisions charge their actual operation
  count instead of `len(a)*len(b)`. The limits themselves are unchanged.

## Results (seconds)

| case | base first | new first | sympy first | base rep | new rep | sympy rep | base/new | sympy/new | same result as sympy |
|---|---|---|---|---|---|---|---|---|---|
| count_5000 | 68.975 | 0.012 | 9.214 | 69.561 | 0.012 | 9.229 | 5649.1x | 754.6x | True |
| count_5000_ab | 69.997 | 0.013 | 9.192 | 70.490 | 0.013 | 10.185 | 5336.1x | 700.7x | True |
| count_5000_rat | 74.024 | 0.018 | 10.570 | 82.585 | 0.018 | 10.266 | 4154.8x | 593.3x | True |
| count_5000_float | 74.587 | 0.017 | 10.883 | 71.495 | 0.017 | 11.110 | 4371.4x | 637.8x | True |
| count_3000_sparse | 9.792 | 0.154 | 4.879 | 9.825 | 0.153 | 4.895 | 63.6x | 31.7x | True |
| count_cheb120 | 0.007 | 0.009 | 0.073 | 0.007 | 0.008 | 0.022 | 0.8x | 8.4x | True |
| count_cheb120_ab | 0.016 | 0.010 | 0.098 | 0.016 | 0.010 | 0.049 | 1.6x | 9.6x | True |
| count_poly_sq | 0.003 | 0.001 | 0.031 | 0.002 | 0.001 | 0.004 | 2.1x | 25.3x | True |
| count_complex_60 | 0.323 | 0.322 | 1.787 | 0.323 | 0.322 | 1.782 | 1.0x | 5.6x | True |
| factor_300 | error | 0.138 | 7.264 | error | 0.135 | 7.090 | - | 52.8x | True |
| factor_200 | 0.267 | 0.253 | 2.259 | 0.270 | 0.251 | 2.233 | 1.1x | 8.9x | True |
| factor_250 | 0.025 | 0.022 | 5.484 | 0.025 | 0.022 | 5.582 | 1.2x | 250.2x | True |
| factor_301 | error | 0.127 | 8.996 | error | 0.126 | 8.992 | - | 70.7x | True |
| factor_expr_300 | 0.161 | 0.140 | 6.776 | 0.162 | 0.138 | 6.852 | 1.2x | 48.6x | True |
| factor_sparse_256 | 0.047 | 0.041 | 3.956 | 0.047 | 0.040 | 4.020 | 1.2x | 97.0x | True |
| factor_sparse_210 | 0.028 | 0.025 | 1.765 | 0.028 | 0.025 | 1.743 | 1.1x | 69.2x | True |

Larger inputs remain refused by the unchanged 2*10**9 work limit after a
bounded time (new build): `factor_list(x**2000 + x + 1)` 3.3 s and
`factor_list(x**4000 + x + 1)` 3.0 s, both "factorization resource limit:
work" (SymPy 1.14 already needs about 30 s at n = 350/400). New build:
`x**600 + x + 1` 1.3 s and `x**1000 + x + 1` 2.1 s, irreducible
(n = 0, 1 mod 3, Selmer).

`count_cheb120*` and `count_poly_sq` exercise the new chain on dense
polynomials with many real roots and repeated factors; `count_complex_60`
is the unchanged rectangle path. The only slower case is `count_cheb120`
(first call 0.007 s -> 0.009 s, repeated 0.007 s -> 0.008 s: the primitive
pseudo-remainder sequence pays content gcds on a degree-120 normal
sequence); it stays about 8x faster than SymPy. `factor_expr_300` already
returned SymPy's result on the base (observed; only `factor_list` of the
same polynomial raised).
