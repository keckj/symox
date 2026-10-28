# Wave 19 polycap2: modular resultant, cyclotomic product formula, decomposition degree

Trees: before = commit 58e50206 with its release `_core.abi3.so`; after =
branch feat-polycap2-20261029 (same tree plus this change), release build.
Harness: `run.sh <before-python> <after-python> raw.jsonl` runs every case
of `bench.py` in 3 fresh processes per tree (first call, then the median of
5 repeated calls); `summarize.py raw.jsonl` gives the table below (median
over the 3 processes). Fingerprints are SHA-256 of `repr(result)`; a refusal
is reported as such and never counted as a speed-up. Runs used caps
`MemoryMax=6G`, `timeout 120`, on a machine shared with other jobs (the
integration/ODE corpora and cargo were running at the same time; before and
after are interleaved per case, so ratios compare like with like).

| case | before first (s) | after first (s) | ratio | before repeat (s) | after repeat (s) | ratio | result |
|---|---|---|---|---|---|---|---|
| res_witness_200_150 | 15.36 | 0.2238 | - | - | 0.2166 | - | before refused (ValueError: pseudo-division resource limit: PRS storage); after ok |
| res_60_45 | 0.08558 | 0.00458 | 18.68x | 0.08517 | 0.003988 | 21.36x | identical |
| res_dense_40_30 | 0.004086 | 0.001346 | 3.04x | 0.00356 | 0.0008867 | 4.01x | identical |
| res_small_zz | 0.0004972 | 0.0004074 | 1.22x | 0.0001816 | 0.0001521 | 1.19x | identical |
| res_small_qq | 0.000503 | 0.0004583 | 1.10x | 0.0001768 | 0.0001624 | 1.09x | identical |
| res_prs_small | 0.0004798 | 0.0004463 | 1.08x | 0.0001681 | 0.0001645 | 1.02x | identical |
| res_bivariate | 0.0006477 | 0.0005863 | 1.10x | 0.0002731 | 0.0002918 | 0.94x | identical |
| poly_res_25_20 | 0.003842 | 0.001133 | 3.39x | 0.00344 | 0.0007332 | 4.69x | identical |
| disc_dense_25 | 0.001438 | 0.0008217 | 1.75x | 0.001052 | 0.0004951 | 2.12x | identical |
| disc_small | 0.0004829 | 0.0003796 | 1.27x | 0.0001703 | 0.0001493 | 1.14x | identical |
| cyclotomic_20001 | 1.406e-05 | 0.06309 | - | - | 0.06401 | - | before refused (ValueError: polynomial structure resource limit: cyclotomic order); after ok |
| cyclotomic_15015 | 0.824 | 0.06237 | 13.21x | 0.8224 | 0.06278 | 13.10x | identical |
| cyclotomic_19999 | 0.2037 | 0.08047 | 2.53x | 0.2059 | 0.08059 | 2.55x | identical |
| factor_x1155_m1 | 0.08885 | 0.07471 | 1.19x | 0.08534 | 0.07471 | 1.14x | identical |
| decompogen_4000 | 0.06439 | 0.06442 | 1.00x | 0.06402 | 0.06503 | 0.98x | identical |
| decompogen_4100 | 0.0004154 | 0.06696 | - | - | 0.06636 | - | before refused (NotImplementedError: decompogen: functional decomposition is limited to degree 4096); after ok |

Slower cases: `res_bivariate` repeated calls 0.94x (the PRS path is
unchanged; the admission check unifies the two operands once more before
declining) and `decompogen_4000` 0.98x (within noise; same algorithm).

SymPy 1.14 reference on the three witnesses (same caps, run from /tmp):
`resultant(expand((x + 1)**200), expand((x - 2)**150 + 1), x)` 15.8 s
(bit length 47549, residue mod 2^61-1 equal to symox's);
`cyclotomic_poly(20001, x)` 1.8-2.7 s (degree 12992, height 2, same
coefficients as symox by fingerprint); `decompogen(x**4100 + x**2, x)`
0.38 s, `[x**2050 + x, x**2]` in both. `sympy_reference.py` is the probe
that produced the SymPy values hard-coded in `tests/test_native_polycap2.py`.

Remaining gap: `resultant(..., includePRS=True)` on the 200/150 witness still
reaches the PRS storage limit (the PRS is returned, so the modular path does
not apply).
