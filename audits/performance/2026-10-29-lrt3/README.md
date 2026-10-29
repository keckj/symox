# Native log_to_real for auxiliary polynomials of degree >= 3, before and after

- Harness: `bench_case.py <case>` (one fresh process, one public
  `ratint(f, x)` case: the first call, then the median of 7 timed batches of
  repeated calls) and `run.sh <base> <new> <out> [n]` (n=5 fresh processes
  per case and side, base and new interleaved, each under a 6G memory cap).
  `summarize.py raw.jsonl` gives the medians in `summary.md` and checks that
  the printed results (srepr, Dummy indices removed) are identical on both
  sides.
- Sweep: `diff.py <out.json>` runs 94 ratint inputs once each (binomials
  `x**n + a`, n = 3..8, and mixed quartic/sextic denominators);
  `cmp.py sweep_base.json sweep_new.json` compares the srepr outputs with
  Dummy indices normalized and prints the cases whose single-call time
  changed by more than 1.3x. Result: 94 identical, 0 different.
- Base: tree of commit 9e2340a1 (wave 19 integration) with its release
  `_core.abi3.so`. New: this branch, release build of the same tree plus
  this change, same machine, same Python 3.10 venv. Other agents' builds and
  this branch's single-job `cargo test` and corpora were running (8 cores),
  so absolute numbers are noisy; base and new are interleaved.
- Raw numbers: `raw.jsonl` (150 lines), `sweep_base.json`, `sweep_new.json`;
  medians: `summary.md`.

Result (medians, ms for the first call, us for repeated calls):

| case | integrand | first base | first new | ratio | repeat base | repeat new | ratio | same |
|---|---|---|---|---|---|---|---|---|
| lrt_x4p1 | 1/(x**4+1) | 6.561 | 3.876 | 1.69x | 5869.2 | 3301.2 | 1.78x | True |
| lrt_x4p4 | 1/(x**4+4) | 5.596 | 2.948 | 1.90x | 4854.5 | 2517.6 | 1.93x | True |
| lrt_x4pq | 1/(x**4+1/4) | 5.493 | 2.938 | 1.87x | 4804.5 | 2442.5 | 1.97x | True |
| lrt_x4mx2p1 | 1/(x**4-x**2+1) | 7.379 | 4.143 | 1.78x | 6645.5 | 3510.3 | 1.89x | True |
| lrt_mixed | 1/((x**2+1)*(x**2+4)*(x-1)) | 10.321 | 5.352 | 1.93x | 9754.3 | 4896.1 | 1.99x | True |
| lrt_x4_x6p1 | x**4/(x**6+1) | 9.526 | 7.642 | 1.25x | 8428.9 | 6982.3 | 1.21x | True |
| lrt_x6p1 | 1/(x**6+1) | 9.452 | 7.548 | 1.25x | 8346.8 | 7001.8 | 1.19x | True |
| lrt_x3p1_x_1 | 1/(x**3+x+1) | 177.124 | 172.611 | 1.03x | 162498.9 | 158135.7 | 1.03x | True |
| lrt_x3m2 | 1/(x**3-2) | 21.282 | 21.725 | 0.98x | 4864.3 | 4796.6 | 1.01x | True |
| lrt_x5p1 | 1/(x**5+1) | 31.335 | 31.103 | 1.01x | 13828.2 | 13677.8 | 1.01x | True |
| lrt_x6m2 | 1/(x**6-2) | 30.644 | 30.937 | 0.99x | 13208.2 | 12781.4 | 1.03x | True |
| lrt_x8p1 | 1/(x**8+1) | 28.394 | 31.354 | 0.91x | 11899.2 | 11991.7 | 0.99x | True |
| lrt_x2_x4px1 | x**2/(x**4+x+1) | 20.861 | 20.561 | 1.01x | 19474.0 | 19322.5 | 1.01x | True |
| lrt_x5p2 | 1/(x**5+2) | 6.215 | 6.054 | 1.03x | 5147.5 | 5100.9 | 1.01x | True |
| lrt_x7p1 | 1/(x**7+1) | 11.068 | 11.148 | 0.99x | 10064.7 | 10144.9 | 0.99x | True |


The first seven cases are computed natively now (1.2x to 2.0x). The other
eight decline to the unchanged Python algorithm by design (EX `c(r_u, v)`,
cube-root or Cardano resultant roots, incomplete real root lists where the
exact RootSum is kept). A declined `(h, q)` is recorded in a small
per-thread memo (8 entries, declines only, never a signal), so the public
`log_to_real` after a declined native ratint driver does not try again, and
an odd-degree `q` whose factors over QQ have a real root of degree >= 2 or
are not symmetric about the mean of their roots declines after one
univariate factorization, before the resultant (exact: such a `q` gives
`R(u)` an irrational real root, which the pre-check always declines for odd
degree). The previous version (commit 5574a340) attempted twice: repeat
ratios 0.71x for 1/(x**8+1) and 0.62x for 1/(x**7+1). Now all declined cases
are within noise (0.99x to 1.03x repeated; repeated calls also hit the memo).
The first call of an even-degree declined case still pays one native
attempt (resultant and factorization of `R`, about 2-3.5 ms): 1/(x**8+1) is
0.91x on the first call. In the sweep (single calls, `diff.py`, file
`sweep_new.json`), the only case below 0.95x is `1/(x**6+8)` (0.73x, was
0.59x); the `x**7 + a` and `x**8 + a` families, which were 0.57x-0.76x,
are no longer below 0.95x.
`1/(x**8-1)` is 8.21 s base and 8.41 s new: an earlier version of this
change tripled it (24 s) by computing the 8 s `roots` of a sextic factor of
`c(1/8, v)` natively before declining; the degree >= 5 factor check now
declines first.
