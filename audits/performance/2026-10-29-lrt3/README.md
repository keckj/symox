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
| lrt_x4p1 | 1/(x**4+1) | 6.779 | 3.972 | 1.71x | 6229.1 | 3440.1 | 1.81x | True |
| lrt_x4p4 | 1/(x**4+4) | 5.730 | 3.139 | 1.83x | 5071.6 | 2640.8 | 1.92x | True |
| lrt_x4pq | 1/(x**4+1/4) | 6.024 | 3.023 | 1.99x | 4827.1 | 2525.2 | 1.91x | True |
| lrt_x4mx2p1 | 1/(x**4-x**2+1) | 7.725 | 4.211 | 1.83x | 6983.8 | 3633.7 | 1.92x | True |
| lrt_mixed | 1/((x**2+1)*(x**2+4)*(x-1)) | 10.682 | 5.388 | 1.98x | 10086.1 | 5017.9 | 2.01x | True |
| lrt_x4_x6p1 | x**4/(x**6+1) | 9.708 | 7.831 | 1.24x | 8649.8 | 7321.7 | 1.18x | True |
| lrt_x6p1 | 1/(x**6+1) | 9.699 | 7.697 | 1.26x | 8411.8 | 7142.6 | 1.18x | True |
| lrt_x3p1_x_1 | 1/(x**3+x+1) | 184.164 | 181.222 | 1.02x | 163618.6 | 165022.3 | 0.99x | True |
| lrt_x3m2 | 1/(x**3-2) | 21.901 | 23.525 | 0.93x | 5012.8 | 5243.9 | 0.96x | True |
| lrt_x5p1 | 1/(x**5+1) | 32.696 | 33.721 | 0.97x | 14344.4 | 15223.4 | 0.94x | True |
| lrt_x6m2 | 1/(x**6-2) | 34.033 | 32.638 | 1.04x | 13751.7 | 14303.7 | 0.96x | True |
| lrt_x8p1 | 1/(x**8+1) | 30.754 | 35.730 | 0.86x | 12427.5 | 17539.0 | 0.71x | True |
| lrt_x2_x4px1 | x**2/(x**4+x+1) | 21.167 | 22.352 | 0.95x | 20405.8 | 20413.3 | 1.00x | True |
| lrt_x5p2 | 1/(x**5+2) | 6.368 | 6.785 | 0.94x | 5324.6 | 5773.0 | 0.92x | True |
| lrt_x7p1 | 1/(x**7+1) | 11.321 | 17.681 | 0.64x | 10332.6 | 16628.7 | 0.62x | True |

The first seven cases are computed natively now (1.2x to 2.0x). The other
eight decline to the unchanged Python algorithm by design (EX `c(r_u, v)`,
cube-root or Cardano resultant roots, incomplete real root lists where the
exact RootSum is kept). They pay the native attempt (resultant,
factorization of `R` and Sturm counts), twice when the native ratint driver
declines and the public `log_to_real` tries again: within noise for most,
but 1/(x**8+1) is 0.71x (+5 ms) and 1/(x**7+1) 0.62x (+6 ms). In the sweep,
the median single-call ratio is 0.98x (range 0.57x to 2.19x); the
`x**7 + a`, `x**8 + a` and `x**6 + 8` families lose 5 to 8 ms each.
`1/(x**8-1)` is 8.21 s base and 8.41 s new: an earlier version of this
change tripled it (24 s) by computing the 8 s `roots` of a sextic factor of
`c(1/8, v)` natively before declining; the degree >= 5 factor check now
declines first.
