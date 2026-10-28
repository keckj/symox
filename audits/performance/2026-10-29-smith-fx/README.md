# Smith normal form over QQ, GF(p) and QQ[x], before and after

- Harness: `bench_case.py <case>` (one fresh process, one case: the first call,
  then the median of 7 timed batches of repeated calls) and
  `run.sh <base> <new> <out> [n]` (n=5 fresh processes per case and side, base
  and new interleaved). `summarize.py raw.jsonl` gives the medians in
  `summary.md` and checks that the printed results are identical on both sides.
- Base: tree of commit 58e50206 (wave 18, native Smith over ZZ only) with its
  release `_core.abi3.so`. New: this branch with its release build, same
  machine, same Python 3.10 venv. A single-job `cargo test` was compiling
  during the run (8 cores).
- Raw numbers: `raw.jsonl` (130 lines); medians: `summary.md`.

Result (medians, ms for the first call, us for repeated calls):

| case | first base | first new | ratio | repeat base | repeat new | ratio | same result |
|---|---|---|---|---|---|---|---|
| zz_snf_doc_3x3 | 1.907 | 1.284 | 1.49x | 50.3 | 51.4 | 0.98x | True |
| zz_invf_doc_3x3 | 1.824 | 1.299 | 1.40x | 44.0 | 45.3 | 0.97x | True |
| zz_decomp_doc_3x3 | 1.900 | 1.404 | 1.35x | 75.0 | 77.8 | 0.96x | True |
| qq_snf_3x3 | 2.461 | 1.304 | 1.89x | 338.2 | 63.1 | 5.36x | True |
| qq_decomp_6x6 | 7.549 | 1.768 | 4.27x | 4944.3 | 384.7 | 12.85x | True |
| qq_invf_10x10 | 11.580 | 1.710 | 6.77x | 9359.0 | 470.3 | 19.90x | True |
| gf7_invf_10x10 | 9.972 | 1.529 | 6.52x | 7915.7 | 269.7 | 29.35x | True |
| gf101_decomp_8x8 | 25.738 | 2.268 | 11.35x | 23942.6 | 733.8 | 32.63x | True |
| gf65537_snf_20x20 | 105.335 | 4.085 | 25.79x | 101395.8 | 2712.3 | 37.38x | True |
| qqx_invf_jordan_4x4 | 2.945 | 1.863 | 1.58x | 681.3 | 121.7 | 5.60x | True |
| qqx_decomp_jordan_4x4 | 5.794 | 1.733 | 3.34x | 3473.0 | 437.4 | 7.94x | True |
| qqx_snf_char_6x6 | 6.778 | 2.329 | 2.91x | 4816.6 | 1112.1 | 4.33x | True |
| qqx_decomp_3x4 | 5.221 | 1.810 | 2.88x | 2958.7 | 483.1 | 6.12x | True |

The `zz_*` first calls include the one-time capture of the admission snapshot
(717 entries, 129 functions, unchanged). Its capture now looks each name up
once per class namespace instead of on every pass, which takes it from about
2.25 ms to 1.6 ms in isolation and the ZZ 3x3 first call from 1.9 ms to 1.3 ms.
Repeated ZZ calls are unchanged within noise (0.96x to 0.98x, +1 to +3 us).
QQ, GF(p) and QQ[x] calls were Python bodies in the base and are now native:
4.3x to 37x faster on repeated calls, 1.6x to 26x on the first call.
