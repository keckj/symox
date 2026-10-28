# Smith normal form / number-field HNF port, before and after

- Harness: `bench_case.py <case>` (one fresh process, one case: first call, then
  the median of 7 timed batches of repeated calls) and `run.sh <base> <new> <out> [n]`
  (n=5 fresh processes per case and side, base and new interleaved).
  `summarize.py raw.jsonl` gives the medians in `summary.md` and checks that the
  printed results are identical on both sides.
- Base: tree of commit 34cb18b with its release `_core.abi3.so`. New: this
  branch with its release build, same machine, same Python 3.10 venv.
- Raw numbers: `raw.jsonl` (130 lines); medians: `summary.md`.

Result (medians, ms for first call, us for repeated calls):

| case | first base | first new | ratio | repeat base | repeat new | ratio |
|---|---|---|---|---|---|---|
| snf_doc_3x3 | 1.062 | 2.315 | 0.46x | 348.2 | 50.4 | 6.91x |
| invf_doc_3x3 | 1.115 | 1.893 | 0.59x | 337.7 | 44.9 | 7.52x |
| decomp_doc_3x3 | 1.984 | 2.403 | 0.83x | 871.4 | 75.9 | 11.47x |
| decomp_6x6 | 6.145 | 2.080 | 2.95x | 5153.5 | 133.7 | 38.54x |
| invf_8x8 | 6.151 | 2.020 | 3.05x | 4972.3 | 88.7 | 56.09x |
| decomp_4x6_big (|a| < 1e20) | 4.248 | 2.686 | 1.58x | 3233.2 | 282.3 | 11.45x |
| invf_12x12 | 16.776 | 2.042 | 8.21x | 15523.6 | 144.5 | 107.41x |
| matrix_snf_4x4 (Matrix API) | 3.151 | 3.280 | 0.96x | 1536.7 | 791.5 | 1.94x |
| round_two_x3m2 | 1.797 | 1.302 | 1.38x | 579.9 | 579.8 | 1.00x |
| round_two_cyclo7 | 2.453 | 1.838 | 1.33x | 659.9 | 638.9 | 1.03x |
| prime_decomp_5 | 3.946 | 3.318 | 1.19x | 2268.7 | 2282.4 | 0.99x |
| nf_hnf_5x8 | 0.212 | 0.158 | 1.34x | 159.4 | 85.5 | 1.86x |
| nf_hnf_mod_4x4 | 0.683 | 0.171 | 3.99x | 107.3 | 72.4 | 1.48x |

Regression: the first Smith call of a process pays the one-time capture of the
admission snapshot (about 1 ms), so the first call on 3x3 inputs is 0.46x to
0.83x of the base (+0.4 to +1.3 ms). Repeated calls and larger inputs are
faster. round_two/prime_decomp are dominated by other Python code: repeat
times are unchanged within noise (0.99x to 1.03x).
