# limrest3: limits held after limrest2, before and after

- Harness: `bench_case.py <case>` (one fresh process, one public call: the
  first call, then the median of 7 timed batches of repeated calls) and
  `run.sh <base_pythonpath> <new_pythonpath> <out.jsonl> [n]` (n = 5 fresh
  processes per case and side, base and new interleaved, each under a 6G
  memory cap). `summarize.py raw.jsonl` gives the medians in `summary.md`
  with the printed results of both sides.
- Base: commit 82d674aa (wave 20 integration) with its release
  `_core.abi3.so`. New: this branch, release build. Same Python tree (only
  the Rust extension differs), same machine and venv. Other agents' builds
  and this branch's single-job `cargo test` were running (8 cores), so
  absolute numbers are noisy.
- Raw numbers: `raw.jsonl` (170 lines); medians: `summary.md`.

Result: the nine first cases were not computed on the base (it returns the
unevaluated `Limit`/`Integral`), so their times are not speed-ups: they are
the cost of a computed answer against the cost of a decline. All new values
equal SymPy 1.14. The Stirling remainder decline cost 83 ms on the base (the
hidden zero sent Gruntz into failing retries); it is now computed in 17 ms.
The eight control cases (already computed on the base: `sin(x)/x`,
`(1 + 1/x)**x`, a Gruntz `exp*sin` difference, `Ei(x)*exp(-x)*x`,
`(1 - erf(x))*exp(x**2)*x`, `harmonic(x, 3)`, `loggamma(x)/(x*log(x))`,
`atan(x)*x - pi*x/2`) are unchanged within noise (0.95x to 1.05x) with
identical results.

See `summary.md` for the table.
