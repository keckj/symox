# Holonomic wave 22 (feat-holo-20261029d): before/after timings

- before: base commit c7289185 Python tree (git archive) with the same
  release `_core.abi3.so`; after: this branch. No Rust changed.
- `run.sh <before-python> <after-python> raw.jsonl`: per case, 5 fresh
  processes timing the first call ("first", median of 5) and one process
  timing the median of 7 repeated calls after a warm-up ("repeat").
- Raw numbers: `raw.jsonl`; table: `summary.txt` (milliseconds, ratio =
  after/before; > 1 is slower).

Slower cases and why:
- `erf_rewrite_erfc` (0.014 -> 0.149 ms repeat): before, no rule existed and
  `erf(x)` was returned unchanged; after, the SymPy rewrite `1 - erfc(x)` is
  built.
- `sin_pow3` (1.24x repeat): `__pow__` now uses SymPy's square-and-multiply
  (Dx seed at x0 = 0, then f, f**2 products), matching SymPy's results.
- `float_cos`, `float_param_sin` (1.3x-2x repeat): Float inputs now follow
  SymPy's RR domain (exact placeholder annihilator, then Float restore and
  marker-based normalization). `float_cos` raised ValueError before.
- `exp_pow_neg` raised NotImplementedError before; it now returns the SymPy
  result.
