# ratintex: QQ<a> logpart, EX log_to_atan, RootSum order, exroots forms

- Harness: `diffcorpus.py sympy|symox OUT.jsonl` (38 cases: `ratint_logpart`
  over one or two algebraic fields, `log_to_atan` over EX, `integrate` of
  rational functions with radicals or Cardano roots, `ratint(..., real=False)`,
  `roots` of decomposable polynomials). Each case is run once for warm-up,
  then timed as the best of three runs. `compare.py ref.jsonl X.jsonl ...`
  compares printed results (same exception class for errors).
- `raw/sympy.jsonl`, `raw/before.jsonl`, `raw/after.jsonl`: first pass of the
  lot (35 cases; before = branch base, after = the WIP commit c436f030).
  `raw/sympy2.jsonl`, `raw/after2.jsonl`: final pass (38 cases, the three
  two-field logpart cases added), table in `raw/compare2.md`.
  `raw/merged.jsonl`: the WIP tree after the merge of main (2e623e66).
  `raw/after3.jsonl`: the final build (53-bit refinement of print-order
  near ties only), same 35/38.
- SymPy 1.14 in `~/symox-ci/venv`; symox: release build of this branch.
  Other agents' builds were running, so timings are noisy (+-10%).
- `raw/sympy_tests/`: SymPy's test_rationaltools.py (12/12) and
  test_polyroots.py (17/32, same 15 failures before and after; imports
  made optional by `soften.py`).

Result (final pass): 35/38 printed results identical to SymPy (before the
lot: 19/35; after the WIP commit: 30/35). New since the WIP commit:

| case | before | now |
|---|---|---|
| logpart QQ<s2>/QQ<s3> 1/(x**2+s3) | TypeError | SymPy's pair |
| logpart QQ<s2>/QQ<s3> (x+s2)/(x**2+s3*x+1) | TypeError | SymPy's pair |
| logpart QQ<s2>/QQ<2**(1/3)> x/(x**3+2**(1/3)) | TypeError | SymPy's pair |
| integrate 1/(x**3+x+1) | 9569 chars, other exact form | SymPy's (16653 chars) |
| integrate 1/(x**3+2*x-1) | 9948 chars, other exact form | SymPy's (17352 chars) |

The three remaining differences are the `nested_binomial` deviation
(awaiting user validation, docs/algorithmic-parity.md, ratintex):
`roots(x**8-2*x**4+x**2-3)`, the same octic in `x**2+1`, and the preimages
of the real Cardano root in `roots(x**6+x**2+1)`; their zero sets are
complete (checked to 1e-20 in tests/test_native_ratintex.py).

Times (ms, best of three; symox / SymPy): `integrate 1/(x**3+x+1)`
336 / 2190, `1/(x**3+2*x-1)` 312 / 2123, two-field logparts 1.7-8.5 /
46-105. Printing the 16653-character primitive: 0.022 s (SymPy 0.85 s);
only near ties of the print order are recomputed at 53 bits (evaluating
every factor at 53 bits took 0.16 s unmemoized, and memoized it pushed the
tests/extra shard past its 6 GiB cgroup).
