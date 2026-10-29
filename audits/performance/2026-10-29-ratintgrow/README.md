# ratint growth family and the non-real field decline, before and after

- Harness: `bench_case.py <case>` (one fresh process, one public
  `ratint(f, x)` case: the first call, then the median of 7 timed batches of
  repeated calls) and `run.sh <base> <new> <out> [n]` (n=3 fresh processes
  per case and side, base and new interleaved, each under a 6G memory cap).
  `summarize.py raw.jsonl` gives the medians in `summary.md` and checks that
  the printed results (srepr, Dummy indices removed) are identical on both
  sides.
- Base: tree of commit 65a3f2af (wave 22 integration) with its release
  `_core.abi3.so`. New: this branch, release build of the same tree plus the
  structural decline in `integrate/log_to_real.rs`. Same machine, same
  Python 3.10 venv; other agents' builds and this branch's suite/corpora were
  running, so absolute numbers are noisy (+-5%).
- Raw numbers: `raw.jsonl`; medians: `summary.md`.

Result (medians, ms for the first call, us for repeated calls):

| case | integrand | first base | first new | ratio | repeat base | repeat new | ratio | same |
|---|---|---|---|---|---|---|---|---|
| x6p8 | 1/(x**6+8) | 15.031 | 11.588 | 1.30x | 9974.6 | 9745.8 | 1.02x | True |
| x4m2x2p4 | 1/(x**4-2*x**2+4) | 9.802 | 8.463 | 1.16x | 7594.9 | 7593.1 | 1.00x | True |
| x4p2 | 1/(x**4+2) | 25.878 | 26.664 | 0.97x | 7316.9 | 7346.1 | 1.00x | True |
| x8p1 | 1/(x**8+1) | 32.295 | 28.901 | 1.12x | 11008.6 | 10911.3 | 1.01x | True |
| x8p16 | 1/(x**8+16) | 32.994 | 29.217 | 1.13x | 11670.4 | 11823.7 | 0.99x | True |
| x4pxp1 | 1/(x**4+x+1) | 19.519 | 19.852 | 0.98x | 18368.4 | 18275.4 | 1.01x | True |
| x2_x6p8 | x**2/(x**6+8) | 1.603 | 1.563 | 1.03x | 1182.5 | 1172.5 | 1.01x | True |
| x4p1 | 1/(x**4+1) | 3.706 | 3.842 | 0.96x | 3113.3 | 3209.9 | 0.97x | True |
| x4mx2p1 | 1/(x**4-x**2+1) | 4.049 | 4.129 | 0.98x | 3482.5 | 3522.0 | 0.99x | True |
| xm1_4p1 | 1/((x-1)**4+1) | 3.629 | 3.787 | 0.96x | 3097.8 | 3259.6 | 0.95x | True |
| x6p27 | 1/(x**6+27) | 7.841 | 7.884 | 0.99x | 7310.7 | 7394.2 | 0.99x | True |
| x6p1 | 1/(x**6+1) | 7.060 | 7.018 | 1.01x | 6456.3 | 6494.3 | 0.99x | True |
| x3pxp1 | 1/(x**3+x+1) | 173.756 | 169.970 | 1.02x | 157420.1 | 152373.2 | 1.03x | True |
| x3m3xp1 | 1/(x**3-3*x+1) | 7.659 | 7.734 | 0.99x | 6623.9 | 6566.4 | 1.01x | True |
| x5p1 | 1/(x**5+1) | 32.181 | 32.152 | 1.00x | 13938.7 | 13872.2 | 1.00x | True |

Reading:

- Only first calls of the declined groups (`1/(x**6+8)`,
  `1/(x**4-2*x**2+4)`, `1/(x**8+1)`, `1/(x**8+16)`) gain (1.12x-1.30x): the
  native kernel used to compute the resultant `R(u)`, isolate its roots and
  only then decline on the biquadratic `Q(sqrt(2), sqrt(6))` pairs. Repeated
  calls were already served by the thread-local decline memo, so they are
  unchanged. `1/(x**4+2)` (root field `Q(2**(1/4), i)`, degree 8) declined
  early before as well; its 0.97x is within noise. The kept native cases
  (`x4p1`, `x4mx2p1`, `xm1_4p1`, `x6p27`, `x6p1`) pay the factorization and
  the cubic resolvent check: 0.95x-1.0x, mostly noise.
- All 15 results are identical between base and new.

SymPy 1.14 reference (same machine, first in-process call of `ratint`, ms):
`1/(x**6+8)` 131.7, `1/(x**3-3*x+1)` 53.6, `1/(x**4+x+1)` 157.2,
`1/(x**4-4*x**2+x+1)` 165.7, `x/(x**3-3*x+1)` 50.3, `1/(x**3-4*x+1)` 75.3,
`1/((x-1)**4+1)` 76.8, `1/(x**4-2*x**2+4)` 80.7, `1/(x**6+27)` 110.5,
`1/(x**4+2)` 218.8, `x**2/(x**6+8)` 27.4, `1/(x**8+16)` 3396.4 (returns 0,
see the intentional correction in docs/algorithmic-parity.md),
`1/(x**4+x**2+1)` 88.4, `1/(x**3+x+1)` about 7200 (`ratint` and
`integrate`). Symox is faster than SymPy on all of them: the "0.73x" quoted
for `1/(x**6+8)` in TODO.md was a symox base-vs-new ratio of an earlier
wave (wasted native attempt on the first call), not a ratio against SymPy.

`1/(x**3+x+1)`: the base already prints 7892 characters (SymPy 16653;
the 82704 of the 2026-10-26 audit predates the ratint-dom EX cancel), and
`integrate` returns the same expression as `ratint`. Verification: the
derivative of the primitive minus the integrand at `x = 3/7` is below
`1e-40` at 50 digits, the definite value `F(2) - F(1)` matches SymPy's form
to 38 digits, and the quadratic log argument is `x**2 + b*x + c` with `b`
equal to the real root of `x**3 + x + 1` (structural check). An exact
symbolic certificate (minimal polynomial of the derivative difference) did
not finish: SymPy's `minimal_polynomial` timed out after 900 s and symox's
ran out of the 6G cap, so no exact proof of the cubic form is claimed here.
