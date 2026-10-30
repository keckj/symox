# markovrust (wave 33): native finite Markov chain computations

Harness: `bench.py` (one public call per process: `first` is the first call,
`median` the median of 5 repeated calls, plus the srepr of the result for the
cross-check). Release builds, every run under
`systemd-run --user --scope -q -p MemoryMax=6G -p MemorySwapMax=0`, SymPy
runs from /tmp.

    python bench.py sympy|symox CALL 5

Before is the base tree (363879cc), after is this branch. Raw lines:
`raw_sympy.jsonl` (SymPy 1.14.0), `raw_before.jsonl`, `raw_after.jsonl`
(seconds). The machine was shared with other builds (load average 4-5), so
differences under about 15% are noise.

Chains: `X` 8 states, dense rational; `Y` 8 states, two absorbing states and
transient states (rational, density 0.5); `Z` 6 states, rational; `C` a 5x5
rational generator.

## Repeated calls (median)

| call | before | after | SymPy 1.14 | after/before | after/SymPy |
|---|---|---|---|---|---|
| `classes8` (`Y.communication_classes()`) | 341.2 us | 60.1 us | 744.0 us | 0.18 | 0.08 |
| `decompose8` | 404.8 us | 55.2 us | 1.21 ms | 0.14 | 0.05 |
| `canonical8` | 524.1 us | 56.6 us | 4.15 ms | 0.11 | 0.01 |
| `stationary8` (`X.stationary_distribution()`) | 443.0 us | 225.7 us | 14.34 ms | 0.51 | 0.02 |
| `fundamental8` | AttributeError | 41.2 us | 4.04 ms | - | 0.01 |
| `absorbing8` | AttributeError | 50.2 us | 5.30 ms | - | 0.01 |
| `regular8` (`X.is_regular()`) | AttributeError | 72.7 us | 914.1 us | - | 0.08 |
| `ctmc_limiting5` | 1.15 ms | 46.6 us | 1.39 ms | 0.04 | 0.03 |
| `P_path6x4` (4 equalities, 6^4 paths) | 143.98 ms | 4.75 ms | 47.23 ms | 0.03 | 0.10 |
| `E_path8x3` (`E(X[1]+X[2]+X[3], Eq(X[0],0))`) | 56.10 ms | 1.91 ms | 2.81 ms (refusal) | 0.03 | - |
| `P_path8x4` (4 `Ne`, 8^4 paths) | 464.51 ms | 15.55 ms | 137.45 ms (refusal) | 0.03 | - |

## First call

| call | before | after | SymPy 1.14 |
|---|---|---|---|
| `classes8` | 526.2 us | 274.9 us | 1.14 ms |
| `decompose8` | 805.7 us | 338.4 us | 1.80 ms |
| `canonical8` | 885.5 us | 330.1 us | 4.89 ms |
| `stationary8` | 1.03 ms | 642.9 us | 40.66 ms |
| `fundamental8` | AttributeError | 573.7 us | 5.33 ms |
| `absorbing8` | AttributeError | 380.3 us | 6.56 ms |
| `regular8` | AttributeError | 171.5 us | 1.36 ms |
| `ctmc_limiting5` | 1.79 ms | 307.9 us | 24.57 ms |
| `P_path6x4` | 153.82 ms | 5.52 ms | 166.27 ms |
| `E_path8x3` | 67.86 ms | 2.39 ms | 25.65 ms (refusal) |
| `P_path8x4` | 463.80 ms | 16.32 ms | 265.90 ms (refusal) |

## Cross-check

- After returns the same srepr as before on every call that existed
  before, and the same srepr as SymPy 1.14 on the 9 calls SymPy answers
  (`fundamental8`, `absorbing8`, `regular8` were missing before).
- `E_path8x3`: SymPy 1.14 raises NotImplementedError ("Mechanism for
  handling ... queries hasn't been implemented yet."). `P_path8x4`: SymPy
  raises ValueError ("The process is in multiple states at X[3] ...").
  symox answered both before this wave (the exact path sum), and the native
  enumeration gives the same value. The SymPy column is the time of its
  refusal, so no ratio is given.
- No slowdown was measured.

## Round 2: Float and symbolic systems, query assembly

`bench_r2.py sympy|symox` (same machine, mean over at least 0.5 s, from
/tmp under a 6G memory cap):

| call | symox | SymPy 1.14 |
|---|---|---|
| Float 6x6 `stationary_distribution` | 0.446 ms | 19.321 ms |
| Float 6x6 `fundamental_matrix` | ValueError (see parity gap: shared Float `det`) | 6.806 ms |
| Float 4x4 `fundamental_matrix` (binary-exact rows) | 0.593 ms | 12.838 ms |
| symbolic 3x3 `stationary_distribution` | 0.279 ms | 7.683 ms |
| CTMC symbolic 3x3 `limiting_distribution` | 0.188 ms | 2.477 ms |
| CTMC Float 4x4 `limiting_distribution` | IndexError | IndexError |
| CTMC `P(Eq(C(1), 1), Eq(C(0), 1))` | 1.929 ms | 16.672 ms |
| CTMC `variance(C(3/2), Eq(C(0), 1))` | 3.597 ms | 41.799 ms |

Before round 2 these went through the Python `linsolve` (the Float and
symbolic stationary vectors, with different last bits and uncancelled
forms) and `simplify` (the queries), so the values differed from SymPy's;
no before/after ratio is given. Values: the 260 round-2 probe results
(tests/test_native_markovrust_r2.py) are SymPy's srepr except the known
`Indexed` parameter form of a non-unique CTMC limiting distribution.
