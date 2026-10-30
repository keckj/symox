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
