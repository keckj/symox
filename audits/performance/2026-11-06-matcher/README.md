# Native reference matcher (wave 43, 2026-11-06)

Differential evidence for `crates/symox-core/src/matches.rs`, the Rust port
of SymPy 1.14's `Expr.match` protocol (`Basic.matches`,
`AssocOp._matches_commutative`, `Mul`/`Add._combine_inverse`,
`Pow.matches`, `Wild.matches`, `Mul.as_base_exp`).

## Files

- `corpus.json.gz`: 2311 unique `expr.match(pattern, old)` calls.
  - `ode`, 1563 calls: every `.match(` made under `sympy/solvers/ode/` by
    `classify_ode` and by `dsolve(hint)` on the 315 cases of
    `audits/performance/2026-11-05-ode`.
  - `test_match`, 248 calls: from SymPy's `core/tests/test_match.py`.
  - `random`, 500 calls: seeded variations, seed 20261106.
  - The cases were harvested with `scripts/harvest.py`.
- `sympy-results.json.gz`: SymPy 1.14's answer for each case.
- `symox-before.json.gz` and `symox-after.json.gz`: symox replays without
  SymPy imports (`scripts/compare.py symox`).
- `report-before.json` and `report-after.json`: per-implementation counts,
  per-source counts, differing ids and timings (`scripts/compare.py report`).

## Implementations compared

- `python_port`: `symox._matching.matches`. This is the Python port, and it
  was the public `Expr.match` before this lot.
- `pattern`: the native backtracking `pattern::pattern_match`. The native
  ODE matchers (`second.rs`, `liouville.rs`, `reducible.rs`) called it
  before this lot.
- `native`: `_core._reference_match`, the new port. Its verdict is either a
  match or *declined*.
- `public`: `Expr.match` after this lot. It uses `native` and falls back to
  `python_port` when `native` declines.

## Results against SymPy

symox prints its `WildFunction` (a `Wild` subclass) as `Wild('g')`. The
comparison normalises that srepr, because the binding is the same.

| implementation | equal | different | declined |
|---|---|---|---|
| python_port (before = after) | 2309 | 2 | – |
| pattern (before = after) | 2059 | 252 | – |
| native (after) | 2226 | **0** | 85 |
| public (before) | 2309 | 2 | – |
| public (after) | **2311** | **0** | – |

Where the answers are equal, the binding insertion order is also equal for
`native` and `public`.

**`pattern` differences by source.** There are 99 ODE, 91 test_match and 62
random. None of the ODE differences comes from a pattern of `second.rs`
(the `a3, b3, c3` patterns, 360 calls), `liouville.rs` (`d, e`, `d, e, k`
and `fy`, 417 calls), `reducible.rs` (`a, b, c, d`, 123 calls) or the Airy
pattern (`a4, b4`, 14 calls). They come from these patterns:

- Bernoulli (`P, Q, n`): 37 calls.
- `eq_high_order_free` / `classify_ode` (`a, c1`): 38 calls.
- Bessel (`a4, c4`, `d4`, `a4, b, b4`, `a, b`): 23 calls.
- `P, Q`: 1 call.

symox did not match those calls natively.

**`python_port` differences.** The port gives `-sqrt(y**2)*cos(4)` for
`(y**2*cos(4)**2).match(d**2)`, where SymPy gives `y*cos(4)`. The second
case is the same with `z**2*cos(y)**2`. The cause is that symox's
`as_base_exp` lacks `Mul.as_base_exp`. The native port implements
`Mul.as_base_exp` and agrees with SymPy on both cases.

**Native declines.** The public call then falls back to the Python port.

| reason | count | where |
|---|---|---|
| noncommutative | 61 | test_match |
| `zoo` / `oo` | 18 | 17 ODE (Bernoulli on `zoo` coefficients), 1 test_match |
| bound `Sum` | 3 | test_match |
| `Or` / relationals | 1 | test_match |
| Float in a multiplicative inverse | 2 | ODE |

## Timing

The runs were single-process. Other agents' builds were running at the
same time, so read the numbers as orders of magnitude.

Corpus totals:

| | SymPy | python_port | pattern | native | public after |
|---|---|---|---|---|---|
| seconds | 25.9 | 14.8 | 0.12 | 1.93 | 1.98 |

On the ODE subset alone, `pattern` takes 0.09 s, `native` 1.49 s and
`python_port` 11.3 s. The native ODE matchers that used `pattern_match`
are therefore slower per call (about 1 ms instead of 0.06 ms). In exchange
they now follow the reference protocol. The Bessel matcher no longer calls
back into Python.
