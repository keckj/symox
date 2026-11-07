# hasmatch (wave 44): reference `has` and the `pattern_match` callers

Evidence for the lot `hasmatch` (branch `feat-hasmatch-20261106`, from
`main` 583315a8). Reference: SymPy 1.14.

## 1. `has`, `has_free` and `has_xfree`

The native `traverse::has` / `has_structural` (Rust) now follow
`Basic._has` and `AssocOp._has_matcher`. The rules:

- A pattern of class `Add`, `Mul`, `Max`, `Min`, `Union`, `Intersection`,
  `And` or `Or` is found in a node of the same class when:
  - its arguments with `is_commutative is True` are a subset of the node's
    arguments;
  - its other arguments are a consecutive run of the node's other
    arguments, in the reference argument order.
- Arguments with unknown commutativity take the ordered rule. These include
  sets other than `Union`/`Intersection`, relationals, `true`/`false`,
  `Tuple`, `Not`, `Xor` and the like.
- Python patterns that `_sympify` rejects are skipped, so `x.has('x')` and
  `x.has([x])` are False. This is an admission rule in the PyO3 binding,
  for both expressions and matrices.
- `has_free` uses the same Rust code through `_core._has_matcher`. It
  admits patterns as the reference does: a single non-Basic iterable
  raises TypeError, and so does an unhashable pattern (`set(patterns)`).
- The duplicate `matches::reference_has` is gone. The matcher now calls
  `has_structural`. The matcher corpus of wave 43
  (`../2026-11-06-matcher`, 2311 calls) gives identical results for all
  four implementations before and after this change.

### Corpus

`corpus.json.gz` holds 1255 cases, built by `scripts/has_corpus.py harvest`:

- hand-written sections: sub-products, sub-sums, noncommutative operands,
  powers, functions and classes, Wild, several arguments, derivatives,
  integrals, relationals/Piecewise/logic, matrices, mixed, `has_free`,
  `has_xfree`, and lattice/boolean operators;
- seeded random sub-collections of the arguments of sums and products,
  over commutative, free-symbol and noncommutative operands.

The SymPy answers are in `sympy-has.json`. The symox answers are in
`symox-has-before.json` (main 583315a8) and `symox-has-after.json`, and the
reports are in `report-has-*.txt`.

| | SAME | DIFF | DECLINED | WRONG |
|---|---|---|---|---|
| before | 1117 | 6 | 0 | 132 |
| after | **1255** | 0 | 0 | **0** |

The 132 wrong answers before the fix:

- commutative sub-products, including under functions, integrals and
  matrices (`(2*x*y*z).has(x*y)`, `(-x*y).has(x*y)`);
- `Max`/`Min`/`Union`/`And`/`Or` sub-collections;
- a noncommutative sub-sum (`(A + B + C).has(A + C)` was True);
- a `has_free` sub-product.

The 6 DIFF cases were `has('x')`, `has([x])` and the `has_free` admission
errors.

`tests/test_has_parity.py` replays the corpus against live SymPy, together
with 16 named reference cases.

### Internal Rust callers of `has`

The audit was dynamic. A temporary instrumented build compared, on every
call, the new answer with an emulation of the old one (for `Add`, `Mul`,
`Max`, `Min` and set patterns). It logged every caller location
(`#[track_caller]`) whose answer changed. The instrumentation was removed
before the commit. Runs and changed answers from internal callers:

| run | changed answers from internal callers |
|---|---|
| copies of the SymPy test files (`scripts/run_sympy_tests.sh`) | 0 |
| ODE corpus, 315 cases | 0 |
| full suite | 0 |

In the full suite, 129 answers changed. All of them came from the public
bindings: `crates/symox-py/src/lib.rs` (`has`: 123 now True, 3 now False)
and `crates/symox-py/src/matrices.rs` (`Matrix.has`: 3 now True).

So the change of `has` reaches public results only through the public
`has` / `has_free`. Two suite tests encoded the old deviations and were
updated to the SymPy answers: `expr.has((Derivative, Symbol))` is False
(it raised TypeError before), and `expr.has('x')` is False (it was True).

## 2. Callers of the old `pattern::pattern_match`

The audit was dynamic (temporary build with `SYMOX_PM_AUDIT`, removed
before the commit). Each call made by the three callers logged the old
answer and the answer of `matches::reference_match`, the native port of
`Expr.match`. Calls that the port declines were replayed in SymPy 1.14.
The scripts are `scripts/pm_audit.py` (classification) and
`scripts/pm_sound.py` (whether an old binding substituted into the pattern
reproduces the expression).

Distinct logged calls, after removing duplicates. "reference" means the
native port answered. "sympy" means the port declined and the call was
replayed in SymPy, with a 10 s limit per call.

| caller | run | SAME | OLD-NONE | OLD-OTHER-BINDING | OLD-SPURIOUS | not replayed |
|---|---|---|---|---|---|---|
| collect | SymPy tests | 25 | 0 | 0 | 0 | 0 |
| collect | ODE corpus | 52 | 0 | 0 | 0 | 0 |
| meijerint | SymPy tests | 1064 | 3 | 0 | 0 | 0 |
| simp | SymPy tests | 16 | 1144 | 0 | 0 | 19 (17 timeouts, 2 unparsed) |
| simp | ODE corpus | 89 | 5148 | 41 | 0 | 38 (timeouts) |

OLD-NONE means the old matcher found no match where the reference matched.
For `simp`, these calls fall into two groups.

- **Degenerate bindings: 1013 + 4508 calls.** The reference binds a wild
  coefficient to 0, so the rule rewrites the expression to itself.
- **Algebraic inversion by `_combine_inverse`: 131 + 640 calls.** The
  reference binds the coefficient to a quotient. An example is
  `cos(2*log(x))**2` against `2*a*cos(b)**2 - a + c`, with `a = 1/2`,
  `c = 1/2`.
  - On these calls, the reference matcher would make symox's heuristic
    rules fire where they do not fire today, for instance power
    reductions.
  - SymPy has no such rules, so this would move simplification outputs
    away from SymPy rather than towards it.

The 41 OLD-OTHER-BINDING bindings are all sound: each makes the rewrite an
identity (`pm_sound.py`). No old answer was spurious. Logs and JSON
classifications are kept out of the tree; `pm_audit.py` regenerates them
from an instrumented build.

### Decisions

- **`poly/collect.rs`, fixed in Rust.**
  - `match_factors` now calls `matches::reference_match`. It falls back to
    `pattern_match` only when the port declines (noncommutative, infinite
    or bound operands).
  - The logged calls agreed. However, probes with Wild patterns showed
    wrong public results:
    - `collect((x - 2)*y + z, -x + a, evaluate=False)` gave
      `{y: x - 2, z: 1}`. SymPy gives `{1: y*(x - 2) + z}`.
    - The `(-x + a)**2` and `a - x` variants were wrong in the same way.
  - `tests/test_collect_match_parity.py` has 10 cases against live SymPy:
    3 failed before the fix and 0 after.
- **meijerint `_condsimp`, fixed in Rust (`integrate/meijerint/conds.rs`).**
  - The port declines relationals as a whole. `cond_match` therefore pairs
    the two arguments of relationals of the same class, as `Basic.matches`
    does: each argument goes through `reference_match`, with the bindings
    of the previous one substituted.
  - The old matcher missed 3 of 1057 distinct calls. All three were the
    rule `And(p**2 > 1, Ne(p**2, 1)) -> p**2 > 1` against
    `re(s) + re(a) > 1` and `cos(Abs(arg(t)))*sqrt(Abs(t**2)) > 1`, where
    SymPy binds `p = sqrt(X)`.
  - With the fix, `And(Ne(X, 1), X > 1)` reduces to `X > 1`, as in SymPy.
    The Rust unit test is
    `condsimp_merges_square_rule_through_relational_matching`.
  - No SymPy test of `test_meijerint` / `test_transforms` /
    `test_integrals` reaches that merge (probed with an instrumented
    `_condsimp`), so the copies of those test files do not change.
- **`simp.rs`, `pattern_match` kept (deliberate).**
  - The patterns are symox's own heuristic trigonometric identities, such
    as `a*sin(b)**2 + a*cos(b)**2 + c -> a + c`, with no SymPy counterpart.
  - Where the old binding differs from the reference, two cases occur:
    - the reference takes a degenerate binding (`a = 0`);
    - the old matcher misses.
  - Every differing old binding was checked sound with SymPy
    (`pm_sound.py`). Each rewrite is therefore an identity whatever binding
    is chosen. Switching the matcher would change which identity fires,
    not correctness.

## 3. SymPy test files

`scripts/mkcopies.py` rewrites a set of SymPy 1.14 test files for symox:

- simplify: `test_simplify`, `test_trigsimp`, `test_fu`, `test_radsimp`;
- meijerint and transforms: `test_meijerint`, `test_transforms`,
  `test_laplace`;
- core: `test_match`, `test_basic`.

`scripts/run_sympy_tests.sh` runs them.

| file | before pass/fail | after pass/fail |
|---|---|---|
| test_core_basic | 11/14 | 12/13 (`test_has` now passes) |
| test_core_match | 41/2 | 41/2 |
| test_integrals_laplace | 2/2 | 2/2 |
| test_integrals_meijerint | 20/13 | 20/13 |
| test_integrals_transforms | 9/6 | 9/6 |
| test_simplify_fu | 29/1 | 29/1 |
| test_simplify_radsimp | 24/2 | 24/2 |
| test_simplify_simplify | 55/9 | 55/9 |
| test_simplify_trigsimp | 26/2 | 26/2 |
| **total** | 217/51 | **218/50** |

The only changed test is `test_core_basic::test_has`, which went from fail
to pass. The `has` tests of SymPy's `test_expr.py` were ported too: 10 of
11 pass. The remaining failure is `test_has_physics`, where
`FockState((x, y)).has(x)` is a Python physics class.

The ODE corpus scores 315/289/270 after the fix, the same as the baseline.

## Remaining gaps

- `FockState(...).has` (Python physics class).
- `srepr(WildFunction('F'))` prints `Wild('F')`. The WildFunction mark is a
  Python predicate, so printing it needs a native flag on `Wild`.
- Not addressed in this lot:
  - the cost of ODE matching;
  - the native matcher declining Float, `zoo`/`oo` and `Sum`.

## Validation

- Suite (`run_suite.py`): SUITE CLEAN, 60 580 passed, 6 xfailed,
  42 skipped, 60 628 collected (baseline 60 601 + 27 new tests).
- `cargo test --release`: 873 + 1 + 0 + 5 passed (one new unit test).
- `cargo clippy --release -- -A clippy::all -D clippy::correctness`: clean.
- ODE corpus: 315/289/270, the same as the baseline.
