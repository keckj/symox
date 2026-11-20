# Wave 47: ordered noncommutative reference matching

This lot extends the Rust reference matcher in
`crates/symox-core/src/matches.rs` to native noncommutative scalar expressions.
It preserves the wave-46 Float inverse convention and the admission guard
for `evaluate(False)`.

## Reference protocol

The implementation follows SymPy 1.14.0 `Mul.matches`, including:

- distinct commutative/noncommutative operand guards;
- `args_cnc` splitting, including the reference's sign splitting and the
  comparison of commutative factor lists when bindings are empty;
- matching the commutative prefix before expanding noncommutative powers;
- positive exact integer power expansion with a bounded allocation;
- the LIFO state agenda, preferring advancement to the next pattern node
  before extending the current Wild's slice;
- shared insertion-ordered Wild ranges and accumulated bindings across
  backtracking, without replacing them by a different search algorithm;
- literal-factor matches, compound-factor matches with fresh local
  bindings, exclusions, property callbacks, and repeated Wild checks;
- the reference's return of `None` when its noncommutative Mul branch
  finishes with empty bindings, even for an identical literal product.

Two details affect the selected bindings. Repeated Wild slices are compared
by zipped prefixes, and a later successful slice overwrites the earlier
binding. Thus `A*B*A*C` matched against `w*B*w` gives `{w: A*C}` in the
reference. The port preserves this behavior rather than substituting a
different matching rule. That binding does **not** reconstruct the target:
substituting it yields `(A*C)*B*(A*C)`, rather than `A*B*A*C`. This is an
existing reference-protocol behavior already present in the Python fallback,
not evidence that the binding is mathematically valid. It remains a known
mathematical gap, tracked separately from a differential `SAME` count.
Also, `Mul.as_base_exp` must refuse to regroup
`A**2*B**2` as `(A*B)**2` for distinct noncommutative bases, while still
regrouping `x**2*A**2` as `(x*A)**2` when `x` commutes. Both cases are covered.

The expansion helper reproduces SymPy's `TypeError` for positive Rational
or Float exponents in the noncommutative factor sequence. It reports an
owned Rust `ProtocolTypeError`, carried by the existing matching error type.
The shared PyO3 wildcard-error converter recognizes that type by downcast;
existing Python callback errors retain their previous conversion. No new
mathematical calculation is performed in Python. The same typed propagation
is tested through public `collect`.

Guard order is deliberate: a failed commutative prefix comparison returns
`None` before any positive noninteger expansion error. Positive integer
expansion is declined before allocation when it exceeds 200,000 factors;
the existing search-step budget also covers the noncommutative agenda.
A declined computation remains a fallback request, never a no-match result.

## Regression coverage and baseline

`tests/test_native_matcher_wave47.py` adds 120 cases: 52 matching pairs under
both `old` flags, three conservative fallback witnesses, eight native
numeric-expansion error witnesses, prefix-guard ordering, shared `collect`
error propagation, symbolic-exponent fallback errors, huge-power admission,
and property callback error propagation. It checks binding insertion order
and exact expression representations, not only whether any match exists.
One new Rust unit test covers the agenda, empty-binding convention,
repeated-Wild choice, typed errors and guard ordering.

The pre-existing `test_noncommutative_declines` regression retains its exact
input `A*B` against `Wild('w')*B`, with a commutative Wild. It is renamed and
updated to require native `None`, public `None`, and agreement with SymPy;
this is the reference's rejection of a commutative Wild consuming an NC
prefix. Separate fallback cases retain the decline contract.

The frozen baseline is the completed wave-46 runtime at
`/tmp/symox-wave47-baseline/python`, extension SHA256
`63651dbe6ec1fd4b925345de1aa08d0061aebe87e8d84840438eb605f29eb837`.
The lightweight [baseline probe](probe_baseline.py) records
[56 witnesses](baseline-witnesses.json):

| Group | Cases | Native served before | Public equal before |
| --- | --- | --- | --- |
| Ordered matching | 52 | 0 | 51 |
| Positive numeric expansion errors | 4 | 0 | 0 |

The additional public matching difference is `x**2*A**2` against `w**2`:
SymPy binds `w` to `x*A`; the old Python protocol binds it to
`sqrt(x**2*A**2)`. The Rust `as_base_exp` handling above covers this witness.
The old Python expansion helper coerces exponents with `int`, so it returns
a match or `None` where SymPy raises `TypeError` for the four numeric error
witnesses. The native typed-error path corrects these admitted calls.

Reproduce the lightweight baseline probe:

```
PYTHONPATH=/tmp/symox-wave47-baseline/python \
  /home/keckj/symox-ci/venv/bin/python \
  audits/performance/2026-11-09-matcher47/probe_baseline.py
```

## Frozen build: passing tests and complete matcher corpus

The final frozen source SHA256 is
`8acf4d47316eb81e8da94393bfa504be776c58c24d48890dd50b6dd07d2751c0`,
and the extension SHA256 is
`4268d4ec5662bea632ceb2cc460847cd7dda2758d995718635e752e23adab1fe`
([build identity](../2026-11-09-wave47/raw/build-input.json)).
The integrated [focused run](../2026-11-09-wave47/raw/focused.log) passed
**596 tests**, including all **120** matcher regressions
([command receipt](../2026-11-09-wave47/raw/focused.command.json)).

The [complete suite](../2026-11-09-wave47/raw/suite-final-summary.json)
records **61,665 passed, 11 xfailed and 41 skipped**. All **61,717** collected
cases across 709 files finished, with collection/execution agreement and no
missing or unexpected cases. These test receipts identify the same frozen
extension.

The corresponding 56-witness probe gives:

| Group | Public equal before | Public equal on frozen build | Native on frozen build |
| --- | --- | --- | --- |
| Ordered matching | 51 / 52 | 52 / 52 | 52 decided results equal to SymPy |
| Positive numeric expansion errors | 0 / 4 | 4 / 4 | 4 direct TypeErrors equal to SymPy |

The four errors are native exceptions, not declines. The probe's
`native_handled` field only counts returned `(True, value)` pairs, so that
field stays zero for this error group; its row-level `native.type` and
`native.message` prove the typed propagation. Public agreement includes both
exception class and message. No public mismatch remains among these 56
witnesses. Receipts: [witness rows](../2026-11-09-wave47/raw/matcher47-after.json),
[probe command](../2026-11-09-wave47/raw/matcher47-probe.command.json),
and [probe log](../2026-11-09-wave47/raw/matcher47-probe.log).

The complete 2,311-call matcher corpus was then replayed on that same frozen
extension:

| Implementation | Wave 46: equal / different / declined | Wave 47: equal / different / declined |
| --- | --- | --- |
| Native reference matcher | 2228 / 0 / 83 | 2289 / 0 / 22 |
| Public `Expr.match` | 2311 / 0 / 0 | 2311 / 0 / 0 |
| Existing Python protocol alone | 2309 / 2 / 0 | 2309 / 2 / 0 |
| Older `pattern_match` | 2059 / 252 / 0 | 2059 / 252 / 0 |

Case-by-case inspection finds exactly **61** changed rows, each changing
only the native field from a decline to a returned result. The public,
Python-protocol and older-matcher fields are unchanged on all 2,311 rows
after removing timings. All 61 newly served cases belong to `test_match`:
its native counts move from 182 equal / 66 declined to 243 equal / 5 declined.
The 1,563 ODE-source cases retain 1,546 native equal / 17 declined; the 500
random cases remain native equal. No native difference or binding-order
difference is reported. Receipts: [wave-46 rows](../2026-11-09-wave46/raw/matcher-after.json),
[final rows](../2026-11-09-wave47/raw/matcher-after.json),
[final report](../2026-11-09-wave47/raw/matcher-report-after.json),
[replay command](../2026-11-09-wave47/raw/matcher.command.json), and
[report command](../2026-11-09-wave47/raw/matcher-report.command.json).

The integrated [Rust checks](../2026-11-09-wave47/raw/cargo-test.log) passed
**880 unit tests** (one ignored), **1 integration test** and **5 doctests**;
the new noncommutative matcher unit test passed. The
[correctness clippy run](../2026-11-09-wave47/raw/clippy.log) completed
successfully, with no correctness error. Both
[Rust-test](../2026-11-09-wave47/raw/cargo-test.command.json) and
[clippy](../2026-11-09-wave47/raw/clippy.command.json) receipts identify the
same frozen extension.

The matcher's mathematical reservation remains independent of these passing
checks: the repeated-Wild example above agrees with the reference while
returning a binding that does not reconstruct the target. The witness
agreement is a protocol result, not a proof of valid substitutions or
solution completeness.

These complete matcher receipts establish the recorded protocol behavior
over the replayed cases; they do not establish global mathematical matching
parity. Other mathematical corpora and resource-boundary controls are
reviewed separately by the integrating agent. Compilation and post-change
testing were centralized by that agent; no build or heavy corpus replay was
run by the matcher agent.

## Remaining gaps

- The repeated-Wild prefix comparison can return a binding that does not
  reconstruct the expression, as in the example above. Compatibility with
  that SymPy result must not be counted as mathematical binding validity;
  this known gap is independent of the differential agreement counter.
- Host classes, matrix nodes, binders, nonfinite values and the previously
  unsupported node classes retain their admission rules.
- Unknown commutativity and nonliteral/nonpositive noninteger powers in an
  NC factor sequence conservatively decline. A symbolic unconstrained
  exponent retains the existing reference-like `TypeError` fallback.
- Positive nonliteral exponents such as `pi`, `sqrt(2)` or a positive Symbol
  remain outside the native expansion model. The existing Python fallback
  coerces with `int`; its values or exception messages can disagree with
  SymPy's sequence-repetition `TypeError`. This gap is recorded, not counted
  as a native no-match or a fixed public case.
- Public `matches` calls carrying nonempty `repl_dict` and
  `evaluate(False)` contexts still use the established Python protocol;
  they are not covered by the new native numeric-expansion error path.
- Power-expansion and search budgets retain complete fallback behavior.
  They are not mathematical completeness claims.

The lot aims to preserve the tested SymPy matching protocol and accepted
error behavior. It makes no global matching-parity or performance claim.
