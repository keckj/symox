# Typed RR/CC composite equality — frozen follow-up

The CC[t] equality failure is a real native API defect. SymPy returns true for both `a == a` and `(a/R.one) == a`, with exactly the same selected domain and dyadic coefficient values. Native returned false even when the two stored expressions were structurally identical. The follow-up fixes this equality dispatch without changing representations, arithmetic, precision or hashes.

New independent stage: `/home/keckj/.cache/symox-validation/native-composite-equality-20261027/work`. Neither main nor the frozen inexact_linear / CC-rounding packages was modified. Patch: `composite-equality.patch` beside this report.

## Provenance and changed files

Base is the complete frozen CC-rounding stage, checked against every entry of its 2893-file manifest before editing. The lineage is immutable r1b -> inexact_linear -> CC-rounding -> this follow-up; no current-main or transitory binary was copied.

- Base manifest: `71aa6054ad10fc3261033189fa93bfbc7b0cd77193a35b5cb8b5656591a012be`.
- Final manifest (2894 files): `114441b755b5f024196fa89430d778a8fe99f390dfb1b92499ff2e8deac58d93`.
- Patch: `8a1e42061a23208c701edd180949ecdcdac11cf1bff5d612dcf2d48656794f39`.
- Private DEBUG native: `16b37897a3f215fb26e0b25f89d4d316995ea63ddaa539c6fb4280403768aa8b`.

Two files only:

1. `crates/symox-py/src/domains.rs`: 17 added Rust lines in `PyDomainElement::__richcmp__`.
2. New `tests/test_inexact_composite_equality.py`: 107 tests.

The binding base hash is `48c4e97ca2d7e547851d2943e93e439e9b714b1c1ee03de38f204e66390aa7db`, final `d3c99d729f09d3bdaf0cb3887f99087997a6be7f593a9454961606dfe9b6da37`. Merge only this hunk after the inexact_linear binding update, retaining other queued BigInt/current binding changes. No core coefficient algorithm or Python production code changed. `identity.json` records all hashes, complete Demeter source verification passes, and `git apply --check` passes on base.

## Reference/API diagnosis

`initial-probe.log` compares p1/53/139/140/211 for CC[t] and CC(t), with `a=(1+I)*t+2`. All coefficients are exactly representable dyades even at p1. In both implementations, a and a/one have equal parents. Native's expressions and hashes were also identical, but both `a==a` and `a/one==a` were false. Reference equality is true in every case.

`protocol-probe.log` extends the investigation to constants, different parents, raw expressions and the separate Python ring wrapper. The new tests extend it to RR as well, using each native ground's exact stored dyadic coefficients to construct the independent reference coefficients. Float input values are read as mantissa/exponent or exact rationals, without a text round-trip.

The first divergence is the binding's equality branch, not a coefficient value/type/precision difference. Every `ElementValue::Inexact` previously entered a scalar-only path. If `exact_numeric_parts` failed because an expression contained parameters, it returned false for Eq and true for Ne, including a comparison with itself. The immutable Rust `DomainElement` already implements reflexive typed Eq and matching Hash.

## Minimal correction and invariants

The original scalar comparison path remains first, preserving numerical comparisons and cross-parent constants. If scalar extraction is unavailable, an inexact composite compared with another native DomainElement of exactly the same domain now delegates Eq/Ne to the existing immutable typed value equality. Other comparisons and foreign-object protocols retain their previous paths. Ordering still returns NotImplemented here.

There is no simplification, cancellation, numerical approximation, precision coercion or change to a stored coefficient. EX domains are not selected by this branch. Distinct coefficient dyades remain distinct, including values whose one-bit printed strings collide. The selected ground, precision, parameters and order are compared as part of the domain.

Hashes remain unchanged: numeric constants keep their existing scalar-compatible hash; nonnumeric composite elements keep the existing hash of domain and immutable value, matching the Eq operation now exposed. Dict/set lookup, equality symmetry, self-equality, copies and pickle reconstruction are covered. No reference mutable-polynomial hash anomaly is copied.

Callbacks are preserved: custom public equality overrides still run; native typed comparison does not newly call public parent/as_expr overrides; foreign sympify callbacks and their established error handling remain in place. Conversion, zero-division and ordering paths were not modified.

## Demeter validation

Private directory: `/home/keckj/symox-ci/wt/composite-equality-20261027`. Private target: `/home/keckj/symox-ci/targets/composite-equality`. Cargo builds with the dedicated Python3.10 venv; only that directory's .so was copied to that directory's Python package. No maturin/shared installation or benchmark.

- `cargo build -p symox-py`: passes (`composite-equality-build.log`). No core source was changed, so no additional core test suite was repeated.
- **107 new tests pass** (`composite-equality-python-v2.log` remotely; included in the final log).
- Final focused run: **729 passed**, consisting of 107 new and 622 existing tests (`composite-equality-python-final.log`). It covers inexact composites, EX Float0 payloads, inexact followup, all 203 frozen inexact_linear tests, domain division and both Heurisch regression files.
- All 2894 final source entries checked on Demeter (`composite-equality-source-check.log`).

The tests cover RR/CC at p1/53/139/140/211, polynomial/fraction fields, one/two parameters, nontrivial fractions, equality and inequality against reference peers, zero/one/scalars, other grounds/precisions/generators/orders, hashes and dict/set keys, conversions, copy/deepcopy/pickle, EX Float0, different coefficient bits, callbacks and errors. Presentation metadata p' is compared with the reference while stored coefficient mantissas remain within p; no assertion replaces calculation precision with presentation precision.

## Explicit remaining gaps and reference limitations

- A native `Domain.convert` element and the separate Python `K.ring.from_expr` wrapper still expose different native parent/protocol representations and compare unequal, whereas reference peers compare equal. This distinct cross-wrapper interoperability gap remains explicit. The patch does not opportunistically coerce either side or bypass foreign protocols to hide it.
- Nonnumeric ordering is unchanged: native `a<a` raises TypeError, while reference `ar<ar` returns false. This follow-up repairs Eq/Ne only.
- The reference mutable PolyElement hash cache can produce unequal hashes for equal ar and ar/one in the initial probe. Native immutable values have equal hashes and retain Python's equal-values/equal-hashes invariant.
- SymPy 1.14 PolyRing deepcopy/pickle raises `RuntimeError: dictionary changed size during iteration` in this environment. The exploratory test log and independent `reference-copy-probe.log` preserve both failures. Native copy/deepcopy/pickle round-trips are verified against the original reference mathematical element and full native stored expression; the reference failure is not imitated or counted as a new speed/capability gain.
- This is equality of existing typed values in the same domain, not a general algebraic-equivalence or cross-parent coercion engine. It makes no universal polynomial/domain parity claim.

No existing test or frozen delivery was edited. The parent retains full integrated release/Demeter validation responsibility.
