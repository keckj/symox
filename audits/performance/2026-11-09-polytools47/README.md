# Structural cancel contracts, wave 47, 2026-11-09

The reference is installed SymPy 1.14.0. This lot ports the structural protocol
of `sympy.polys.polytools.cancel`, rather than inventing a `deep` cancellation
mode. Selected regressions do not establish global polynomial parity.

## Reference protocol and implementation

`cancel` first calls `allowed_flags(args, ['polys'])`. This validates unknown
option names as well as forbidden flags, before numbers, relations, booleans or
containers return. Thus **`deep=True` is an invalid cancel option** and raises
`OptionError`, even on a number or Piecewise; `auto=False` raises `FlagError`.
The private `_signsimp` argument is a separate control, consumed before this
validation and tested by truth value afterwards. The binding retains `None`,
zero and custom `__bool__` behavior, including exceptions.

The expression path applies optional `signsimp`, then
`factor_terms(radical=True)`, before polynomial construction. Piecewise forces
the reference's structural fallback. A sum or product is partitioned into
commutative Piecewise-free arguments and the remaining arguments. The first
partition is combined before cancellation; this permits cancellation across
several factors outside a Piecewise. Other nodes use the reference's descendant
protocol: cancel the next expression as a whole and skip descendants once
handled; traverse relational/boolean/container nodes without treating them as
scalar expressions. Branch conditions and nested Piecewise expressions keep
the native reconstruction and canonicalization paths.

These steps are in `crates/symox-core/src/poly/api.rs`. The binding handles
Piecewise expressions and plain numerator/denominator tuples before specialized
polynomial-domain paths, consistent with the reference's forced structural
fallback. A Piecewise tuple returns a transformed two-element tuple, rather
than an ordinary `(coefficient, numerator, denominator)` polynomial result.

Nonfinite atoms have a different role. SymPy's `sring` accepts `oo`, `-oo`,
`zoo` and `nan` as formal polynomial generators. The previous blanket refusal
prevented cancellation and polynomial reconstruction on expressions containing
them. `convert::from_exprs_with_nonfinite` now exposes that behavior only to
`cancel`: arithmetic on the polynomial coefficients remains exact, and native
expression constructors evaluate the nonfinite values when the polynomial is
reconstructed. Other polynomial readers retain their rejection of these atoms.
There is no blanket replacement of an expression containing `zoo` by `nan`.

No Python arithmetic is added. The old Rust `poly::cancel` three-argument entry
point remains available to PRDE callers. Existing algebraic, Gaussian,
expression-domain, finite-field and inexact polynomial paths are retained.

## Evidence and reproduction

`tests/test_polytools_wave47.py` contains **91 differential cases and one exact
branch-value certificate (92 tests)**. The cases compare structural values,
Float precision, Poly domain/generators and exception classes. Controls cover
Gaussian cancellation, an explicit algebraic extension, modulus 2 and Float
coefficients at 80 decimal digits. Structural cases cover Piecewise products,
sums, functions, powers, reciprocals, nested branches and conditions, tuples,
noncommutative factors, held nonfinite products, ordinary functions and
derivatives. The nonfinite group includes the rational expression encountered
by the `abaco2_similar` Lie heuristic, without invoking the full ODE solver.

The frozen wave-46 baseline at `/tmp/symox-wave47-baseline/python` gives
**47/91 SAME** in [before.json](before.json). Its identity is stored alongside
that baseline. The snapshot is produced by this lot's small probe:

```sh
PYTHONPATH=/tmp/symox-wave47-baseline/python /home/keckj/symox-ci/venv/bin/python \
  audits/performance/2026-11-09-polytools47/probe.py \
  audits/performance/2026-11-09-polytools47/reproduced-before.json
PYTHONPATH=python /home/keckj/symox-ci/venv/bin/python \
  audits/performance/2026-11-09-polytools47/probe.py \
  audits/performance/2026-11-09-polytools47/after.json
PYTHONPATH=python /home/keckj/symox-ci/venv/bin/python -m pytest -q \
  tests/test_polytools_wave47.py tests/test_polytools_wave46.py
```

On the first release build, the [post-build probe](../2026-11-09-wave47/raw/attempt1-polytools47-after.json)
and [probe log](../2026-11-09-wave47/raw/attempt1-polytools47-probe.log) give
**91/91 SAME**, versus **47/91** on the frozen baseline: **44 newly equal cases,
zero regressed cases**. Case names and all 91 frozen reference snapshots are
unchanged. The 92 polytools regressions are among the passing cases of the
[intermediate focused run](../2026-11-09-wave47/raw/focused-attempt1.log);
that combined run also recorded six failures in the separate series lot, so it
is not a passing global validation result.

The combined suite on the third build found two strict XPASS results in
`test_public_heuristics_match_sympy[False-abaco2_similar]` and
`test_public_heuristics_match_sympy[True-abaco2_similar]`. Their existing xfail
reason was cancellation with `zoo`; the same inputs and assertions now pass.
Only their two obsolete xfail markers were removed. The final-source focused
and complete checks confirm these two resolved expected failures.

## Final focused and suite results

The [final probe](../2026-11-09-wave47/raw/polytools47-after.json) and
[probe log](../2026-11-09-wave47/raw/polytools47-probe.log) again give
**91/91 SAME**. All 91 case names and frozen reference snapshots match the
baseline: **44 newly equal cases, zero regressions**. The 92 polytools tests are
included in the [596 passing focused tests](../2026-11-09-wave47/raw/focused.log).

The [complete final suite](../2026-11-09-wave47/raw/suite-final-summary.json)
is clean: **61,665 passed, 11 xfailed, 41 skipped and zero failures**. All
**61,717 collected cases** in 709 files were executed exactly once, with none
missing or unexpected. The
[focused](../2026-11-09-wave47/raw/focused-identity.json),
[lot](../2026-11-09-wave47/raw/lots-identity.json) and
[suite](../2026-11-09-wave47/raw/suite-identity.json) identities agree:
source manifest `8acf4d47316eb81e8da94393bfa504be776c58c24d48890dd50b6dd07d2751c0`,
extension `4268d4ec5662bea632ceb2cc460847cd7dda2758d995718635e752e23adab1fe`.

Final Rust checks pass: 880 unit tests (one ignored), one integration test,
five doctests and correctness clippy. All 249 algfactor, 139 algpoly and
75 ratint2 recorded payloads are unchanged against wave 46; the per-case
comparison is retained in `../2026-11-09-wave47/raw/corpora-deltas-early.json`.
The remaining corpus replay and resource controls are complete;
their results and explicit unresolved resource gaps belong to the
[combined wave-47 audit](../2026-11-09-wave47/README.md). The suite and selected
probes establish the scoped results above, not global mathematical parity.
There is no speedup or asymptotic-complexity claim for this correctness lot.

## Remaining scope and limits

- `together(deep=True)` and other functions' valid `deep` controls are outside
  this lot; they must not be inferred from cancel's `OptionError` contract.
- Cancellation does not certify values at excluded poles. The branch certificate
  samples exact rational values away from the original denominators' zeros.
- Every combination of domain, generator, noncommutative and Piecewise input
  has not been audited. The regression controls preserve established paths;
  they do not establish all options/domain combinations.
- A specific pre-existing error-contract gap remains outside the 91 probes:
  `cancel(zoo, domain='invalid-domain')` raises `OptionError` in SymPy 1.14,
  while both the frozen wave-46 baseline and this lot return `zoo`. Unknown
  option **names** and disallowed flags are now validated before early returns;
  that does not establish validation of every option **value** on nonfinite
  inputs. This result is recorded as an open gap, not a reference match or a
  regression introduced by the structural cancellation changes.
- General host-defined expression classes, matrices/containers, interruption
  and resource limits remain outside this structural scalar lot.
- `same_root`, private ring APIs, algebraic-field representatives and the
  remaining polynomial gaps tracked by the global contract are unchanged.
