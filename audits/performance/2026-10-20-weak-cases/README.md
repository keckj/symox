# Follow-up fixes for weak performance cases

Reference: SymPy 1.14.0. The predecessor is the
[first patch round](../2026-10-20/README.md), source `0a6b26f` with package
sources from `2b1dab9`. The new wheel records production source `f301412`.
The validation source additionally includes `8a0c220`, which changes one test
assertion and does not change package or Rust sources.

## Changes

- Standard whitespace no longer causes an extra tokenization pass before the
  parser. Unusual whitespace still follows the token-aware normalization that
  protects quoted strings.
- Lexical tokens for at most 16 sources of at most 4,096 characters are reused.
  They are immutable; each parse gets a fresh mutable list. Namespace lookup,
  transformations, user functions and mathematical evaluation run on every
  parse. Transformations also avoid rebuilding token streams when the relevant
  operator/notation is absent. Inputs outside the cache limit retain the parser.
- Concrete native scalar expansion calls its native method directly, using the
  same options and work/term budgets. Array, tensor, subclass and Dirac delta
  dispatch retains its established path. Default calls avoid empty hint copies.
- Native numeric conversion avoids trying algebraic coefficient conversion
  before primitive Python numbers. Exact built-in numeric operands avoid sparse
  marker lookup. Numeric subclasses retain their marker/priority protocols.
- Scalar wrappers avoid repeated operation-name parsing and unnecessary numeric
  priority lookups. Type classification uses identity, without calling foreign
  metaclass equality.
- Missing differential-field registries are checked without creating an
  attribute exception for every scalar result. Registered fields retain the
  free-symbol check and their callable subtype.
- Exact `Add` reuse builds an injective byte key from Python object addresses
  in one native pass. Retained owners prevent address recycling; no user hash
  or conversion hook is invoked. Cache limits and mathematical eligibility
  are unchanged.
- Powers of concrete involutive permutations reduce by exponent parity after
  verifying every image satisfies `p[p[i]] == i`. This covers arbitrary signed
  integer exponents. Subclasses retain inversion and constructor callbacks.

## Evidence

`after.json.gz`, `after.md` and `after-summary.json` contain the new run, readable
timing report and machine-readable summary. `comparison.json` compares it with
the predecessor's retained `after.json.gz`, checking identical inputs, timing
sources and timing configuration. `after-package.json` contains all 611 package
file hashes and the wheel/build fingerprint.

The final report contains 589 paired trials across 137 workloads, with seven
fresh-process trials for 46 workloads. The refinement threshold is 5×, so the
newly improved cases retain seven trials too. All 132 comparable workloads
remain verified; the same five inconclusive/timeout cases remain excluded.
Seventeen workloads have a first-call or repeated-call speedup at most 1.5×.
The core Python suite passes 30,045 tests with five documented skips and no
failures or xfails; the final focused suite passes 1,618 tests.

The final logs are `weak-core.log.gz`, `weak-regressions-final.log.gz`,
`weak-build.log.gz`, `weak-screening.log.gz` and `weak-refinement.log.gz`.
The core Python suite includes Trisplyne and excludes `tests/extra`; optional
extra backend and downstream package qualification were not rerun. The Rust
core sources are unchanged from the predecessor's 528-pass core run.

An initial focused run exposed an overstrict new test requiring Python object
identity for NaN. The predecessor produces equal but distinct wrappers too.
The assertion now checks the established value semantics; production code did
not change for this correction. `weak-regressions.log.gz` retains that first
run, separately from the final successful focused run.

Reproduce using the [same runner](../../../benches/README.md) with the immutable
package matching `after-package.json`. Timings are serial on CPU 0; tests run
on CPU 1. The shared host is not isolated. Mathematical result checks precede
timing rankings, and excluded validations/timeouts remain excluded. These
measurements do not renew the full API or downstream parity certificates.
