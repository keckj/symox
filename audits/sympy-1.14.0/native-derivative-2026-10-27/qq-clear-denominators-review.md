# Native QQ DomainMatrix denominator preparation — frozen delivery

Isolated source: `../qq_clear_denoms`; immutable source base: `../qq_clear_denoms_base`.
Base is the previously delivered Float0 EX v18 source (including FF), not the
transient main extension and not main's subsequent Derivative/conditional-integral
changes. Integrate only the five-file patch; do not replace the source tree.
Main and all earlier frozen stages were left unchanged. No commit, local Cargo,
maturin, shared installation, or timing campaign was performed.

## Provenance

- Base: 2777 source files, manifest SHA256 `27690544f92b547bc5292395dd58ff2e16fe837e626340fcc1390809b7a9d62b`.
- Final: 2778 source files, manifest SHA256 `8e7a4c81ee29ad5eb8904f5b61c8e4bfa44d37548a4f42eb0a645247a92b95f4`.
- The final manifest was independently reproduced on Demeter with the same count/hash.
- Patch SHA256 `fbfc73ba702cff6d4c77084d5b3652f8d5b2dc6720c7913e027b50c4a4673d86`; `git apply --check` succeeds on the immutable base.
- DEBUG extension SHA256 `94fed4c591da847b41b4276ee3b4bad45076ea3a2cbd223d807766e4e93ddd78`.
- Private Demeter tree: `/home/keckj/symox-ci/wt/qq-clear-denoms-20261026`.
- Built with `PYO3_PYTHON=/home/keckj/symox-ci/venv/bin/python`, private Cargo target `.../targets/expression-sparse`; the resulting library was copied only into that private tree's Python package.
- Build v2 used `source-v2.sha256`; subsequent changes were Python protocol guards and tests only. Both Rust files are unchanged since that successful build. The two core tests used the identical final core file.
- Per-file base/final hashes are in `file-hashes.json`.

## Change

The Rust `fraction_free::clear_denominators` kernel computes exact positive LCMs
and exact scaled integer numerators. It accepts ordered occupied rows rather
than a rectangular dense copy. Global mode shares one LCM across all rows;
rowwise mode retains one factor per occupied row. Empty matrices, empty dense
rows, arbitrary-size integer numerators/denominators and interruptions are covered.

The binding validates exact QQ scalar objects, builtin row containers, legal
indices and matching shapes before reading through CAPI. It preserves sparse
row/column insertion order and the original builtin integer key objects.
When a factor is1 and convert=False, source coefficient PyObjects are reused.
Other numerator elements are newly constructed in QQ or ZZ as required.
There is no floating conversion or added Python mathematics.

The Python adapter verifies stock class/method/operator/descriptor/domain
protocols without invoking metadata descriptors. It retains the original
DomainMatrix factories for diagonal denominators and numerator assembly.
These remain authoritative for output classes, formats, zero objects, domains,
and construction behavior. Custom classes, instance overrides, descriptors,
operators, conversion hooks, truth-valued convert arguments and foreign rows
retain the complete previous implementation. Signals are checked throughout
parsing, arithmetic and publication, and errors propagate directly.

## Validation on Demeter

- Final new file: **97 passed**, `qq-clear-denoms-final.log`.
- Exact comparison against local SymPy1.14, plus certificates D*A=B after domain
  unification: dense/sparse, rowwise/global, convert=False/True, empty/zero
  dimensions and matrices, signed fractions, noncoprime denominators and
  coefficients up to251 numerator bits / arbitrary denominator sizes.
- Native route spies prove both methods execute the native preparation and
  eight whole-call AUTO/CD rref/rref_den cases use it.
- Source coefficient aliases, complete output alias partitions, builtin key
  identities and noncanonical sparse insertion order match the old path.
- Callback traces match the old implementation for matrix/representation/
  domain/conversion/coefficient-truth/convert-truth/subclass protocols.
  Metadata descriptors are not probed during eligibility; callback exception
  object identity is preserved.
- Global empty sparse shapes `(2**90,2**91)` retain the old complete fallback.
- Existing selected FF + DomainMatrix tests: **376 passed**, in combined v2
  (454 passed overall plus two faulty new test witnesses, subsequently fixed).
- Two targeted Rust tests passed: rowwise/global exact preparation including
  empty input, and direct interruption propagation.
- No whole-suite or performance claim is made. The parent owns the integrated
  release and full Demeter suite.

The archived intermediate combined run contains two new witness errors:
a setterless rep property made the old factory fail as well, and a global
convert=True test incorrectly expected a conversion callback although the old
path performs none. Final tests use a read/write metadata property and an
actual convert=False conversion respectively. No algorithm/result assertion
was relaxed, and final97 all pass.

## Deliberate retained scope and gaps

Native eligibility is stock canonical QQ only. ZZ, RR, finite fields,
polynomial/fraction domains, subclasses, foreign or manually poisoned entries,
non-bool convert, custom protocols and unsupported native index capacities use
the original complete code; they are not counted as native gains.

An explicit empty sparse row keeps its original rowwise ValueError. Explicit
stored sparse zeros also stay on the original path: SymPy and the old Symox
code apply dense dup_strip during some conversions, which can shift a later
entry to an earlier saved column. For example sparse row `{2:QQ(0),0:QQ(1,2)}`
under convert=False produces the coefficient at column2. This is a preexisting
noncanonical-representation issue, not claimed as D*A=B parity for this port.
Canonical zero matrices/zero rows are certified normally.

The AUTO algorithm selector still performs density/denominator inspection and
its selection LCM in Python. Clearing denominators is now native, but the
whole AUTO call is not claimed to be entirely native. No speedup ratio was
measured or promised. No pivot, FF elimination or accepted-domain algorithm
was changed.
