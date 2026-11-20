# Polytools edge contracts, wave 46, 2026-11-09

This lot fixes ground-domain and list polynomial operations in Rust, retaining
exact rational arithmetic and inexact intermediate precision. The reference is
the installed SymPy 1.14.0. It does not establish global polynomial parity.

## Scope

- Ground `gcd`, `lcm` and `cofactors` distinguish the numerical-domain fallback
  from polynomial normalization with generators. Thus `gcd(3., 9.)` is a Float
  one, `lcm(3., 9.)` is a Float 27, and `cofactors(8, 6)` is `(2, 4, 3)`.
  Rational scalar cofactors remain exact. `cofactors(0, 0)` raises
  `ZeroDivisionError`, whereas the polynomial operation with explicit generators
  returns three zeros, following the reference.
- The arithmetic precision is the greatest input Float precision. Ground output
  follows the reference RealField decimal-digits/bit-precision roundtrip only
  after arithmetic. Regressions include 1, 2, 52, 54, 80 and 100 input bits,
  a Float at 80 decimal digits with exponent 400, and a mixed large rational.
- Native list functions jointly infer generators and domains, preserve requested
  `Poly` output, support explicit domains, `field=True`, modulus and selected
  algebraic extensions, and distinguish empty and generator-free requests.
  A pair of native `Poly` inputs retains its coefficient domains, including
  finite fields, when generators are inferred, identical or permuted. Distinct
  finite fields in such a pair raise `UnificationFailed`. Changed generators,
  mixed expression/Poly lists and lists of three or more Polys reconstruct
  their expression domains, as in the reference.
- `gcdex` and `half_gcdex` with inferred symbol generators now retain `auto`,
  `domain` and `polys` options. In particular, `auto=False` over ZZ raises a
  domain error instead of silently computing in QQ.

`python/symox/__init__.py` only admits the two list operations to the new native
exports. The list binding reuses `_poly_constructor._options` for existing input
and option validation. Arithmetic and polynomial construction are in Rust;
existing specialized algebraic and expression-domain paths are retained.

## Evidence and reproduction

`tests/test_polytools_wave46.py` has 111 public differential cases, 35 precision
cases, one exact cofactor/Bézout certificate and two iterator-admission cases.
Snapshots compare values, Float precision, Poly generators/domain and exception
classes. The finite-field cases compare both the polynomial and its domain.
The total is **149 tests: 111 differential + 35 precision + 1 certificate +
2 iterator cases**. The differential set consists of the original 91 cases and
20 generator/domain cases added during review.

`before.json` records **24/91 SAME** on the wave-45 extension. The probe's
`--baseline` switch restores only the old Python admission
`gcd(list(seq), ...)` / `lcm(list(seq), ...)`, since the new native exports are
absent from that extension; every calculation still runs on the baseline.

```sh
PYTHONPATH=python /home/keckj/symox-ci/venv/bin/python \
  audits/performance/2026-11-09-polytools/probe.py \
  audits/performance/2026-11-09-polytools/reproduced-before.json --baseline --initial-only
PYTHONPATH=python /home/keckj/symox-ci/venv/bin/python \
  audits/performance/2026-11-09-polytools/probe.py \
  audits/performance/2026-11-09-polytools/after.json
PYTHONPATH=python /home/keckj/symox-ci/venv/bin/python -m pytest -q \
  tests/test_polytools_wave46.py
```

The initial implementation passed all original 129 tests and 298 focused tests
on the coordinating agent's merged build. A subsequent targeted review found
that preserving original Poly objects ignored explicit generator changes.
`before-generators.json` freezes the intermediate wave-46 implementation:
**6/20 SAME**, including failures with coefficient or finite-field changes and
with generator order. The initial `before.json` is retained unchanged.
These two snapshots describe different builds and must not be combined into a
single pre-fix score.

The fix follows the reference's special case for exactly two Poly inputs.
Identical generators or a permutation retain their domains and apply the
requested order; replacing or adding generators reconstructs expressions in the
new generators. Lists of three or more Polys also reconstruct, as the reference
does. The 20 added cases cover QQ/GF with identical, replacement and additional
generators, a multivariate QQ permutation, a domain coercion error, and triples.
The complete focused file now contains 149 tests.

```sh
PYTHONPATH=python /home/keckj/symox-ci/venv/bin/python \
  audits/performance/2026-11-09-polytools/probe.py \
  audits/performance/2026-11-09-polytools/after-generators.json --generators-only
```

The final [probe results](after.json) and [probe log](../2026-11-09-wave46/raw/polytools-probe.log)
are **111/111 SAME**. Comparing case names and frozen reference snapshots gives:

| Case group | Its baseline | Final build | Newly equal | Regressed |
| --- | --- | --- | --- | --- |
| Original 91, [wave-45 baseline](before.json) | 24/91 | 91/91 | 67 | 0 |
| Additional 20, [intermediate wave-46 baseline](before-generators.json) | 6/20 | 20/20 | 14 | 0 |

Every reference snapshot is unchanged, the two case groups are disjoint, and
the final probe contains each of their 111 names once. The baseline columns
still describe separate builds; the improvements are case comparisons, not a
count of independent algorithmic bugs.

All **149 polytools tests** are included in the final
[318-test focused check](../2026-11-09-wave46/raw/focused.log). The
[complete final suite](../2026-11-09-wave46/raw/suite-run3-summary.json) passed
**61,384 tests**, with **14 xfailed, 41 skipped and zero failures**. All 61,439
collected cases were executed exactly once, with none missing or unexpected.
The focused and full-suite
[identities](../2026-11-09-wave46/raw/suite-identity.json) record the same final
binary. [Cargo tests](../2026-11-09-wave46/raw/cargo-test.log) passed 879 unit
tests (one ignored), one integration test and five doctests;
[correctness clippy](../2026-11-09-wave46/raw/clippy.log) succeeded for core and
bindings. The [combined audit](../2026-11-09-wave46/README.md) distinguishes the
earlier intermediate checks and the final-source validation.

No speedup or asymptotic-complexity claim is made for this correctness lot.

## Remaining gaps

- `cancel(deep=True)`, Piecewise cancellation, `same_root`, private ring APIs,
  QQ<a> representatives and interruption/resource control are outside this lot.
- All combinations of list construction options (`expand=False`, non-greedy
  inference, generator ordering and composite domains) have not been audited.
- With a list of native Poly inputs and no generators/options, SymPy 1.14 raises
  an `AttributeError` in its numerical-domain detection. Symox supplies a valid
  Poly result; this difference is recorded, not counted as a reference match.
  Supplying `polys=True`/`False` or a domain exercises the compared paths.
- The reference itself loses a singleton Poly's coefficient domain when it
  reconstructs the list. This lot follows the reconstructed domain; it preserves
  the complete input domains for exactly two Poly inputs whose explicit
  generators, if given, have the same elements as the original generators.
- Public numerical expressions that are not represented as ordinary numeric
  atoms (for example general exact algebraic combinations) retain their existing
  paths. Selected examples do not prove parity of every numerical ground.
