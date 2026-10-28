# Certified native CC products and quotients — frozen patch

The reported CC140/211 differences exposed real native rounding defects. They are corrected in a private patch using exact dyadic rounding certificates, without changing main, the frozen inexact_linear delivery, or any other received package. This is a numerical correctness repair, never a claimed speed gain over an old wrong result.

Stage: `/home/keckj/.cache/symox-validation/native-complex-rounding-20261027/work`.
Patch: `complex-rounding.patch` beside this report. Initial diagnosis is preserved in `DIAGNOSIS.md`; full before/after exact traces and reproducer accompany it.

## Source and binary identities

The base is a new independent copy of the frozen inexact_linear delivery, whose underlying 328 native build inputs match immutable r1b db859db2. Its complete 2891-file source manifest was verified before editing. No transitory main extension was used as a baseline.

- Base manifest: `62d8dbd2086febd4708dca0fd8d18baf1c4809045feff11561ef4616621541c5`.
- Final manifest (2893 files): `71aa6054ad10fc3261033189fa93bfbc7b0cd77193a35b5cb8b5656591a012be`.
- Patch: `12a9e4e69ef6ab3c7bfc00649859a007ee1ca8abb3034ca22e946470d359365b`.
- Final private DEBUG native: `527baf7753f13fc7ed576bfdb23392294757669ce6142b95df957690a592311b`.

`identity.json` records base/final hashes of the three changed files. Every final manifest entry was checked on Demeter after validation (`complex-rounding-source-check.log`), and `git apply --check` passes on base.

Three changed files only:

1. `crates/symox-core/src/poly/domains.rs`: private module declaration and dispatch of plain CC `*` and `/` after existing compatibility/pole checks.
2. New `crates/symox-core/src/poly/domains/complex_rounding.rs`: native certificates, exact fallback, five Rust tests.
3. New `tests/test_complex_domain_division_rounding.py`: 364 correctness/regression tests, including an independent integer rounding oracle.

There are no new Python production algorithms or binding changes. Addition/subtraction, symbolic expression simplification, numeric conversion, RR and parameter-domain algorithms retain their existing paths. Apply only these hunks after the independently delivered inexact_linear patch; do not replace unrelated current files.

## First faulty operations, measured against the exact same stored inputs

`audit.py` reads native `.as_expr()` and reference `._mpc_` mantissa/exponent payloads and separately verifies equality to each domain's `to_sympy` payload. Reference inputs are rebuilt from those exact dyades. There is no string-to-number conversion or tolerance. For each primitive, the exact rational complex result is computed using Fraction; bit lengths, integer comparison, integer division and the exact remainder implement independent nearest/ties-to-even rounding at p.

| Precision | First old native faulty primitive | Exact error from correctly rounded component |
| --- | --- | --- |
| 140 | second pivot reciprocal, imaginary component | `-2**(-140)`, exactly -1 ulp |
| 211 | first pivot reciprocal, real component | `-2**(-209)`, exactly -1 ulp |

The operands of those first errors were identical across implementations. The reference reciprocals were correctly rounded. Exact expected versus old native mantissas/exponents and complete rational inputs are in `DIAGNOSIS.md`, `complex-rounding-audit-v2.json` and `classification.json`.

The old field division constructed an expression using `ops::gaussian_reciprocal`: it rounded `r*r+i*i` and that norm's inverse at p, then multiplied the components. Exact integer simulation of those steps reproduces the wrong bits. Removing only eval_mp's p+4 final double rounding would not fix these reciprocal errors.

After fixing division, further incorrect native products became visible on the now-correct inputs. At p140 one imaginary product was +1 ulp wrong; at p211 two products had a -1 ulp component error. Exact rounding through p+4 and then p reproduces these errors (`complex-rounding-audit-after-v1.json`). Field multiplication had also been routed through expression evaluation. The patch therefore corrects both proven faulty primitives, multiplication and division; it does not broaden the change to addition/subtraction.

The reference p140 product `pivot1.eliminate[0,1].product` is itself +1 ulp wrong. Direct mpmath `mpc_mul` reproduces that output, while independent exact-rational rounding gives the other neighbor. The old native product also fails when given those particular reference operands, although its earlier erroneous reciprocal had previously caused it to take different operands there. The final native product is correctly rounded and intentionally differs from that reference bug. The new test retains this precise witness.

## Native proof mechanism and large exponents

Products and quotient numerators/denominators remain sums of exact dyadic terms, obtained by multiplying integer mantissas and adding arbitrary-precision exponents. No sum is flattened across an enormous exponent gap.

To determine a sparse sum's sign, the kernel merges the two terms with the greatest leading bit positions only when their magnitudes are close. Their alignment shifts are consequently bounded by already stored mantissa widths and the small term count. Otherwise the leading term exceeds twice an exact upper bound for the remaining tail, certifying its sign without materializing the gap. This also handles exact cancellation before a tiny remaining term. The same mechanism gives an exact magnitude using comparisons with powers of two.

A quick native approximation is only a candidate. Its lower and upper rounding-cell midpoints are constructed exactly, with the different spacing below a binade boundary and the correct tie parity, including p=1. Each component is accepted only after exact sparse comparisons of `N-midpoint*D`. Agreement between two approximate precisions is never used as a proof.

If the candidate fails, exact comparisons determine the quotient's binary magnitude, then binary search its bounded p-bit mantissa and apply the exact midpoint/tie test. There is no fixed guard budget, exponent-gap quota, fabricated zero or Unsupported-as-success path. A helper decline retains the complete previous expression path; incompatible domains and poles retain their earlier checks/errors. The tests exercise deliberately wrong/zero candidates to cover the exact search, and exponents +/-10**12 and +/-2**190 without dense shifts.

The selected domain p remains authoritative. Every result is built from components rounded to p and tagged at p; no presentation precision replaces calculation precision. Each public operator still rounds separately, so `a-b*c` retains its product rounding before subtraction. Python scalar/protocol calls remain in their original bindings and installed wrappers. Negative powers continue their existing operator sequence and receive the corrected final division.

## Demeter validation

Private tree: `/home/keckj/symox-ci/wt/complex-rounding-20261027`.
Private target: `/home/keckj/symox-ci/targets/complex-rounding`.
Cargo used the dedicated Python3.10 venv through PYO3_PYTHON; only the resulting .so was copied to this private package. No maturin/shared installation, main mutation, timing or benchmark.

- `cargo build -p symox-py`: pass (`complex-rounding-build-final.log`).
- `cargo test -p symox complex_rounding`: **5 passed** (`complex-rounding-rust-final.log`).
- Final focused command: **973 passed**, including **364 new** and **609 existing** (`complex-rounding-python-final.log`). It runs the new file, all 203 frozen inexact_linear tests, inexact followup, complex multivariate cancellation, composite arithmetic, domain division, and both Heurisch regression files. No existing test was modified.
- Complete final exact trace: **20/20 native primitives correctly rounded at p140 and 20/20 at p211**, independently certified against each operation's actual stored dyadic inputs (`complex-rounding-audit-after-final.json`). At p211 the full trace now matches the reference bitwise. At p140 the first divergence is exactly the reference's incorrect product, with native correctly rounded.

The new tests cover p1/53/139/140/211, additional tie precisions2/3, random exact-stored inputs, reciprocal/negative-power protocols, real/imaginary operands, zeros, near poles, cancellation, underflow-resistant tiny values, exponent gaps, all output components/declared precisions/domains, existing scalar callback paths, poles/errors and parameter-domain reconstruction. Addition/subtraction and separate operator grouping are checked as regressions. Full integration/release suites remain the parent's responsibility.

## Deliberate limits and remaining gaps

This patch proves rounded products/quotients for covered finite literal plain CC elements, preserving the full fallback outside that representation. It does not establish universal CC/SymPy parity or replace parameter-field algorithms.

General expression evalf/construction retains its established rounding conventions, including p+4 double rounding; the repair is specifically typed field arithmetic. The general legacy BigFloat/MpComplex algorithms are not globally replaced: their proposed candidate must pass exact certificates here. An old symbolic expression result is not silently redefined to match a field result.

A separate pre-existing representation-equality gap remains: in CC[t], `(a/one)==a` is false for `a=(1+I)*t+2` at p53/140/211 although the selected domains agree and the expression difference simplifies exactly to zero. The same probe is false in the frozen baseline and final code; reconstruction/domain tests retain the complete arithmetic, without claiming that equality API is repaired. The initial overly strong new equality assertion was replaced by the exact reconstruction certificate and this explicit gap report.

The two previously documented Heurisch unknown outcomes remain unknown; numerical precision repairs are not counted as newfound solution completeness. No claimed speedup derives from returning a formerly wrong or unsupported result.
