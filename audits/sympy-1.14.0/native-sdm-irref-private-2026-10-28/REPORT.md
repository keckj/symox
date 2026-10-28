# Private sparse QQ Gauss–Jordan candidate

This candidate ports the complete selected `sdm_irref` field algorithm to Rust. It is based on accepted R6r1 (`80c98ada`, native `cc6d0210`), contains exactly five source/test changes and has not been integrated into main. It has no release benchmark or complete-suite certification. `final-validation.json` identifies the final source, binary and every accepted private gate.

The Rust kernel performs stable first-column sorting followed by reverse row processing, previous-pivot cancellation, exact QQ inversion/normalization, back elimination and final pivot publication. It preserves all three results, including nonzero-column insertion order and historical empty sets. Original pivot/column key objects and large row-number aliases are retained; the input containers and scalar identities remain unchanged. Sparse dictionaries never become dimension-sized dense storage.

Persistent CPython sets reproduce the established representation and iteration sequence. They perform no scalar algebra. Every membership result is checked against independent Rust set state, and temporary set slots are reused. The Python helper performs protocol selection only. Unsupported domains, explicit stored zeros, empty/malformed rows, custom containers/scalars, unsupported coordinate representations and consumed callbacks retain the complete prior algorithm. Exceptions and interruptions are propagated rather than reported as zero rank. The sparse field algorithm remains distinct from fraction-free reduction.

## Final validation

- Local: **226 focused tests**, then **4,689 tests in 50 coupled files**, all passing.
- Demeter: **226 focused tests**, then the same **4,689 tests in 50 coupled files**, all passing.
- Collection/finish Counters prove every collected node ran exactly once in each gate. All five source hashes and the identical DEBUG binary are checked before and after each gate.
- Five changed Rust unit tests pass locally. Workspace `cargo clippy --all-targets -- -D clippy::correctness` exits zero; all 335 native inputs remain unchanged. Existing warnings remain visible. These are local Rust gates; complete Rust validation on Demeter remains part of integration.
- A fresh process rejects all SymPy imports, runs the selected public routine and checks the entire result against independent standard-library Fraction elimination, with unchanged input scalar identities.

The 226 focused tests include 120 seeded comparisons against both the unchanged old path and SymPy 1.14, with exact row/pivot/rank/row-space certificates and ordered metadata/alias checks; 30 selected public AUTO/GJ checks; large sparse coordinates and 270-row metadata identity checks; callbacks, import order, protocol overrides and fallbacks.

## Bootstrap correction found during final review

The initial selector could capture a Python replacement of `builtins.min` before importing Symox and treat it as stock. The selected path then omitted four observable old-path callbacks. The final selector authenticates initial C builtin function ownership/type/name and static CPython factory types before trusting identity snapshots. Heap factory types with copied builtin names/modules decline. Custom `all`/`any` helpers captured by the selector itself also decline without invoking those callbacks; the general algorithm does not consume them.

Nine new regressions all fail on the archived initial helper and pass in the final focused and coupled gates. The `len` case tests initial registration independently because replacing `len` before importing the existing package prevents dataclasses from importing; that preexisting package import failure is preserved in `initial-new-regressions-focused-validated.log`. No assertions in the original 217 tests were removed or weakened. Earlier 217-test and 224-test cohorts are archived separately and do not stand in for the final gates.

## Delivery

`freeze-manifest.json` contains base/final hashes for each file and identifies the DEBUG binary. `frozen-source/` contains immutable copies of the five final files. `candidate.patch` passes `git apply --check` against the accepted baseline. Apply only after checking the five current base hashes, rebuild the combined release and rerun the complete Python/Rust suites on Demeter. New domain floor/ceil or LLL code has not been included in this private base; the integrated combination still needs its own gate. Other sparse domains and the general Python algorithms remain available; this candidate does not establish global mathematical parity.
