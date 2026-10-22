# Versioned SymPy parity audit

Reference: **SymPy 1.14.0**. Global algorithmic parity is **not certified**.

The catalog examines every documented module export and public callable class method.
A binding probe checks Python argument binding, including keyword names and positional options.
A passing binding probe does not establish that an option is implemented. A rejected probe
can also expose an overly permissive reference signature whose implementation rejects that call.

Documentation examples execute the reference and native libraries in separate processes,
with timeouts and an immutable native package snapshot. Reference imports and desktop viewers
are blocked in native workers. Output differences require review: different printing is not
automatically a mathematical error. Equal exceptions are recorded separately from successful examples.

## Evidence

- Module paths: 310.
- Export/method bindings, including aliases: 67547.
- Distinct reference callables: 6077.
- Documented objects replayed: 2916.
- Matching example observations: 1749 objects.
- Different observations to review: 1163 objects.
- Incomplete replays: 2 objects.
- Named-call objects replayed: 2916.
- Matching named-call observations: 1489 objects.
- Different named-call observations to review: 1425 objects.
- Incomplete named-call replays: 0 objects.

Reference Python source fingerprint: `972084bcdf2aeda02be0f83a39ce771094bd7657efa5a0e81ac4d75cad8737a9`.
Installed source files matching distribution RECORD: 1533; mismatches: 0.
The separate [published-wheel integrity record](../audits/sympy-1.14.0/reference-wheel-integrity.json)
compares these sources with the [SymPy 1.14.0 PyPI release](https://pypi.org/project/sympy/1.14.0/#files), without modifying the reference installation.

Catalog native fingerprint: `c6198033b800a80410459efc0220add61ad8a6868a9376220e56666c485b90d0`.

Example snapshot fingerprint: `dfb31b50e5d39540a91f64fee9bfe81fb2e4aa97e793be3cac10d7b86a7244b6`.

Named-call snapshot fingerprint: `mixed; consult each record`.

Named-call replays use actual reference argument binding to move eligible positional inputs
to their public keyword names, preserving positional-only and variadic calls. They record
69903 reference calls and 68692 native calls, with
69785 and 67922 successful calls respectively.
Cascading exceptions and repeated calls affect these counts; they are not a compliance percentage.

These fingerprints may differ while fixes are being made. Each example record keeps its
own fingerprint; stale evidence must be rerun before approving a changed implementation.
Targeted supplementary runs replace older observations by reference callable identity,
including aliases. Other records retain their original snapshots and stale status.
Ordinary evidence files: `audits/sympy-1.14.0/wave4-examples.json.gz`.
Named-call evidence files: `audits/sympy-1.14.0/fifth-keyword-examples.json.gz`, `audits/sympy-1.14.0/batch29-final-keywords.json.gz`, `audits/sympy-1.14.0/extra-isolation-keywords.json.gz`, `audits/sympy-1.14.0/native-protocol-keywords.json.gz`.
Scoped review source hashes are also checked against the working tree, so edited Rust
cannot retain a fresh review merely because an older compiled extension is still installed.

Scoped source-equivalence checks for container algorithms are recorded in
[source-certificates.json](../audits/sympy-1.14.0/source-certificates.json).
These compare signatures, normalized control flow and closed standard-library dependencies,
including error propagation and time/space behavior. Symbolic input-object protocols are a separate obligation.

Scoped invariant, exact-arithmetic and complexity reviews are recorded in
[permutation-review.json](../audits/sympy-1.14.0/permutation-review.json) and
[group-review.json](../audits/sympy-1.14.0/group-review.json).
Complete elementary solution-lattice invariants are recorded separately in
[solver-review.json](../audits/sympy-1.14.0/solver-review.json).
Radical controls and principal complex-branch invariants are recorded in
[radical-review.json](../audits/sympy-1.14.0/radical-review.json).
Scoped gamma-moment convergence and public Piecewise branch protocols are recorded in
[integration-review.json](../audits/sympy-1.14.0/integration-review.json).
Their finite-domain assumptions and unresolved complexity obligations are retained in each record;
they do not certify an entire constructor, inherited protocol or unsupported optional interface.

## Module traversal

“Instance” means that a method could not be inspected on its native constructor; it needs
an instance probe. This is not a claim that the method is absent on returned objects.

| Module | Exports | Method bindings | Binding gaps | Instance probes |
| --- | ---: | ---: | ---: | ---: |
| `sympy` | 921 | 28719 | 2285 | 21213 |
| `sympy.abc` | 0 | 0 | 0 | 0 |
| `sympy.algebras` | 1 | 0 | 0 | 0 |
| `sympy.assumptions` | 12 | 65 | 2 | 62 |
| `sympy.assumptions.ask` | 6 | 0 | 0 | 0 |
| `sympy.assumptions.assume` | 6 | 37 | 2 | 36 |
| `sympy.assumptions.predicates` | 0 | 0 | 0 | 0 |
| `sympy.assumptions.refine` | 9 | 0 | 0 | 0 |
| `sympy.calculus` | 17 | 0 | 0 | 0 |
| `sympy.calculus.euler` | 1 | 0 | 0 | 0 |
| `sympy.calculus.finite_diff` | 3 | 0 | 0 | 0 |
| `sympy.calculus.singularities` | 7 | 0 | 0 | 0 |
| `sympy.calculus.util` | 9 | 0 | 0 | 0 |
| `sympy.categories` | 11 | 276 | 2 | 265 |
| `sympy.categories.diagram_drawing` | 5 | 0 | 0 | 0 |
| `sympy.codegen` | 13 | 423 | 3 | 380 |
| `sympy.codegen.algorithms` | 2 | 0 | 0 | 0 |
| `sympy.codegen.approximations` | 2 | 6 | 0 | 6 |
| `sympy.codegen.ast` | 43 | 888 | 5 | 778 |
| `sympy.codegen.cfunctions` | 11 | 952 | 311 | 24 |
| `sympy.codegen.cnodes` | 11 | 259 | 0 | 232 |
| `sympy.codegen.cutils` | 1 | 0 | 0 | 0 |
| `sympy.codegen.cxxnodes` | 1 | 29 | 0 | 26 |
| `sympy.codegen.fnodes` | 32 | 1428 | 0 | 1388 |
| `sympy.codegen.futils` | 1 | 0 | 1 | 0 |
| `sympy.codegen.matrix_nodes` | 1 | 106 | 1 | 102 |
| `sympy.codegen.pyutils` | 1 | 0 | 1 | 0 |
| `sympy.codegen.rewriting` | 18 | 4 | 0 | 4 |
| `sympy.combinatorics` | 33 | 421 | 1 | 296 |
| `sympy.combinatorics.galois` | 21 | 0 | 0 | 0 |
| `sympy.combinatorics.generators` | 6 | 0 | 0 | 0 |
| `sympy.combinatorics.graycode` | 6 | 0 | 0 | 0 |
| `sympy.combinatorics.group_constructs` | 1 | 0 | 0 | 0 |
| `sympy.combinatorics.group_numbers` | 4 | 0 | 0 | 0 |
| `sympy.combinatorics.named_groups` | 8 | 0 | 0 | 0 |
| `sympy.combinatorics.partitions` | 7 | 0 | 0 | 0 |
| `sympy.combinatorics.perm_groups` | 4 | 0 | 0 | 0 |
| `sympy.combinatorics.permutations` | 4 | 95 | 0 | 93 |
| `sympy.combinatorics.polyhedron` | 1 | 0 | 0 | 0 |
| `sympy.combinatorics.prufer` | 1 | 0 | 0 | 0 |
| `sympy.combinatorics.subsets` | 2 | 0 | 0 | 0 |
| `sympy.combinatorics.tensor_can` | 12 | 0 | 0 | 0 |
| `sympy.combinatorics.testutil` | 2 | 0 | 0 | 0 |
| `sympy.combinatorics.util` | 0 | 0 | 0 | 0 |
| `sympy.core` | 101 | 0 | 13 | 0 |
| `sympy.core.add` | 1 | 0 | 1 | 0 |
| `sympy.core.assumptions` | 8 | 13 | 1 | 0 |
| `sympy.core.basic` | 3 | 0 | 1 | 0 |
| `sympy.core.cache` | 5 | 0 | 0 | 0 |
| `sympy.core.containers` | 5 | 9 | 0 | 2 |
| `sympy.core.evalf` | 44 | 2 | 7 | 0 |
| `sympy.core.expr` | 6 | 4 | 2 | 0 |
| `sympy.core.exprtools` | 7 | 16 | 0 | 0 |
| `sympy.core.function` | 28 | 229 | 77 | 5 |
| `sympy.core.intfunc` | 10 | 0 | 0 | 0 |
| `sympy.core.kind` | 3 | 3 | 1 | 1 |
| `sympy.core.mod` | 1 | 0 | 0 | 0 |
| `sympy.core.mul` | 4 | 0 | 1 | 0 |
| `sympy.core.multidimensional` | 4 | 0 | 0 | 0 |
| `sympy.core.numbers` | 34 | 1692 | 736 | 71 |
| `sympy.core.power` | 1 | 0 | 0 | 0 |
| `sympy.core.random` | 4 | 0 | 0 | 0 |
| `sympy.core.relational` | 14 | 0 | 0 | 0 |
| `sympy.core.singleton` | 3 | 2 | 3 | 0 |
| `sympy.core.symbol` | 8 | 28 | 2 | 26 |
| `sympy.core.sympify` | 4 | 0 | 0 | 0 |
| `sympy.core.traversal` | 7 | 0 | 0 | 0 |
| `sympy.crypto` | 56 | 0 | 0 | 0 |
| `sympy.crypto.crypto` | 57 | 2 | 0 | 0 |
| `sympy.diffgeom` | 24 | 1082 | 15 | 1021 |
| `sympy.discrete` | 11 | 0 | 0 | 0 |
| `sympy.discrete.convolutions` | 8 | 0 | 0 | 0 |
| `sympy.discrete.transforms` | 8 | 0 | 0 | 0 |
| `sympy.functions` | 178 | 98 | 3 | 98 |
| `sympy.functions.elementary` | 0 | 0 | 0 | 0 |
| `sympy.functions.elementary.complexes` | 14 | 0 | 0 | 0 |
| `sympy.functions.special` | 0 | 0 | 0 | 0 |
| `sympy.functions.special.bessel` | 21 | 392 | 0 | 388 |
| `sympy.functions.special.beta_functions` | 4 | 0 | 2 | 0 |
| `sympy.functions.special.elliptic_integrals` | 4 | 0 | 0 | 0 |
| `sympy.functions.special.error_functions` | 21 | 196 | 66 | 4 |
| `sympy.functions.special.gamma_functions` | 9 | 0 | 0 | 0 |
| `sympy.functions.special.hyper` | 17 | 1306 | 3 | 1300 |
| `sympy.functions.special.mathieu_functions` | 5 | 98 | 0 | 97 |
| `sympy.functions.special.polynomials` | 14 | 98 | 32 | 2 |
| `sympy.functions.special.singularity_functions` | 1 | 0 | 1 | 0 |
| `sympy.functions.special.zeta_functions` | 6 | 0 | 0 | 0 |
| `sympy.geometry` | 30 | 0 | 10 | 0 |
| `sympy.geometry.curve` | 1 | 0 | 0 | 0 |
| `sympy.geometry.ellipse` | 2 | 0 | 1 | 0 |
| `sympy.geometry.entity` | 5 | 92 | 5 | 61 |
| `sympy.geometry.line` | 12 | 195 | 26 | 93 |
| `sympy.geometry.plane` | 1 | 0 | 0 | 0 |
| `sympy.geometry.point` | 3 | 0 | 3 | 0 |
| `sympy.geometry.polygon` | 5 | 0 | 0 | 0 |
| `sympy.geometry.util` | 9 | 0 | 0 | 0 |
| `sympy.holonomic` | 9 | 16 | 6 | 0 |
| `sympy.integrals` | 30 | 0 | 3 | 0 |
| `sympy.integrals.intpoly` | 22 | 0 | 0 | 0 |
| `sympy.integrals.meijerint` | 3 | 0 | 0 | 0 |
| `sympy.integrals.meijerint_doc` | 0 | 0 | 0 | 0 |
| `sympy.integrals.transforms` | 32 | 398 | 2 | 396 |
| `sympy.interactive` | 3 | 0 | 0 | 0 |
| `sympy.interactive.printing` | 1 | 0 | 0 | 0 |
| `sympy.interactive.session` | 6 | 0 | 0 | 0 |
| `sympy.liealgebras` | 1 | 0 | 1 | 0 |
| `sympy.logic` | 20 | 0 | 1 | 0 |
| `sympy.logic.inference` | 7 | 6 | 1 | 0 |
| `sympy.matrices` | 73 | 45 | 27 | 26 |
| `sympy.matrices.expressions` | 40 | 105 | 8 | 103 |
| `sympy.matrices.expressions.blockmatrix` | 16 | 0 | 1 | 0 |
| `sympy.matrices.immutable` | 6 | 249 | 3 | 249 |
| `sympy.matrices.kind` | 1 | 0 | 0 | 0 |
| `sympy.matrices.matrixbase` | 7 | 1 | 0 | 0 |
| `sympy.matrices.sparse` | 3 | 174 | 2 | 174 |
| `sympy.matrices.sparsetools` | 1 | 0 | 0 | 0 |
| `sympy.ntheory` | 73 | 0 | 8 | 0 |
| `sympy.ntheory.bbp_pi` | 1 | 0 | 0 | 0 |
| `sympy.ntheory.continued_fraction` | 5 | 0 | 0 | 0 |
| `sympy.ntheory.digits` | 3 | 0 | 0 | 0 |
| `sympy.ntheory.ecm` | 2 | 3 | 0 | 0 |
| `sympy.ntheory.egyptian_fraction` | 6 | 0 | 0 | 0 |
| `sympy.ntheory.factor_` | 38 | 12 | 6 | 0 |
| `sympy.ntheory.generate` | 11 | 0 | 0 | 0 |
| `sympy.ntheory.modular` | 5 | 0 | 0 | 0 |
| `sympy.ntheory.multinomial` | 4 | 0 | 0 | 0 |
| `sympy.ntheory.partitions_` | 1 | 0 | 0 | 0 |
| `sympy.ntheory.primetest` | 12 | 0 | 0 | 0 |
| `sympy.ntheory.qs` | 4 | 2 | 0 | 0 |
| `sympy.ntheory.residue_ntheory` | 16 | 0 | 3 | 0 |
| `sympy.parsing` | 1 | 0 | 0 | 0 |
| `sympy.parsing.sym_expr` | 1 | 5 | 0 | 0 |
| `sympy.parsing.sympy_parser` | 21 | 19 | 0 | 0 |
| `sympy.physics` | 5 | 0 | 0 | 0 |
| `sympy.physics.biomechanics` | 15 | 810 | 22 | 762 |
| `sympy.physics.biomechanics.activation` | 3 | 0 | 0 | 0 |
| `sympy.physics.biomechanics.curve` | 9 | 0 | 8 | 0 |
| `sympy.physics.biomechanics.musculotendon` | 3 | 0 | 0 | 0 |
| `sympy.physics.continuum_mechanics` | 4 | 60 | 5 | 24 |
| `sympy.physics.continuum_mechanics.arch` | 1 | 0 | 0 | 0 |
| `sympy.physics.continuum_mechanics.beam` | 2 | 38 | 1 | 29 |
| `sympy.physics.continuum_mechanics.cable` | 1 | 0 | 0 | 0 |
| `sympy.physics.continuum_mechanics.truss` | 1 | 0 | 0 | 0 |
| `sympy.physics.control` | 33 | 337 | 13 | 278 |
| `sympy.physics.control.lti` | 16 | 0 | 0 | 0 |
| `sympy.physics.hep` | 0 | 0 | 0 | 0 |
| `sympy.physics.hep.gamma_matrices` | 5 | 0 | 0 | 0 |
| `sympy.physics.hydrogen` | 4 | 0 | 0 | 0 |
| `sympy.physics.matrices` | 4 | 0 | 0 | 0 |
| `sympy.physics.mechanics` | 74 | 239 | 24 | 114 |
| `sympy.physics.mechanics.actuator` | 7 | 0 | 0 | 0 |
| `sympy.physics.mechanics.joint` | 7 | 0 | 1 | 0 |
| `sympy.physics.mechanics.kane` | 1 | 0 | 0 | 0 |
| `sympy.physics.mechanics.lagrange` | 1 | 0 | 0 | 0 |
| `sympy.physics.mechanics.linearize` | 1 | 0 | 0 | 0 |
| `sympy.physics.mechanics.pathway` | 4 | 0 | 1 | 0 |
| `sympy.physics.mechanics.wrapping_geometry` | 3 | 0 | 1 | 0 |
| `sympy.physics.optics` | 39 | 1640 | 14 | 1612 |
| `sympy.physics.optics.gaussopt` | 16 | 0 | 0 | 0 |
| `sympy.physics.optics.medium` | 1 | 0 | 0 | 0 |
| `sympy.physics.optics.polarization` | 11 | 0 | 0 | 0 |
| `sympy.physics.optics.utils` | 10 | 0 | 0 | 0 |
| `sympy.physics.optics.waves` | 1 | 0 | 0 | 0 |
| `sympy.physics.paulialgebra` | 1 | 0 | 0 | 0 |
| `sympy.physics.qho_1d` | 3 | 0 | 0 | 0 |
| `sympy.physics.quantum` | 44 | 2636 | 51 | 2410 |
| `sympy.physics.quantum.anticommutator` | 1 | 0 | 1 | 0 |
| `sympy.physics.quantum.cartesian` | 15 | 1075 | 29 | 967 |
| `sympy.physics.quantum.cg` | 5 | 380 | 4 | 372 |
| `sympy.physics.quantum.circuitplot` | 7 | 210 | 6 | 186 |
| `sympy.physics.quantum.commutator` | 1 | 0 | 1 | 0 |
| `sympy.physics.quantum.constants` | 2 | 0 | 0 | 0 |
| `sympy.physics.quantum.dagger` | 1 | 0 | 1 | 0 |
| `sympy.physics.quantum.gate` | 29 | 1523 | 54 | 1343 |
| `sympy.physics.quantum.grover` | 5 | 202 | 6 | 178 |
| `sympy.physics.quantum.hilbert` | 8 | 0 | 1 | 0 |
| `sympy.physics.quantum.innerproduct` | 1 | 0 | 0 | 0 |
| `sympy.physics.quantum.operator` | 6 | 0 | 3 | 0 |
| `sympy.physics.quantum.operatorset` | 2 | 0 | 0 | 0 |
| `sympy.physics.quantum.piab` | 3 | 293 | 8 | 261 |
| `sympy.physics.quantum.qapply` | 1 | 0 | 0 | 0 |
| `sympy.physics.quantum.qft` | 4 | 305 | 10 | 269 |
| `sympy.physics.quantum.qubit` | 11 | 394 | 12 | 350 |
| `sympy.physics.quantum.represent` | 6 | 0 | 0 | 0 |
| `sympy.physics.quantum.shor` | 6 | 102 | 3 | 89 |
| `sympy.physics.quantum.spin` | 25 | 1570 | 38 | 1410 |
| `sympy.physics.quantum.state` | 13 | 0 | 1 | 0 |
| `sympy.physics.quantum.tensorproduct` | 2 | 0 | 0 | 0 |
| `sympy.physics.secondquant` | 37 | 1669 | 42 | 1519 |
| `sympy.physics.sho` | 2 | 0 | 0 | 0 |
| `sympy.physics.units` | 361 | 243 | 5 | 223 |
| `sympy.physics.units.dimensions` | 2 | 0 | 0 | 0 |
| `sympy.physics.units.prefixes` | 2 | 95 | 0 | 95 |
| `sympy.physics.units.quantities` | 2 | 98 | 0 | 92 |
| `sympy.physics.units.unitsystem` | 1 | 0 | 0 | 0 |
| `sympy.physics.units.util` | 3 | 0 | 0 | 0 |
| `sympy.physics.vector` | 27 | 0 | 2 | 0 |
| `sympy.physics.vector.functions` | 14 | 0 | 0 | 0 |
| `sympy.physics.vector.point` | 1 | 0 | 0 | 0 |
| `sympy.physics.wigner` | 11 | 98 | 0 | 96 |
| `sympy.plotting` | 11 | 3 | 0 | 0 |
| `sympy.plotting.plot` | 12 | 17 | 4 | 5 |
| `sympy.plotting.pygletplot` | 1 | 0 | 0 | 0 |
| `sympy.polys` | 176 | 0 | 9 | 0 |
| `sympy.polys.domains` | 32 | 0 | 4 | 0 |
| `sympy.polys.fields` | 4 | 0 | 0 | 0 |
| `sympy.polys.matrices` | 2 | 91 | 9 | 46 |
| `sympy.polys.matrices._dfm` | 1 | 57 | 10 | 4 |
| `sympy.polys.matrices._typing` | 1 | 0 | 0 | 0 |
| `sympy.polys.matrices.ddm` | 1 | 71 | 29 | 4 |
| `sympy.polys.matrices.dense` | 17 | 0 | 0 | 0 |
| `sympy.polys.matrices.domainmatrix` | 2 | 0 | 0 | 0 |
| `sympy.polys.matrices.sdm` | 13 | 70 | 42 | 4 |
| `sympy.polys.numberfields` | 10 | 0 | 0 | 0 |
| `sympy.polys.numberfields.modules` | 13 | 94 | 6 | 78 |
| `sympy.polys.numberfields.subfield` | 3 | 0 | 0 | 0 |
| `sympy.polys.polyconfig` | 4 | 0 | 0 | 0 |
| `sympy.polys.polyoptions` | 2 | 0 | 0 | 0 |
| `sympy.polys.rings` | 4 | 0 | 0 | 0 |
| `sympy.polys.solvers` | 5 | 179 | 2 | 175 |
| `sympy.printing` | 41 | 0 | 0 | 0 |
| `sympy.printing.aesaracode` | 4 | 3 | 1 | 1 |
| `sympy.printing.c` | 10 | 23 | 0 | 11 |
| `sympy.printing.codeprinter` | 11 | 2 | 0 | 0 |
| `sympy.printing.conventions` | 2 | 0 | 0 | 0 |
| `sympy.printing.cxx` | 3 | 18 | 0 | 9 |
| `sympy.printing.fortran` | 3 | 6 | 0 | 3 |
| `sympy.printing.gtk` | 1 | 0 | 0 | 0 |
| `sympy.printing.jscode` | 4 | 6 | 0 | 3 |
| `sympy.printing.julia` | 5 | 6 | 0 | 3 |
| `sympy.printing.lambdarepr` | 6 | 25 | 0 | 10 |
| `sympy.printing.latex` | 7 | 5 | 0 | 3 |
| `sympy.printing.maple` | 3 | 5 | 0 | 2 |
| `sympy.printing.mathematica` | 3 | 5 | 0 | 2 |
| `sympy.printing.mathml` | 6 | 16 | 0 | 8 |
| `sympy.printing.octave` | 5 | 6 | 0 | 3 |
| `sympy.printing.precedence` | 12 | 0 | 0 | 0 |
| `sympy.printing.pretty` | 6 | 0 | 0 | 0 |
| `sympy.printing.pretty.pretty` | 19 | 4 | 0 | 2 |
| `sympy.printing.pretty.pretty_symbology` | 25 | 0 | 0 | 0 |
| `sympy.printing.pretty.stringpict` | 2 | 29 | 0 | 0 |
| `sympy.printing.preview` | 3 | 0 | 0 | 0 |
| `sympy.printing.printer` | 3 | 3 | 0 | 1 |
| `sympy.printing.pycode` | 7 | 15 | 0 | 6 |
| `sympy.printing.python` | 3 | 6 | 0 | 3 |
| `sympy.printing.rcode` | 4 | 6 | 0 | 3 |
| `sympy.printing.repr` | 2 | 4 | 0 | 2 |
| `sympy.printing.rust` | 5 | 297 | 3 | 294 |
| `sympy.printing.smtlib` | 2 | 3 | 0 | 1 |
| `sympy.printing.str` | 4 | 5 | 0 | 2 |
| `sympy.printing.tree` | 4 | 0 | 0 | 0 |
| `sympy.series` | 18 | 0 | 2 | 0 |
| `sympy.series.sequences` | 10 | 155 | 13 | 112 |
| `sympy.sets` | 26 | 0 | 7 | 0 |
| `sympy.sets.conditionset` | 3 | 0 | 1 | 0 |
| `sympy.sets.fancysets` | 12 | 371 | 114 | 215 |
| `sympy.sets.powerset` | 1 | 0 | 0 | 0 |
| `sympy.sets.sets` | 22 | 92 | 40 | 52 |
| `sympy.simplify` | 32 | 0 | 1 | 0 |
| `sympy.simplify.combsimp` | 1 | 0 | 0 | 0 |
| `sympy.simplify.cse_main` | 11 | 7 | 0 | 0 |
| `sympy.simplify.epathtools` | 2 | 0 | 0 | 0 |
| `sympy.simplify.fu` | 31 | 0 | 1 | 0 |
| `sympy.simplify.hyperexpand` | 37 | 214 | 3 | 188 |
| `sympy.simplify.hyperexpand_doc` | 0 | 0 | 0 | 0 |
| `sympy.simplify.powsimp` | 2 | 0 | 0 | 0 |
| `sympy.simplify.radsimp` | 17 | 0 | 0 | 0 |
| `sympy.simplify.ratsimp` | 2 | 0 | 0 | 0 |
| `sympy.simplify.simplify` | 21 | 0 | 0 | 0 |
| `sympy.simplify.sqrtdenest` | 6 | 1 | 0 | 0 |
| `sympy.simplify.trigsimp` | 5 | 0 | 0 | 0 |
| `sympy.solvers` | 46 | 0 | 2 | 0 |
| `sympy.solvers.inequalities` | 8 | 0 | 0 | 0 |
| `sympy.solvers.ode` | 11 | 0 | 1 | 0 |
| `sympy.solvers.ode.ode` | 23 | 0 | 1 | 0 |
| `sympy.solvers.pde` | 9 | 0 | 0 | 0 |
| `sympy.solvers.recurr` | 4 | 0 | 0 | 0 |
| `sympy.solvers.simplex` | 6 | 2 | 0 | 0 |
| `sympy.solvers.solvers` | 15 | 0 | 1 | 0 |
| `sympy.solvers.solveset` | 14 | 1 | 0 | 0 |
| `sympy.stats` | 157 | 1389 | 13 | 1335 |
| `sympy.stats.crv` | 9 | 264 | 7 | 230 |
| `sympy.stats.crv_types` | 54 | 0 | 0 | 0 |
| `sympy.stats.frv` | 9 | 283 | 8 | 236 |
| `sympy.stats.frv_types` | 11 | 0 | 0 | 0 |
| `sympy.stats.rv` | 41 | 633 | 8 | 610 |
| `sympy.tensor` | 19 | 0 | 3 | 0 |
| `sympy.tensor.array` | 16 | 62 | 2 | 50 |
| `sympy.tensor.array.expressions` | 16 | 597 | 8 | 567 |
| `sympy.tensor.index_methods` | 3 | 1 | 0 | 0 |
| `sympy.tensor.indexed` | 4 | 1 | 1 | 0 |
| `sympy.tensor.tensor` | 35 | 831 | 30 | 661 |
| `sympy.tensor.toperators` | 1 | 99 | 4 | 83 |
| `sympy.testing` | 2 | 0 | 0 | 0 |
| `sympy.testing.pytest` | 16 | 4 | 1 | 0 |
| `sympy.testing.randtest` | 0 | 0 | 0 | 0 |
| `sympy.testing.runtests` | 21 | 39 | 3 | 2 |
| `sympy.utilities` | 25 | 0 | 1 | 0 |
| `sympy.utilities.autowrap` | 9 | 9 | 5 | 3 |
| `sympy.utilities.codegen` | 16 | 32 | 0 | 16 |
| `sympy.utilities.decorator` | 10 | 0 | 0 | 0 |
| `sympy.utilities.enumerative` | 6 | 13 | 0 | 0 |
| `sympy.utilities.exceptions` | 3 | 1 | 0 | 0 |
| `sympy.utilities.iterables` | 59 | 0 | 1 | 0 |
| `sympy.utilities.lambdify` | 3 | 0 | 1 | 0 |
| `sympy.utilities.memoization` | 2 | 0 | 0 | 0 |
| `sympy.utilities.misc` | 13 | 1 | 0 | 0 |
| `sympy.utilities.source` | 2 | 0 | 0 | 0 |
| `sympy.utilities.timeutils` | 2 | 0 | 0 | 0 |
| `sympy.vector` | 43 | 1997 | 17 | 1802 |

## Function and constructor binding gaps

Aliases share their reference identity below. The compressed catalog retains every public binding,
all method signatures and concrete rejected-call witnesses. These rows are an open work queue.

| Reference callable | Reference signature | Native signature |
| --- | --- | --- |
| `sympy.Add` | `(*args, evaluate=None, _sympify=True)` | `(*args, evaluate=None)` |
| `sympy.Adjoint` | `(*args, **kwargs)` | `(arg)` |
| `sympy.AlgebraicNumber` | `(expr, coeffs=None, alias=None, **args)` | `(expr, coeffs=None, alias=None, *, gen=None)` |
| `sympy.Atom` | `(*args)` | `()` |
| `sympy.AtomicExpr` | `(*args)` | `()` |
| `sympy.BlockMatrix` | `(*args, **kwargs)` | `(blocks)` |
| `sympy.Circle` | `(*args, **kwargs)` | `(center, radius=None, third=None, **kwargs)` |
| `sympy.Complement` | `(a, b, evaluate=True)` | `(a, b)` |
| `sympy.Contains` | `(x, s, evaluate=None)` | `(x, s)` |
| `sympy.DenseNDimArray` | `(*args, **kwargs)` | `(iterable=None, shape=None, **options)` |
| `sympy.DiagonalMatrix` | `(*args, **kwargs)` | `(matrix)` |
| `sympy.DiagonalOf` | `(*args, **kwargs)` | `(matrix)` |
| `sympy.Expr` | `(*args)` | `(obj)` |
| `sympy.FF` | `(mod, symmetric=True)` | `(modulus, symmetric=True)` |
| `sympy.FiniteSet` | `(*args, **kwargs)` | `(*args)` |
| `sympy.FractionField` | `(domain_or_field, symbols=None, order=None)` | `(dom, *gens, order='lex')` |
| `sympy.Function` | `(*args)` | `(name: 'str', *, nargs=None, **assumptions)` |
| `sympy.GroebnerBasis` | `(F, *gens, **args)` | `(f, *gens, order='lex', domain=None, modulus=None, method=None, polys=None, extension=None)` |
| `sympy.ImmutableDenseMatrix` | `(*args, **kwargs)` | `(*args)` |
| `sympy.ImmutableSparseMatrix` | `(*args, **kwargs)` | `(*args)` |
| `sympy.Indexed` | `(base, *args, **kw_args)` | `(base, *indices)` |
| `sympy.Intersection` | `(*args, evaluate=None)` | `(*args)` |
| `sympy.InverseLaplaceTransform` | `(*args)` | `(F, s, x, plane, **opts)` |
| `sympy.InverseMellinTransform` | `(*args)` | `(F, s, x, a, b, **opts)` |
| `sympy.Line` | `(*args, **kwargs)` | `(p1, p2=None, **kwargs)` |
| `sympy.Line2D` | `(p1, pt=None, slope=None, **kwargs)` | `(p1, p2=None, **kwargs)` |
| `sympy.Line3D` | `(p1, pt=None, direction_ratio=(), **kwargs)` | `(p1, p2=None, **kwargs)` |
| `sympy.MatPow` | `(base, exp, evaluate=False, **options)` | `(base, exp, evaluate=False)` |
| `sympy.Matrix` | `(*args, **kwargs)` | `(*args)` |
| `sympy.MatrixExpr` | `(*args, **kwargs)` | `(value)` |
| `sympy.Mul` | `(*args, evaluate=None, _sympify=True)` | `(*args, evaluate=None)` |
| `sympy.MutableSparseMatrix` | `(*args, **kwargs)` | `(*args)` |
| `sympy.Not` | `(arg)` | `(a)` |
| `sympy.Number` | `(*obj)` | `()` |
| `sympy.Parabola` | `(focus=None, directrix=None, **kwargs)` | `(focus, directrix, **kwargs)` |
| `sympy.Point` | `(*args, **kwargs)` | `(*args, dim=None, evaluate=True, on_morph='ignore')` |
| `sympy.Point2D` | `(*args, _nocheck=False, **kwargs)` | `(*args, dim=None, evaluate=True, on_morph='ignore')` |
| `sympy.Point3D` | `(*args, _nocheck=False, **kwargs)` | `(*args, dim=None, evaluate=True, on_morph='ignore')` |
| `sympy.PolynomialRing` | `(domain_or_ring, symbols=None, order=None)` | `(dom, *gens, order='lex')` |
| `sympy.Predicate` | `(*args, **kwargs)` | `(name)` |
| `sympy.ProductSet` | `(*sets, **assumptions)` | `(*args)` |
| `sympy.Rational` | `(p, q=None, gcd=None)` | `(p, q=None)` |
| `sympy.Ray2D` | `(p1, pt=None, angle=None, **kwargs)` | `(p1, p2=None, **kwargs)` |
| `sympy.Ray3D` | `(p1, pt=None, direction_ratio=(), **kwargs)` | `(p1, p2=None, **kwargs)` |
| `sympy.S` | `(locals=None, convert_xor=True, strict=False, rational=False, evaluate=None)` | `(a, locals=None, convert_xor=True, strict=False, rational=False, evaluate=None)` |
| `sympy.SeqAdd` | `(*args, **kwargs)` | `(*arguments, evaluate=True)` |
| `sympy.SeqMul` | `(*args, **kwargs)` | `(*arguments, evaluate=True)` |
| `sympy.SingularityFunction` | `(variable, offset, exponent)` | `(x, a, n)` |
| `sympy.SparseNDimArray` | `(*args, **kwargs)` | `(iterable=None, shape=None, **options)` |
| `sympy.SymmetricDifference` | `(a, b, evaluate=True)` | `(a, b)` |
| `sympy.Transpose` | `(*args, **kwargs)` | `(arg)` |
| `sympy.Union` | `(*args, **kwargs)` | `(*args)` |
| `sympy.WildFunction` | `(*args)` | `(name, **assumptions)` |
| `sympy.ZeroMatrix` | `(m, n)` | `(n, m)` |
| `sympy.betainc` | `(*args)` | `(a, b, x1, x2)` |
| `sympy.betainc_regularized` | `(*args)` | `(a, b, x1, x2)` |
| `sympy.dsolve` | `(eq, func=None, hint='default', simplify=True, ics=None, xi=None, eta=None, x0=0, n=6, **kwargs)` | `(eq, func=None, hint='default', simplify=True, ics=None, **kwargs)` |
| `sympy.expand` | `(e, deep=True, modulus=None, power_base=True, power_exp=True, mul=True, log=True, multinomial=True, basic=True, **hints)` | `(e, deep=True, **hints)` |
| `sympy.expand_log` | `(expr, deep=True, force=False, factor=False)` | `(expr, deep=True, force=False)` |
| `sympy.eye` | `(*args, **kwargs)` | `(n, m=None, **_kwargs)` |
| `sympy.groebner` | `(F, *gens, **args)` | `(f, *gens, order='lex', domain=None, modulus=None, method=None, polys=None, extension=None)` |
| `sympy.hyperexpand` | `(f, allow_hyper=False, rewrite='default', place=None)` | `(e)` |
| `sympy.integrate` | `(*args, meijerg=None, conds='piecewise', risch=None, heurisch=None, manual=None, **kwargs)` | `(f, *args, **_kwargs)` |
| `sympy.is_zero_dimensional` | `(F, *gens, **args)` | `(F, *gens, order='lex', domain=None, modulus=None, method=None, extension=None)` |
| `sympy.lambdify` | `(args, expr, modules=None, printer=None, use_imps=True, dummify=False, cse=False, docstring_limit=1000)` | `(args, expr, modules=None, **kwargs)` |
| `sympy.nroots` | `(f, n=15, maxsteps=50, cleanup=True)` | `(f, n=15, **_kwargs)` |
| `sympy.nsolve` | `(*args, dict=False, **kwargs)` | `(f, x, x0=None, **_kwargs)` |
| `sympy.ones` | `(*args, **kwargs)` | `(r, c=None, **_kwargs)` |
| `sympy.product` | `(*args, **kwargs)` | `(f, *limits)` |
| `sympy.summation` | `(f, *symbols, **kwargs)` | `(f, *limits)` |
| `sympy.together` | `(expr, deep=False, fraction=True)` | `(f, deep=False, **_kwargs)` |
| `sympy.zeros` | `(*args, **kwargs)` | `(r, c=None, **_kwargs)` |
| `sympy.assumptions.assume.UndefinedPredicate` | `(name, handlers=None)` | `(name)` |
| `sympy.categories.Diagram` | `(*args)` | `(premises=(), conclusions=())` |
| `sympy.categories.Object` | `(name, **assumptions)` | `(name)` |
| `sympy.codegen.ast.Comment` | `(*args, **kwargs)` | `(text: str) -> None` |
| `sympy.codegen.cfunctions.exp2` | `(arg)` | `(*args, evaluate=True)` |
| `sympy.codegen.cfunctions.expm1` | `(arg)` | `(*args, evaluate=True)` |
| `sympy.codegen.cfunctions.isinf` | `(arg)` | `(*args, evaluate=True)` |
| `sympy.codegen.cfunctions.isnan` | `(arg)` | `(*args, evaluate=True)` |
| `sympy.codegen.cfunctions.log10` | `(arg)` | `(*args, evaluate=True)` |
| `sympy.codegen.cfunctions.log1p` | `(arg)` | `(*args, evaluate=True)` |
| `sympy.codegen.cfunctions.log2` | `(arg)` | `(*args, evaluate=True)` |
| `sympy.codegen.futils.render_as_module` | `(definitions, name, declarations=(), printer_settings=None)` | `(content, name, Printer=<class 'symox.printing.fortran.FCodePrinter'>, settings=None)` |
| `sympy.codegen.matrix_nodes.MatrixSolve` | `(*args, **kwargs)` | `(A, b)` |
| `sympy.codegen.pyutils.render_as_module` | `(content, standard='python3')` | `(content, Printer=<class 'symox.printing.pycode.PythonCodePrinter'>, settings=None)` |
| `sympy.core.assumptions.ManagedProperties` | `(*args, **kwargs)` | `(name, bases, namespace, **kwargs)` |
| `sympy.core.evalf.evalf_add` | `(v: "'Add'", prec: 'int', options: 'OPT_DICT') -> 'TMP_RES'` | `(expr, prec, options)` |
| `sympy.core.evalf.evalf_alg_num` | `(a: "'AlgebraicNumber'", prec: 'int', options: 'OPT_DICT') -> 'TMP_RES'` | `(expr, prec, options)` |
| `sympy.core.evalf.evalf_atan` | `(v: "'atan'", prec: 'int', options: 'OPT_DICT') -> 'TMP_RES'` | `(expr, prec, options)` |
| `sympy.core.evalf.evalf_mul` | `(v: "'Mul'", prec: 'int', options: 'OPT_DICT') -> 'TMP_RES'` | `(expr, prec, options)` |
| `sympy.core.evalf.evalf_pow` | `(v: "'Pow'", prec: 'int', options) -> 'TMP_RES'` | `(expr, prec, options)` |
| `sympy.core.evalf.evalf_symbol` | `(x: 'Expr', prec: 'int', options: 'OPT_DICT') -> 'TMP_RES'` | `(expr, prec, options)` |
| `sympy.core.evalf.evalf_trig` | `(v: 'Expr', prec: 'int', options: 'OPT_DICT') -> 'TMP_RES'` | `(expr, prec, options)` |
| `sympy.core.function.UndefSageHelper` | `()` | `(function)` |
| `sympy.core.kind.KindMeta` | `(clsname, bases, dct)` | `(name, bases, namespace, **kwargs)` |
| `sympy.core.numbers.RationalConstant` | `()` | `(p, q=None)` |
| `sympy.core.singleton.Singleton` | `(*args, **kwargs)` | `(name, bases, namespace, **kwargs)` |
| `sympy.core.symbol.Str` | `(name, **kwargs)` | `(name)` |
| `sympy.functions.special.error_functions.FresnelIntegral` | `(z)` | `(*args)` |
| `sympy.functions.special.error_functions.TrigonometricIntegral` | `(z)` | `(*args)` |
| `sympy.functions.special.hyper.TupleArg` | `(*args, **kwargs)` | `(*args)` |
| `sympy.liealgebras.CartanType` | `(*args)` | `(specification)` |
| `sympy.logic.inference.entails` | `(expr, formula_set=None)` | `(expr, formula_set)` |
| `sympy.matrices.dotprodsimp` | `(x)` | `(expr, withsimp=False)` |
| `sympy.ntheory.divisor_sigma` | `(n, k=1)` | `(*args)` |
| `sympy.ntheory.jacobi_symbol` | `(m, n)` | `(*args)` |
| `sympy.ntheory.legendre_symbol` | `(a, p)` | `(*args)` |
| `sympy.ntheory.mobius` | `(n)` | `(*args)` |
| `sympy.ntheory.primenu` | `(n)` | `(*args)` |
| `sympy.ntheory.primeomega` | `(n)` | `(*args)` |
| `sympy.ntheory.reduced_totient` | `(n)` | `(*args)` |
| `sympy.ntheory.totient` | `(n)` | `(*args)` |
| `sympy.ntheory.factor_.udivisor_sigma` | `(n, k=1)` | `(*args)` |
| `sympy.physics.biomechanics.CharacteristicCurveFunction` | `()` | `(value, *constants)` |
| `sympy.physics.biomechanics.FiberForceLengthActiveDeGroote2016` | `(l_M_tilde, c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11)` | `(value, *constants)` |
| `sympy.physics.biomechanics.FiberForceLengthPassiveDeGroote2016` | `(l_M_tilde, c0, c1)` | `(value, *constants)` |
| `sympy.physics.biomechanics.FiberForceLengthPassiveInverseDeGroote2016` | `(fl_M_pas, c0, c1)` | `(value, *constants)` |
| `sympy.physics.biomechanics.FiberForceVelocityDeGroote2016` | `(v_M_tilde, c0, c1, c2, c3)` | `(value, *constants)` |
| `sympy.physics.biomechanics.FiberForceVelocityInverseDeGroote2016` | `(fv_M, c0, c1, c2, c3)` | `(value, *constants)` |
| `sympy.physics.biomechanics.TendonForceLengthDeGroote2016` | `(l_T_tilde, c0, c1, c2, c3)` | `(value, *constants)` |
| `sympy.physics.biomechanics.TendonForceLengthInverseDeGroote2016` | `(fl_T, c0, c1, c2, c3)` | `(value, *constants)` |
| `sympy.physics.continuum_mechanics.Beam` | `(length, elastic_modulus, second_moment, area=A, variable=x, base_char='C', ild_variable=a)` | `(length, elastic_modulus, second_moment, area=A, variable=x, base_char='C')` |
| `sympy.physics.mechanics.JointsMethod` | `(newtonion, *joints)` | `(newtonian, *joints)` |
| `sympy.physics.mechanics.WrappingGeometryBase` | `()` | `(radius, point)` |
| `sympy.physics.mechanics.WrappingPathway` | `(attachment_1, attachment_2, geometry)` | `(attachment1, attachment2, geometry)` |
| `sympy.physics.mechanics.is_conservative` | `(field)` | `(field, frame)` |
| `sympy.physics.mechanics.is_solenoidal` | `(field)` | `(field, frame)` |
| `sympy.physics.mechanics.joint.Joint` | `(name, parent, child, coordinates=None, speeds=None, parent_point=None, child_point=None, parent_interframe=None, child_interframe=None, parent_axis=None, child_axis=None, parent_joint_pos=None, child_joint_pos=None)` | `()` |
| `sympy.physics.quantum.AntiCommutator` | `(A, B)` | `(a, b)` |
| `sympy.physics.quantum.Commutator` | `(A, B)` | `(a, b)` |
| `sympy.physics.quantum.Dagger` | `(arg)` | `(expr)` |
| `sympy.physics.quantum.DifferentialOperator` | `(*args, **kwargs)` | `(expression, function)` |
| `sympy.physics.quantum.IdentityOperator` | `(*args, **hints)` | `(N=oo)` |
| `sympy.physics.quantum.OuterProduct` | `(*args, **old_assumptions)` | `(ket, bra)` |
| `sympy.physics.quantum.TensorPowerHilbertSpace` | `(*args)` | `(space, power)` |
| `sympy.physics.quantum.Wavefunction` | `(*args)` | `(expression, *variables)` |
| `sympy.physics.quantum.cartesian.PositionBra3D` | `(*args, **kwargs)` | `(x, y, z)` |
| `sympy.physics.quantum.cartesian.PositionKet3D` | `(*args, **kwargs)` | `(x, y, z)` |
| `sympy.physics.quantum.cartesian.PositionState3D` | `(*args, **kwargs)` | `(x, y, z)` |
| `sympy.physics.quantum.cartesian.PxKet` | `(*args, **kwargs)` | `(p)` |
| `sympy.physics.quantum.cartesian.PxOp` | `(*args, **kwargs)` | `(label='Px')` |
| `sympy.physics.quantum.cartesian.XKet` | `(*args, **kwargs)` | `(x)` |
| `sympy.physics.quantum.cartesian.XOp` | `(*args, **kwargs)` | `(label='X')` |
| `sympy.physics.quantum.cartesian.YOp` | `(*args, **kwargs)` | `(label='X')` |
| `sympy.physics.quantum.cartesian.ZOp` | `(*args, **kwargs)` | `(label='X')` |
| `sympy.physics.quantum.cg.CG` | `(j1, m1, j2, m2, j3, m3)` | `(*args)` |
| `sympy.physics.quantum.cg.Wigner3j` | `(j1, m1, j2, m2, j3, m3)` | `(*args)` |
| `sympy.physics.quantum.cg.Wigner6j` | `(j1, j2, j12, j3, j, j23)` | `(*args)` |
| `sympy.physics.quantum.cg.Wigner9j` | `(j1, j2, j12, j3, j4, j34, j13, j24, j)` | `(*args)` |
| `sympy.physics.quantum.circuitplot.Mx` | `(*args, **kwargs)` | `(target)` |
| `sympy.physics.quantum.circuitplot.Mz` | `(*args, **kwargs)` | `(target)` |
| `sympy.physics.quantum.gate.CGate` | `(*args, **kwargs)` | `(controls, gate)` |
| `sympy.physics.quantum.gate.CGateS` | `(*args, **kwargs)` | `(controls, gate)` |
| `sympy.physics.quantum.gate.CNOT` | `(*args, **kwargs)` | `(control, target)` |
| `sympy.physics.quantum.gate.Gate` | `(*args, **kwargs)` | `(*targets)` |
| `sympy.physics.quantum.gate.H` | `(*args, **kwargs)` | `(target)` |
| `sympy.physics.quantum.gate.IdentityGate` | `(*args, **kwargs)` | `(target)` |
| `sympy.physics.quantum.gate.OneQubitGate` | `(*args, **kwargs)` | `(target)` |
| `sympy.physics.quantum.gate.Phase` | `(*args, **kwargs)` | `(target)` |
| `sympy.physics.quantum.gate.SWAP` | `(*args, **kwargs)` | `(first, second)` |
| `sympy.physics.quantum.gate.T` | `(*args, **kwargs)` | `(target)` |
| `sympy.physics.quantum.gate.TwoQubitGate` | `(*args, **kwargs)` | `(first, second)` |
| `sympy.physics.quantum.gate.UGate` | `(*args, **kwargs)` | `(targets, U)` |
| `sympy.physics.quantum.gate.X` | `(*args, **kwargs)` | `(target)` |
| `sympy.physics.quantum.gate.Y` | `(*args, **kwargs)` | `(target)` |
| `sympy.physics.quantum.gate.Z` | `(*args, **kwargs)` | `(target)` |
| `sympy.physics.quantum.grover.OracleGate` | `(*args, **kwargs)` | `(nqubits, f)` |
| `sympy.physics.quantum.grover.WGate` | `(*args, **kwargs)` | `(nqubits)` |
| `sympy.physics.quantum.piab.PIABBra` | `(*args, **kwargs)` | `(n)` |
| `sympy.physics.quantum.piab.PIABKet` | `(*args, **kwargs)` | `(n)` |
| `sympy.physics.quantum.qft.IQFT` | `(*args, **kwargs)` | `(start, finish)` |
| `sympy.physics.quantum.qft.QFT` | `(*args, **kwargs)` | `(start, finish)` |
| `sympy.physics.quantum.qft.Rk` | `(*args)` | `(target, k)` |
| `sympy.physics.quantum.qubit.IntQubit` | `(*args, **kwargs)` | `(n, nqubits=None)` |
| `sympy.physics.quantum.qubit.IntQubitBra` | `(*args, **kwargs)` | `(n, nqubits=None)` |
| `sympy.physics.quantum.qubit.Qubit` | `(*args, **kwargs)` | `(*bits)` |
| `sympy.physics.quantum.qubit.QubitBra` | `(*args, **kwargs)` | `(*bits)` |
| `sympy.physics.quantum.shor.CMod` | `(*args, **kwargs)` | `(t, a, N)` |
| `sympy.physics.quantum.spin.JxBraCoupled` | `(j, m, jn, *jcoupling)` | `(j, m, jn, jcoupling=None)` |
| `sympy.physics.quantum.spin.JxKetCoupled` | `(j, m, jn, *jcoupling)` | `(j, m, jn, jcoupling=None)` |
| `sympy.physics.quantum.spin.JyBraCoupled` | `(j, m, jn, *jcoupling)` | `(j, m, jn, jcoupling=None)` |
| `sympy.physics.quantum.spin.JyKetCoupled` | `(j, m, jn, *jcoupling)` | `(j, m, jn, jcoupling=None)` |
| `sympy.physics.quantum.spin.JzBraCoupled` | `(j, m, jn, *jcoupling)` | `(j, m, jn, jcoupling=None)` |
| `sympy.physics.quantum.spin.JzKetCoupled` | `(j, m, jn, *jcoupling)` | `(j, m, jn, jcoupling=None)` |
| `sympy.physics.quantum.spin.Rotation` | `(*args, **kwargs)` | `(alpha, beta, gamma)` |
| `sympy.physics.quantum.spin.WignerD` | `(*args, **hints)` | `(j, m, mp, alpha, beta, gamma)` |
| `sympy.physics.secondquant.AnnihilateBoson` | `(k)` | `(state)` |
| `sympy.physics.secondquant.AnnihilateFermion` | `(k)` | `(state)` |
| `sympy.physics.secondquant.Bd` | `(k)` | `(state)` |
| `sympy.physics.secondquant.BosonicBasis` | `()` | `(states)` |
| `sympy.physics.secondquant.BosonicOperator` | `(k)` | `(state)` |
| `sympy.physics.secondquant.CreateFermion` | `(k)` | `(state)` |
| `sympy.physics.secondquant.Dagger` | `(arg)` | `(expr)` |
| `sympy.physics.secondquant.NO` | `(arg)` | `(expression)` |
| `sympy.physics.secondquant.substitute_dummies` | `(expr, new_indices=False, pretty_indices={})` | `(expression, new_indices=False, pretty_indices=None)` |
| `sympy.plotting.plot.MatplotlibBackend` | `(*series, **kwargs)` | `(parent)` |
| `sympy.plotting.plot.TextBackend` | `(*args, **kwargs)` | `(parent)` |
| `sympy.polys.numberfields.modules.PowerBasisElement` | `(module, col, denom=1)` | `(module, vector)` |
| `sympy.polys.numberfields.modules.Submodule` | `(parent, matrix, denom=1, mult_tab=None)` | `(parent, matrix)` |
| `sympy.polys.solvers.RawMatrix` | `(*args, **kwargs)` | `(*args)` |
| `sympy.printing.aesaracode.AesaraPrinter` | `(*args, **kwargs)` | `(cache=None, settings=None, **kwargs)` |
| `sympy.printing.rust.TypeCast` | `(expr, type_) -> 'None'` | `(value, type)` |
| `sympy.printing.rust.float_ceiling` | `(arg)` | `(value)` |
| `sympy.printing.rust.float_floor` | `(arg)` | `(value)` |
| `sympy.series.sequences.SeqBase` | `(*args)` | `()` |
| `sympy.series.sequences.SeqExpr` | `(*args)` | `(gen, limits)` |
| `sympy.simplify.fu.sincos_to_sum` | `(expr)` | `(rv, *args, **kw)` |
| `sympy.simplify.hyperexpand.MeijerShiftB` | `(bi)` | `(ai)` |
| `sympy.simplify.hyperexpand.MeijerShiftD` | `(bi)` | `(ai)` |
| `sympy.stats.FiniteDistributionHandmade` | `(*args)` | `(density)` |
| `sympy.stats.crv.ContinuousDomain` | `(symbols, *args)` | `(symbols, set)` |
| `sympy.stats.crv.ContinuousPSpace` | `(*args)` | `(domain, density)` |
| `sympy.stats.crv.SingleContinuousPSpace` | `(s, distribution)` | `(symbol, distribution)` |
| `sympy.stats.frv.ConditionalFiniteDomain` | `(domain, condition)` | `(fulldomain, condition)` |
| `sympy.stats.frv.FiniteDomain` | `(symbols, *args)` | `(elements)` |
| `sympy.stats.frv.SingleFinitePSpace` | `(s, distribution)` | `(symbol, distribution)` |
| `sympy.stats.rv.PSpace` | `(*args)` | `(variables)` |
| `sympy.stats.rv.RandomDomain` | `(symbols, *args)` | `(symbols, set)` |
| `sympy.stats.rv.RandomSymbol` | `(symbol, pspace=None)` | `(name, distribution, group=None)` |
| `sympy.stats.rv.SinglePSpace` | `(s, distribution)` | `(symbol, distribution)` |
| `sympy.stats.rv.sampling_E` | `(expr, given_condition=None, library='scipy', numsamples=1, evalf=True, seed=None, **kwargs)` | `(expr, condition=None, numsamples=1, **kwargs)` |
| `sympy.stats.rv.sampling_P` | `(condition, given_condition=None, library='scipy', numsamples=1, evalf=True, seed=None, **kwargs)` | `(expr, condition=None, numsamples=1, **kwargs)` |
| `sympy.tensor.array.expressions.ArrayAdd` | `(*args, **kwargs)` | `(*args)` |
| `sympy.tensor.array.expressions.ArrayContraction` | `(expr, *contraction_indices, **kwargs)` | `(expr, *groups)` |
| `sympy.tensor.array.expressions.ArrayDiagonal` | `(expr, *diagonal_indices, **kwargs)` | `(expr, *groups)` |
| `sympy.tensor.array.expressions.ArrayElementwiseApplyFunc` | `(function, element)` | `(function, expr)` |
| `sympy.tensor.array.expressions.ArraySymbol` | `(symbol, shape: 'typing.Iterable') -> "'ArraySymbol'"` | `(name, shape)` |
| `sympy.tensor.array.expressions.ArrayTensorProduct` | `(*args, **kwargs)` | `(*args)` |
| `sympy.tensor.array.expressions.PermuteDims` | `(expr, permutation=None, index_order_old=None, index_order_new=None, **kwargs)` | `(expr, permutation=None, index_order_old=None, index_order_new=None)` |
| `sympy.tensor.tensor.TensorType` | `(*args, **kwargs)` | `(index_types, symmetry, **kwargs)` |
| `sympy.testing.pytest.ExceptionInfo` | `(value)` | `(excinfo: 'tuple[type[E], E, TracebackType] \| None', striptext: 'str' = '', traceback: 'Traceback \| None' = None, *, _ispytest: 'bool' = False) -> 'None'` |
| `sympy.utilities.autowrap.CythonCodeWrapper` | `(*args, **kwargs)` | `(generator, filepath=None, flags=None, verbose=False, **kwargs)` |
| `sympy.utilities.autowrap.F2PyCodeWrapper` | `(*args, **kwargs)` | `(generator, filepath=None, flags=None, verbose=False, **kwargs)` |
| `sympy.utilities.autowrap.UfuncifyCodeWrapper` | `(*args, **kwargs)` | `(generator, filepath=None, flags=None, verbose=False, **kwargs)` |
| `sympy.utilities.autowrap.autowrap` | `(expr, language=None, backend='f2py', tempdir=None, args=None, flags=None, verbose=False, helpers=None, code_gen=None, **kwargs)` | `(expr, language=None, backend='f2py', tempdir=None, args=None, flags=None, verbose=False, helpers=None, **kwargs)` |
| `sympy.utilities.iterables.ordered` | `(*args, **kwargs)` | `(seq, keys=None, default=True, warn=False)` |
| `sympy.vector.Laplacian` | `(expr)` | `(*args)` |

## Reproduce and resume

```sh
python tools/sympy_parity.py catalog --out audits/sympy-1.14.0/catalog.json.gz
python tools/sympy_parity.py examples --catalog audits/sympy-1.14.0/catalog.json.gz \
  --out /tmp/symox-parity/nightly.json --timeout 15
python tools/sympy_parity.py report --catalog audits/sympy-1.14.0/catalog.json.gz \
  --examples /tmp/symox-parity/nightly.json --out docs/sympy-parity-audit.md
python tools/sympy_stdlib_certificates.py --out audits/sympy-1.14.0/source-certificates.json
python tools/sympy_parity.py examples --catalog audits/sympy-1.14.0/catalog.json.gz --keywords \
  --out /tmp/symox-parity/keywords.json --timeout 20
```

Use `--module sympy.utilities.iterables` to target a surface, `--limit N` to bound a batch,
`--retry` to rerun existing evidence, and `--pending-only` to continue traversal while retaining
explicitly stale records. Normally a changed native fingerprint causes previous examples to be replayed.

See [the capability contract and remaining algorithmic gaps](algorithmic-parity.md).

## Triage wave 4 (2026-10-22)

This section is maintained by hand and is not produced by `report`; keep it when regenerating.

A full ordinary-example replay against native fingerprint
`dfb31b50e5d39540a91f64fee9bfe81fb2e4aa97e793be3cac10d7b86a7244b6` covered all 2916 documented objects
(30 s per-worker timeout, each replay under `MemoryMax=4G`, `MemorySwapMax=0`, at most two concurrent replays).
Evidence: `audits/sympy-1.14.0/wave4-examples.json.gz`. Result: 1749 matching, 1163 different,
2 incomplete (native timeouts in `Poly.which_all_roots` and `stats.coskewness`), and 2 with no successful
example (`sympy.test` and `sympy.testing.test`, whose examples are all `+SKIP`). Compared with the previous overlay
(1551 matching, 1363 different), 202 objects now match and 4 previously matching objects differ:
`WildFunction` and `core.expr.unchanged` (wrong results), `stats.Zeta` and `PIDController.eval_frequency` (equivalent forms).

Each of the 1167 non-matching objects has exactly one category in
[triage-wave4.json](../audits/sympy-1.14.0/triage-wave4.json), with evidence, the decisive example and,
for g, e, d and h, a minimal reproducer. 512 objects were reviewed by hand (every candidate for g, h, d and e,
every uncertain case and a sample of the rest). The other 655 were classified by string, parse and
numeric-sampling comparison, then spot-checked. When several examples differ, the most severe category wins
(g > e > d > c > f > h > b > i > a).

| Category | Objects |
| --- | ---: |
| a: printing/form only | 312 |
| b: equivalent result in another form | 137 |
| c: missing domain/feature (raises, missing attribute, unevaluated) | 557 |
| d: missing solution / incomplete result | 3 |
| e: insufficient precision or exactness | 4 |
| f: different error type/message | 45 |
| g: wrong result | 71 |
| h: reference wrong, symox right (correction candidate) | 15 |
| i: nondeterministic/environment | 23 |

| Family | a | b | c | d | e | f | g | h | i | Total |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| sympy.stats | 13 | 45 | 33 | 0 | 0 | 0 | 1 | 8 | 0 | 100 |
| sympy.combinatorics | 41 | 1 | 36 | 0 | 0 | 2 | 5 | 0 | 0 | 85 |
| sympy.polys.matrices | 30 | 0 | 46 | 0 | 0 | 2 | 0 | 0 | 7 | 85 |
| sympy.geometry | 2 | 3 | 61 | 0 | 0 | 0 | 3 | 0 | 0 | 69 |
| sympy.codegen | 12 | 4 | 43 | 0 | 2 | 3 | 2 | 0 | 0 | 66 |
| sympy.functions | 8 | 6 | 43 | 0 | 0 | 0 | 4 | 2 | 0 | 63 |
| sympy.physics.continuum_mechanics | 8 | 2 | 35 | 0 | 0 | 0 | 2 | 2 | 0 | 49 |
| sympy.physics.quantum | 13 | 2 | 23 | 0 | 0 | 0 | 11 | 0 | 0 | 49 |
| sympy.core | 8 | 1 | 19 | 0 | 0 | 7 | 12 | 0 | 1 | 48 |
| sympy.matrices | 14 | 3 | 15 | 0 | 0 | 4 | 4 | 1 | 0 | 41 |
| sympy.physics.mechanics | 15 | 2 | 20 | 0 | 0 | 0 | 4 | 0 | 0 | 41 |
| sympy.physics.vector | 26 | 4 | 9 | 0 | 0 | 0 | 2 | 0 | 0 | 41 |
| sympy.physics.secondquant | 8 | 2 | 23 | 0 | 0 | 0 | 3 | 0 | 0 | 36 |
| sympy.physics.control | 7 | 5 | 16 | 0 | 1 | 1 | 4 | 0 | 0 | 34 |
| sympy.vector | 14 | 8 | 6 | 0 | 0 | 1 | 2 | 1 | 0 | 32 |
| sympy.solvers | 7 | 6 | 5 | 1 | 0 | 7 | 2 | 0 | 0 | 28 |
| sympy.algebras | 4 | 2 | 18 | 0 | 0 | 0 | 0 | 0 | 0 | 24 |
| sympy.tensor | 6 | 0 | 16 | 1 | 0 | 0 | 0 | 0 | 0 | 23 |
| sympy.diffgeom | 16 | 5 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 22 |
| sympy.simplify | 6 | 6 | 6 | 0 | 0 | 0 | 0 | 0 | 0 | 18 |
| sympy.physics.optics | 3 | 8 | 3 | 0 | 0 | 0 | 2 | 1 | 0 | 17 |
| sympy.utilities | 5 | 2 | 5 | 0 | 1 | 3 | 0 | 0 | 0 | 16 |
| sympy.printing | 3 | 0 | 10 | 0 | 0 | 1 | 1 | 0 | 0 | 15 |
| sympy.polys.numberfields | 2 | 0 | 12 | 0 | 0 | 0 | 0 | 0 | 0 | 14 |
| sympy.polys.polytools | 2 | 1 | 4 | 0 | 0 | 6 | 0 | 0 | 0 | 13 |
| sympy.holonomic | 7 | 1 | 4 | 0 | 0 | 0 | 1 | 0 | 0 | 13 |
| sympy.integrals | 3 | 4 | 3 | 0 | 0 | 1 | 1 | 0 | 0 | 12 |
| sympy.categories | 7 | 0 | 5 | 0 | 0 | 0 | 0 | 0 | 0 | 12 |
| sympy.sets | 3 | 0 | 7 | 0 | 0 | 0 | 1 | 0 | 0 | 11 |
| sympy.plotting | 2 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 8 | 11 |
| sympy.polys.domains | 1 | 0 | 8 | 0 | 0 | 1 | 0 | 0 | 0 | 10 |
| Other families (26) | 16 | 14 | 21 | 1 | 0 | 6 | 4 | 0 | 7 | 69 |

Wrong results (g), by area. The reproducers are in the JSON records.

- Core and logic: `zoo + zoo` (zoo instead of nan); `comp(pi.n(4), 3.142)`; `all_close` with a relative tolerance;
  `check_assumptions` (True instead of None); `WildFunction.match`; `unchanged(Piecewise, ...)`;
  `iterfreeargs` and `Integral.free_symbols` include bound function variables; `Factors.gcd/lcm/quo/rem`
  with negative exponents; `And.equals` on equivalent formulas.
- Functions and calculus: `genocchi(4, x)` (wrong polynomial); `rf` and `ff` on a `Poly`;
  `Ei(polar_lift(-1)).n()` drops `I*pi`; `Sum(f[n]**2, (n, 0, 3)).doit()`; `limit_seq` of a nested sum
  (0 instead of 4); a Piecewise double integral that leaves the integration variable free;
  `ConditionSet.subs` merging an undecided element.
- Solvers: `solve` with expression or `Indexed` targets returns `[]`; `solve_linear` on an identity.
- Matrices: `Matrix.diag` does not unpack a list; `extract` with boolean masks; `is_symmetric` without
  simplification; `MatrixSymbol.kind` is a string.
- Geometry and vectors: `Point.is_concyclic`, `Ray.contains(Ray)`, the sign of `Plane.angle_between`,
  `BaseVector.equals` across rotated frames, `projection(scalar=True)`.
- Combinatorics and holonomic: `Cycle` iteration, `Permutation.atoms`, `Subset.subset_indices`,
  `Collector.map_relation` and `subword_index`, `HolonomicFunction.evalf` step sizes.
- Physics: quantum (`TensorProduct` of noncommutative symbols, `IntQubit(1, 1)`, Grover oracle and iteration,
  `represent(..., basis=Jx)`, `Wigner3j(6,0,4,0,2,0)`, `enumerate_states`, `Wavefunction` derivative,
  the `TimeDepKet` and `TimeDepBra` labels), secondquant (`NO.doit`, `wicks`, `PermutationOperator.get_permuted`),
  mechanics (`Body` default inertia, `dynamicsymbols` assumptions, `msubs`, `CoordinateSym` equality,
  `WrappingCylinder.geodesic_length`, `WrappingPathway.to_loads`), control (`phase_margin`,
  `pole_zero_numerical_data` order, `StateSpace.dsolve` output, `TransferFunctionMatrix` slicing),
  continuum mechanics (`Beam3D` bending moment and deflection, `Cable.apply_length`), optics (`deviation` sign,
  `fresnel_coefficients` arity), `wigner.real_gaunt` sign, `stats.PoissonProcess` conditional probability.
- Code generation and printing: `codegen.ast.Element` indices, strides and offset; `Node.attr_params`;
  `maple_code(..., assign_to=...)` emits `=` instead of `:=`.

Precision or exactness (e): `expm1(1e-99).evalf()` and `expand_log(log1p(1e-99))` give 0; `lambdify(..., 'mpmath')`
returns a double at higher `mp.dps`; `physics.control.gbt` returns Floats instead of exact coefficients.

Incomplete results (d): `substitution` drops the `2*n*I*pi` ImageSet families; `get_contraction_structure` omits
nested contraction entries; `MusculotendonDeGroote2016.p` lacks four tendon constants.

Correction candidates (h), which need validation before they are recorded as intentional:
`hn1/hn2.rewrite(hankel1/hankel2)` (the order should be `nu + 1/2`); `DiagonalOf.diagonal_length`
(`Min(3, n)`); `ImplicitRegion.regular_point` (SymPy returns a singular point); `gaussian_conj`
(the `(1 - s/f)**2` formula); `Cable.solve` and `Cable.plot_tension` (SymPy reports unbalanced equilibria);
stats CDFs without upper-support branches or exceeding 1 (`BetaNoncentral`, `BoundedPareto`, `PowerFunction`,
`Reciprocal`, `Levy`, `LogitNormal`); the `Dirichlet` marginal; the `GaussianUnitaryEnsemble` normalization.
