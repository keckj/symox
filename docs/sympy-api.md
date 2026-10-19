# SymPy module API inventory

[README](../README.md)

Measured on 2026-10-19 against SymPy 1.14.0.

Source: [official Sphinx module inventory](https://docs.sympy.org/latest/objects.inv) (282 documented modules), extended with parent packages and the README import surfaces.

All 310 resulting module paths are listed below. 0 native module paths are missing; 0 native imports failed; 0 reference imports could not be measured.

Names come from `__all__`, or public callables defined in the module when `__all__` is absent, plus documented functions, classes and data actually bound in the installed reference module. Imported helper names are excluded unless explicitly exported or documented. Constants absent from these sources and individual class methods are outside this module export measure.

The native column measures the corresponding `symox.*` module path. An API available only at the package root is still absent at its compatible module path. Existing README rows with an explicit root mapping retain their original meaning; they can differ from these module rows.

`missing_module` means the module or a parent package is absent. `native_import_error` means an import failed for another reason. `reference_import_error` means the reference requires an unavailable dependency or the upstream documentation names an unimportable module. `no_exports` has no countable names under the stated rule. N/A is never treated as 100%.

Known documentation inventory corrections are recorded in `tools/sympy_public_modules.json`: `Arch` is a class exported by `sympy.physics.continuum_mechanics.arch`, not a module named `Arch`.

These are name-presence measurements. They do not check signatures, domains, algorithms or behavior. Modules overlap; rows must not be added together. This inventory is bounded by the documented modules and listed parents, not a claim to enumerate every internal SymPy module.

Reproduce offline after installing the pinned reference version:

```bash
python tools/sympy_coverage.py --corpus --json /tmp/symox-coverage.json --update-readme
```

The JSON includes missing names and import errors for every module.

<details>
<summary>root</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| [sympy](https://docs.sympy.org/latest/index.html#module-sympy) | 921 / 921 | 0 | 100.00% | available |

</details>

<details>
<summary>abc</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| [sympy.abc](https://docs.sympy.org/latest/modules/abc.html#module-sympy.abc) | 0 / 0 | 0 | N/A | no_exports |

</details>

<details>
<summary>algebras</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| [sympy.algebras](https://docs.sympy.org/latest/modules/algebras.html#module-sympy.algebras) | 1 / 1 | 0 | 100.00% | available |

</details>

<details>
<summary>assumptions</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| [sympy.assumptions](https://docs.sympy.org/latest/modules/assumptions/index.html#module-sympy.assumptions) | 12 / 12 | 0 | 100.00% | available |
| [sympy.assumptions.ask](https://docs.sympy.org/latest/modules/assumptions/ask.html#module-sympy.assumptions.ask) | 6 / 6 | 0 | 100.00% | available |
| [sympy.assumptions.assume](https://docs.sympy.org/latest/modules/assumptions/assume.html#module-sympy.assumptions.assume) | 6 / 6 | 0 | 100.00% | available |
| [sympy.assumptions.predicates](https://docs.sympy.org/latest/modules/assumptions/predicates.html#module-sympy.assumptions.predicates) | 0 / 0 | 0 | N/A | no_exports |
| [sympy.assumptions.refine](https://docs.sympy.org/latest/modules/assumptions/refine.html#module-sympy.assumptions.refine) | 9 / 9 | 0 | 100.00% | available |

</details>

<details>
<summary>calculus</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| [sympy.calculus](https://docs.sympy.org/latest/modules/calculus/index.html#module-sympy.calculus) | 17 / 17 | 0 | 100.00% | available |
| [sympy.calculus.euler](https://docs.sympy.org/latest/modules/calculus/index.html#module-sympy.calculus.euler) | 1 / 1 | 0 | 100.00% | available |
| [sympy.calculus.finite_diff](https://docs.sympy.org/latest/modules/calculus/index.html#module-sympy.calculus.finite_diff) | 3 / 3 | 0 | 100.00% | available |
| [sympy.calculus.singularities](https://docs.sympy.org/latest/modules/calculus/index.html#module-sympy.calculus.singularities) | 7 / 7 | 0 | 100.00% | available |
| [sympy.calculus.util](https://docs.sympy.org/latest/modules/calculus/index.html#module-sympy.calculus.util) | 9 / 9 | 0 | 100.00% | available |

</details>

<details>
<summary>categories</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| [sympy.categories](https://docs.sympy.org/latest/modules/categories.html#module-sympy.categories) | 11 / 11 | 0 | 100.00% | available |
| [sympy.categories.diagram_drawing](https://docs.sympy.org/latest/modules/categories.html#module-sympy.categories.diagram_drawing) | 5 / 5 | 0 | 100.00% | available |

</details>

<details>
<summary>codegen</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| `sympy.codegen` | 13 / 13 | 0 | 100.00% | available |
| [sympy.codegen.algorithms](https://docs.sympy.org/latest/modules/codegen.html#module-sympy.codegen.algorithms) | 2 / 2 | 0 | 100.00% | available |
| [sympy.codegen.approximations](https://docs.sympy.org/latest/modules/codegen.html#module-sympy.codegen.approximations) | 2 / 2 | 0 | 100.00% | available |
| [sympy.codegen.ast](https://docs.sympy.org/latest/modules/codegen.html#module-sympy.codegen.ast) | 43 / 43 | 0 | 100.00% | available |
| [sympy.codegen.cfunctions](https://docs.sympy.org/latest/modules/codegen.html#module-sympy.codegen.cfunctions) | 11 / 11 | 0 | 100.00% | available |
| [sympy.codegen.cnodes](https://docs.sympy.org/latest/modules/codegen.html#module-sympy.codegen.cnodes) | 11 / 11 | 0 | 100.00% | available |
| [sympy.codegen.cutils](https://docs.sympy.org/latest/modules/codegen.html#module-sympy.codegen.cutils) | 1 / 1 | 0 | 100.00% | available |
| [sympy.codegen.cxxnodes](https://docs.sympy.org/latest/modules/codegen.html#module-sympy.codegen.cxxnodes) | 1 / 1 | 0 | 100.00% | available |
| [sympy.codegen.fnodes](https://docs.sympy.org/latest/modules/codegen.html#module-sympy.codegen.fnodes) | 32 / 32 | 0 | 100.00% | available |
| [sympy.codegen.futils](https://docs.sympy.org/latest/modules/codegen.html#module-sympy.codegen.futils) | 1 / 1 | 0 | 100.00% | available |
| [sympy.codegen.matrix_nodes](https://docs.sympy.org/latest/modules/codegen.html#module-sympy.codegen.matrix_nodes) | 1 / 1 | 0 | 100.00% | available |
| [sympy.codegen.pyutils](https://docs.sympy.org/latest/modules/codegen.html#module-sympy.codegen.pyutils) | 1 / 1 | 0 | 100.00% | available |
| [sympy.codegen.rewriting](https://docs.sympy.org/latest/modules/codegen.html#module-sympy.codegen.rewriting) | 18 / 18 | 0 | 100.00% | available |

</details>

<details>
<summary>combinatorics</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| `sympy.combinatorics` | 33 / 33 | 0 | 100.00% | available |
| [sympy.combinatorics.galois](https://docs.sympy.org/latest/modules/combinatorics/galois.html#module-sympy.combinatorics.galois) | 21 / 21 | 0 | 100.00% | available |
| [sympy.combinatorics.generators](https://docs.sympy.org/latest/modules/combinatorics/permutations.html#module-sympy.combinatorics.generators) | 6 / 6 | 0 | 100.00% | available |
| [sympy.combinatorics.graycode](https://docs.sympy.org/latest/modules/combinatorics/graycode.html#module-sympy.combinatorics.graycode) | 6 / 6 | 0 | 100.00% | available |
| [sympy.combinatorics.group_constructs](https://docs.sympy.org/latest/modules/combinatorics/group_constructs.html#module-sympy.combinatorics.group_constructs) | 1 / 1 | 0 | 100.00% | available |
| [sympy.combinatorics.group_numbers](https://docs.sympy.org/latest/modules/combinatorics/group_numbers.html#module-sympy.combinatorics.group_numbers) | 4 / 4 | 0 | 100.00% | available |
| [sympy.combinatorics.named_groups](https://docs.sympy.org/latest/modules/combinatorics/named_groups.html#module-sympy.combinatorics.named_groups) | 8 / 8 | 0 | 100.00% | available |
| [sympy.combinatorics.partitions](https://docs.sympy.org/latest/modules/combinatorics/partitions.html#module-sympy.combinatorics.partitions) | 7 / 7 | 0 | 100.00% | available |
| [sympy.combinatorics.perm_groups](https://docs.sympy.org/latest/modules/combinatorics/perm_groups.html#module-sympy.combinatorics.perm_groups) | 4 / 4 | 0 | 100.00% | available |
| [sympy.combinatorics.permutations](https://docs.sympy.org/latest/modules/combinatorics/permutations.html#module-sympy.combinatorics.permutations) | 4 / 4 | 0 | 100.00% | available |
| [sympy.combinatorics.polyhedron](https://docs.sympy.org/latest/modules/combinatorics/polyhedron.html#module-sympy.combinatorics.polyhedron) | 1 / 1 | 0 | 100.00% | available |
| [sympy.combinatorics.prufer](https://docs.sympy.org/latest/modules/combinatorics/prufer.html#module-sympy.combinatorics.prufer) | 1 / 1 | 0 | 100.00% | available |
| [sympy.combinatorics.subsets](https://docs.sympy.org/latest/modules/combinatorics/subsets.html#module-sympy.combinatorics.subsets) | 2 / 2 | 0 | 100.00% | available |
| [sympy.combinatorics.tensor_can](https://docs.sympy.org/latest/modules/combinatorics/tensor_can.html#module-sympy.combinatorics.tensor_can) | 12 / 12 | 0 | 100.00% | available |
| [sympy.combinatorics.testutil](https://docs.sympy.org/latest/modules/combinatorics/testutil.html#module-sympy.combinatorics.testutil) | 2 / 2 | 0 | 100.00% | available |
| [sympy.combinatorics.util](https://docs.sympy.org/latest/modules/combinatorics/util.html#module-sympy.combinatorics.util) | 0 / 0 | 0 | N/A | no_exports |

</details>

<details>
<summary>core</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| `sympy.core` | 101 / 101 | 0 | 100.00% | available |
| [sympy.core.add](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.add) | 1 / 1 | 0 | 100.00% | available |
| [sympy.core.assumptions](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.assumptions) | 8 / 8 | 0 | 100.00% | available |
| [sympy.core.basic](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.basic) | 3 / 3 | 0 | 100.00% | available |
| [sympy.core.cache](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.cache) | 5 / 5 | 0 | 100.00% | available |
| [sympy.core.containers](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.containers) | 5 / 5 | 0 | 100.00% | available |
| [sympy.core.evalf](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.evalf) | 44 / 44 | 0 | 100.00% | available |
| [sympy.core.expr](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.expr) | 6 / 6 | 0 | 100.00% | available |
| [sympy.core.exprtools](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.exprtools) | 7 / 7 | 0 | 100.00% | available |
| [sympy.core.function](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.function) | 28 / 28 | 0 | 100.00% | available |
| [sympy.core.intfunc](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.intfunc) | 10 / 10 | 0 | 100.00% | available |
| [sympy.core.kind](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.kind) | 3 / 3 | 0 | 100.00% | available |
| [sympy.core.mod](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.mod) | 1 / 1 | 0 | 100.00% | available |
| [sympy.core.mul](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.mul) | 4 / 4 | 0 | 100.00% | available |
| [sympy.core.multidimensional](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.multidimensional) | 4 / 4 | 0 | 100.00% | available |
| [sympy.core.numbers](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.numbers) | 34 / 34 | 0 | 100.00% | available |
| [sympy.core.power](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.power) | 1 / 1 | 0 | 100.00% | available |
| [sympy.core.random](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.random) | 4 / 4 | 0 | 100.00% | available |
| [sympy.core.relational](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.relational) | 14 / 14 | 0 | 100.00% | available |
| [sympy.core.singleton](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.singleton) | 3 / 3 | 0 | 100.00% | available |
| [sympy.core.symbol](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.symbol) | 8 / 8 | 0 | 100.00% | available |
| [sympy.core.sympify](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.sympify) | 4 / 4 | 0 | 100.00% | available |
| [sympy.core.traversal](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.traversal) | 7 / 7 | 0 | 100.00% | available |

</details>

<details>
<summary>crypto</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| `sympy.crypto` | 56 / 56 | 0 | 100.00% | available |
| [sympy.crypto.crypto](https://docs.sympy.org/latest/modules/crypto.html#module-sympy.crypto.crypto) | 57 / 57 | 0 | 100.00% | available |

</details>

<details>
<summary>diffgeom</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| [sympy.diffgeom](https://docs.sympy.org/latest/modules/diffgeom.html#module-sympy.diffgeom) | 24 / 24 | 0 | 100.00% | available |

</details>

<details>
<summary>discrete</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| [sympy.discrete](https://docs.sympy.org/latest/modules/discrete.html#module-sympy.discrete) | 11 / 11 | 0 | 100.00% | available |
| [sympy.discrete.convolutions](https://docs.sympy.org/latest/modules/discrete.html#module-sympy.discrete.convolutions) | 8 / 8 | 0 | 100.00% | available |
| [sympy.discrete.transforms](https://docs.sympy.org/latest/modules/discrete.html#module-sympy.discrete.transforms) | 8 / 8 | 0 | 100.00% | available |

</details>

<details>
<summary>functions</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| [sympy.functions](https://docs.sympy.org/latest/modules/functions/index.html#module-sympy.functions) | 178 / 178 | 0 | 100.00% | available |
| `sympy.functions.elementary` | 0 / 0 | 0 | N/A | no_exports |
| `sympy.functions.elementary.complexes` | 14 / 14 | 0 | 100.00% | available |
| `sympy.functions.special` | 0 / 0 | 0 | N/A | no_exports |
| [sympy.functions.special.bessel](https://docs.sympy.org/latest/modules/functions/special.html#module-sympy.functions.special.bessel) | 21 / 21 | 0 | 100.00% | available |
| [sympy.functions.special.beta_functions](https://docs.sympy.org/latest/modules/functions/special.html#module-sympy.functions.special.beta_functions) | 4 / 4 | 0 | 100.00% | available |
| [sympy.functions.special.elliptic_integrals](https://docs.sympy.org/latest/modules/functions/special.html#module-sympy.functions.special.elliptic_integrals) | 4 / 4 | 0 | 100.00% | available |
| [sympy.functions.special.error_functions](https://docs.sympy.org/latest/modules/functions/special.html#module-sympy.functions.special.error_functions) | 21 / 21 | 0 | 100.00% | available |
| [sympy.functions.special.gamma_functions](https://docs.sympy.org/latest/modules/functions/special.html#module-sympy.functions.special.gamma_functions) | 9 / 9 | 0 | 100.00% | available |
| `sympy.functions.special.hyper` | 17 / 17 | 0 | 100.00% | available |
| [sympy.functions.special.mathieu_functions](https://docs.sympy.org/latest/modules/functions/special.html#module-sympy.functions.special.mathieu_functions) | 5 / 5 | 0 | 100.00% | available |
| [sympy.functions.special.polynomials](https://docs.sympy.org/latest/modules/functions/special.html#module-sympy.functions.special.polynomials) | 14 / 14 | 0 | 100.00% | available |
| [sympy.functions.special.singularity_functions](https://docs.sympy.org/latest/modules/functions/special.html#module-sympy.functions.special.singularity_functions) | 1 / 1 | 0 | 100.00% | available |
| [sympy.functions.special.zeta_functions](https://docs.sympy.org/latest/modules/functions/special.html#module-sympy.functions.special.zeta_functions) | 6 / 6 | 0 | 100.00% | available |

</details>

<details>
<summary>geometry</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| `sympy.geometry` | 30 / 30 | 0 | 100.00% | available |
| [sympy.geometry.curve](https://docs.sympy.org/latest/modules/geometry/curves.html#module-sympy.geometry.curve) | 1 / 1 | 0 | 100.00% | available |
| [sympy.geometry.ellipse](https://docs.sympy.org/latest/modules/geometry/ellipses.html#module-sympy.geometry.ellipse) | 2 / 2 | 0 | 100.00% | available |
| [sympy.geometry.entity](https://docs.sympy.org/latest/modules/geometry/entities.html#module-sympy.geometry.entity) | 5 / 5 | 0 | 100.00% | available |
| [sympy.geometry.line](https://docs.sympy.org/latest/modules/geometry/lines.html#module-sympy.geometry.line) | 12 / 12 | 0 | 100.00% | available |
| [sympy.geometry.plane](https://docs.sympy.org/latest/modules/geometry/plane.html#module-sympy.geometry.plane) | 1 / 1 | 0 | 100.00% | available |
| [sympy.geometry.point](https://docs.sympy.org/latest/modules/geometry/points.html#module-sympy.geometry.point) | 3 / 3 | 0 | 100.00% | available |
| [sympy.geometry.polygon](https://docs.sympy.org/latest/modules/geometry/polygons.html#module-sympy.geometry.polygon) | 5 / 5 | 0 | 100.00% | available |
| [sympy.geometry.util](https://docs.sympy.org/latest/modules/geometry/utils.html#module-sympy.geometry.util) | 9 / 9 | 0 | 100.00% | available |

</details>

<details>
<summary>holonomic</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| [sympy.holonomic](https://docs.sympy.org/latest/modules/holonomic/index.html#module-sympy.holonomic) | 9 / 9 | 0 | 100.00% | available |

</details>

<details>
<summary>integrals</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| [sympy.integrals](https://docs.sympy.org/latest/modules/integrals/integrals.html#module-sympy.integrals) | 30 / 30 | 0 | 100.00% | available |
| [sympy.integrals.intpoly](https://docs.sympy.org/latest/modules/integrals/integrals.html#module-sympy.integrals.intpoly) | 22 / 22 | 0 | 100.00% | available |
| [sympy.integrals.meijerint](https://docs.sympy.org/latest/modules/integrals/g-functions.html#module-sympy.integrals.meijerint) | 3 / 3 | 0 | 100.00% | available |
| [sympy.integrals.meijerint_doc](https://docs.sympy.org/latest/modules/integrals/g-functions.html#module-sympy.integrals.meijerint_doc) | 0 / 0 | 0 | N/A | no_exports |
| [sympy.integrals.transforms](https://docs.sympy.org/latest/modules/integrals/integrals.html#module-sympy.integrals.transforms) | 32 / 32 | 0 | 100.00% | available |

</details>

<details>
<summary>interactive</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| [sympy.interactive](https://docs.sympy.org/latest/modules/interactive.html#module-sympy.interactive) | 3 / 3 | 0 | 100.00% | available |
| [sympy.interactive.printing](https://docs.sympy.org/latest/modules/interactive.html#module-sympy.interactive.printing) | 1 / 1 | 0 | 100.00% | available |
| [sympy.interactive.session](https://docs.sympy.org/latest/modules/interactive.html#module-sympy.interactive.session) | 6 / 6 | 0 | 100.00% | available |

</details>

<details>
<summary>liealgebras</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| [sympy.liealgebras](https://docs.sympy.org/latest/modules/liealgebras/index.html#module-sympy.liealgebras) | 1 / 1 | 0 | 100.00% | available |

</details>

<details>
<summary>logic</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| [sympy.logic](https://docs.sympy.org/latest/modules/logic.html#module-sympy.logic) | 20 / 20 | 0 | 100.00% | available |
| [sympy.logic.inference](https://docs.sympy.org/latest/modules/logic.html#module-sympy.logic.inference) | 7 / 7 | 0 | 100.00% | available |

</details>

<details>
<summary>matrices</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| [sympy.matrices](https://docs.sympy.org/latest/modules/matrices/index.html#module-sympy.matrices) | 73 / 73 | 0 | 100.00% | available |
| [sympy.matrices.expressions](https://docs.sympy.org/latest/modules/matrices/expressions.html#module-sympy.matrices.expressions) | 40 / 40 | 0 | 100.00% | available |
| [sympy.matrices.expressions.blockmatrix](https://docs.sympy.org/latest/modules/matrices/expressions.html#module-sympy.matrices.expressions.blockmatrix) | 16 / 16 | 0 | 100.00% | available |
| [sympy.matrices.immutable](https://docs.sympy.org/latest/modules/matrices/immutablematrices.html#module-sympy.matrices.immutable) | 6 / 6 | 0 | 100.00% | available |
| [sympy.matrices.kind](https://docs.sympy.org/latest/modules/matrices/kind.html#module-sympy.matrices.kind) | 1 / 1 | 0 | 100.00% | available |
| [sympy.matrices.matrixbase](https://docs.sympy.org/latest/modules/matrices/matrices.html#module-sympy.matrices.matrixbase) | 7 / 7 | 0 | 100.00% | available |
| [sympy.matrices.sparse](https://docs.sympy.org/latest/modules/matrices/sparse.html#module-sympy.matrices.sparse) | 3 / 3 | 0 | 100.00% | available |
| [sympy.matrices.sparsetools](https://docs.sympy.org/latest/modules/matrices/sparsetools.html#module-sympy.matrices.sparsetools) | 1 / 1 | 0 | 100.00% | available |

</details>

<details>
<summary>ntheory</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| `sympy.ntheory` | 73 / 73 | 0 | 100.00% | available |
| [sympy.ntheory.bbp_pi](https://docs.sympy.org/latest/modules/ntheory.html#module-sympy.ntheory.bbp_pi) | 1 / 1 | 0 | 100.00% | available |
| [sympy.ntheory.continued_fraction](https://docs.sympy.org/latest/modules/ntheory.html#module-sympy.ntheory.continued_fraction) | 5 / 5 | 0 | 100.00% | available |
| [sympy.ntheory.digits](https://docs.sympy.org/latest/modules/ntheory.html#module-sympy.ntheory.digits) | 3 / 3 | 0 | 100.00% | available |
| [sympy.ntheory.ecm](https://docs.sympy.org/latest/modules/ntheory.html#module-sympy.ntheory.ecm) | 2 / 2 | 0 | 100.00% | available |
| [sympy.ntheory.egyptian_fraction](https://docs.sympy.org/latest/modules/ntheory.html#module-sympy.ntheory.egyptian_fraction) | 6 / 6 | 0 | 100.00% | available |
| [sympy.ntheory.factor_](https://docs.sympy.org/latest/modules/ntheory.html#module-sympy.ntheory.factor_) | 38 / 38 | 0 | 100.00% | available |
| [sympy.ntheory.generate](https://docs.sympy.org/latest/modules/ntheory.html#module-sympy.ntheory.generate) | 11 / 11 | 0 | 100.00% | available |
| [sympy.ntheory.modular](https://docs.sympy.org/latest/modules/ntheory.html#module-sympy.ntheory.modular) | 5 / 5 | 0 | 100.00% | available |
| [sympy.ntheory.multinomial](https://docs.sympy.org/latest/modules/ntheory.html#module-sympy.ntheory.multinomial) | 4 / 4 | 0 | 100.00% | available |
| [sympy.ntheory.partitions_](https://docs.sympy.org/latest/modules/ntheory.html#module-sympy.ntheory.partitions_) | 1 / 1 | 0 | 100.00% | available |
| [sympy.ntheory.primetest](https://docs.sympy.org/latest/modules/ntheory.html#module-sympy.ntheory.primetest) | 12 / 12 | 0 | 100.00% | available |
| [sympy.ntheory.qs](https://docs.sympy.org/latest/modules/ntheory.html#module-sympy.ntheory.qs) | 4 / 4 | 0 | 100.00% | available |
| [sympy.ntheory.residue_ntheory](https://docs.sympy.org/latest/modules/ntheory.html#module-sympy.ntheory.residue_ntheory) | 16 / 16 | 0 | 100.00% | available |

</details>

<details>
<summary>parsing</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| [sympy.parsing](https://docs.sympy.org/latest/modules/parsing.html#module-sympy.parsing) | 1 / 1 | 0 | 100.00% | available |
| [sympy.parsing.sym_expr](https://docs.sympy.org/latest/modules/parsing.html#module-sympy.parsing.sym_expr) | 1 / 1 | 0 | 100.00% | available |
| `sympy.parsing.sympy_parser` | 21 / 21 | 0 | 100.00% | available |

</details>

<details>
<summary>physics</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| [sympy.physics](https://docs.sympy.org/latest/reference/public/physics/index.html#module-sympy.physics) | 5 / 5 | 0 | 100.00% | available |
| [sympy.physics.biomechanics](https://docs.sympy.org/latest/modules/physics/biomechanics/api/index.html#module-sympy.physics.biomechanics) | 15 / 15 | 0 | 100.00% | available |
| [sympy.physics.biomechanics.activation](https://docs.sympy.org/latest/modules/physics/biomechanics/api/activation.html#module-sympy.physics.biomechanics.activation) | 3 / 3 | 0 | 100.00% | available |
| [sympy.physics.biomechanics.curve](https://docs.sympy.org/latest/modules/physics/biomechanics/api/curve.html#module-sympy.physics.biomechanics.curve) | 9 / 9 | 0 | 100.00% | available |
| [sympy.physics.biomechanics.musculotendon](https://docs.sympy.org/latest/modules/physics/biomechanics/api/musculotendon.html#module-sympy.physics.biomechanics.musculotendon) | 3 / 3 | 0 | 100.00% | available |
| `sympy.physics.continuum_mechanics` | 4 / 4 | 0 | 100.00% | available |
| [sympy.physics.continuum_mechanics.arch](https://docs.sympy.org/latest/modules/physics/continuum_mechanics/arches.html) | 1 / 1 | 0 | 100.00% | available |
| [sympy.physics.continuum_mechanics.beam](https://docs.sympy.org/latest/modules/physics/continuum_mechanics/beam.html#module-sympy.physics.continuum_mechanics.beam) | 2 / 2 | 0 | 100.00% | available |
| [sympy.physics.continuum_mechanics.cable](https://docs.sympy.org/latest/modules/physics/continuum_mechanics/cable.html#module-sympy.physics.continuum_mechanics.cable) | 1 / 1 | 0 | 100.00% | available |
| [sympy.physics.continuum_mechanics.truss](https://docs.sympy.org/latest/modules/physics/continuum_mechanics/truss.html#module-sympy.physics.continuum_mechanics.truss) | 1 / 1 | 0 | 100.00% | available |
| [sympy.physics.control](https://docs.sympy.org/latest/modules/physics/control/index.html#module-sympy.physics.control) | 33 / 33 | 0 | 100.00% | available |
| [sympy.physics.control.lti](https://docs.sympy.org/latest/modules/physics/control/lti.html#module-sympy.physics.control.lti) | 16 / 16 | 0 | 100.00% | available |
| `sympy.physics.hep` | 0 / 0 | 0 | N/A | no_exports |
| [sympy.physics.hep.gamma_matrices](https://docs.sympy.org/latest/modules/physics/hep/index.html#module-sympy.physics.hep.gamma_matrices) | 5 / 5 | 0 | 100.00% | available |
| [sympy.physics.hydrogen](https://docs.sympy.org/latest/modules/physics/hydrogen.html#module-sympy.physics.hydrogen) | 4 / 4 | 0 | 100.00% | available |
| [sympy.physics.matrices](https://docs.sympy.org/latest/modules/physics/matrices.html#module-sympy.physics.matrices) | 4 / 4 | 0 | 100.00% | available |
| [sympy.physics.mechanics](https://docs.sympy.org/latest/modules/physics/mechanics/api/index.html#module-sympy.physics.mechanics) | 74 / 74 | 0 | 100.00% | available |
| [sympy.physics.mechanics.actuator](https://docs.sympy.org/latest/modules/physics/mechanics/api/actuator.html#module-sympy.physics.mechanics.actuator) | 7 / 7 | 0 | 100.00% | available |
| [sympy.physics.mechanics.joint](https://docs.sympy.org/latest/modules/physics/mechanics/api/joint.html#module-sympy.physics.mechanics.joint) | 7 / 7 | 0 | 100.00% | available |
| [sympy.physics.mechanics.kane](https://docs.sympy.org/latest/modules/physics/mechanics/api/kane_lagrange.html#module-sympy.physics.mechanics.kane) | 1 / 1 | 0 | 100.00% | available |
| [sympy.physics.mechanics.lagrange](https://docs.sympy.org/latest/modules/physics/mechanics/api/kane_lagrange.html#module-sympy.physics.mechanics.lagrange) | 1 / 1 | 0 | 100.00% | available |
| [sympy.physics.mechanics.linearize](https://docs.sympy.org/latest/modules/physics/mechanics/api/linearize.html#module-sympy.physics.mechanics.linearize) | 1 / 1 | 0 | 100.00% | available |
| [sympy.physics.mechanics.pathway](https://docs.sympy.org/latest/modules/physics/mechanics/api/pathway.html#module-sympy.physics.mechanics.pathway) | 4 / 4 | 0 | 100.00% | available |
| [sympy.physics.mechanics.wrapping_geometry](https://docs.sympy.org/latest/modules/physics/mechanics/api/wrapping_geometry.html#module-sympy.physics.mechanics.wrapping_geometry) | 3 / 3 | 0 | 100.00% | available |
| `sympy.physics.optics` | 39 / 39 | 0 | 100.00% | available |
| [sympy.physics.optics.gaussopt](https://docs.sympy.org/latest/modules/physics/optics/gaussopt.html#module-sympy.physics.optics.gaussopt) | 16 / 16 | 0 | 100.00% | available |
| [sympy.physics.optics.medium](https://docs.sympy.org/latest/modules/physics/optics/medium.html#module-sympy.physics.optics.medium) | 1 / 1 | 0 | 100.00% | available |
| [sympy.physics.optics.polarization](https://docs.sympy.org/latest/modules/physics/optics/polarization.html#module-sympy.physics.optics.polarization) | 11 / 11 | 0 | 100.00% | available |
| [sympy.physics.optics.utils](https://docs.sympy.org/latest/modules/physics/optics/utils.html#module-sympy.physics.optics.utils) | 10 / 10 | 0 | 100.00% | available |
| [sympy.physics.optics.waves](https://docs.sympy.org/latest/modules/physics/optics/waves.html#module-sympy.physics.optics.waves) | 1 / 1 | 0 | 100.00% | available |
| [sympy.physics.paulialgebra](https://docs.sympy.org/latest/modules/physics/paulialgebra.html#module-sympy.physics.paulialgebra) | 1 / 1 | 0 | 100.00% | available |
| [sympy.physics.qho_1d](https://docs.sympy.org/latest/modules/physics/qho_1d.html#module-sympy.physics.qho_1d) | 3 / 3 | 0 | 100.00% | available |
| `sympy.physics.quantum` | 44 / 44 | 0 | 100.00% | available |
| [sympy.physics.quantum.anticommutator](https://docs.sympy.org/latest/modules/physics/quantum/anticommutator.html#module-sympy.physics.quantum.anticommutator) | 1 / 1 | 0 | 100.00% | available |
| [sympy.physics.quantum.cartesian](https://docs.sympy.org/latest/modules/physics/quantum/cartesian.html#module-sympy.physics.quantum.cartesian) | 15 / 15 | 0 | 100.00% | available |
| [sympy.physics.quantum.cg](https://docs.sympy.org/latest/modules/physics/quantum/cg.html#module-sympy.physics.quantum.cg) | 5 / 5 | 0 | 100.00% | available |
| [sympy.physics.quantum.circuitplot](https://docs.sympy.org/latest/modules/physics/quantum/circuitplot.html#module-sympy.physics.quantum.circuitplot) | 7 / 7 | 0 | 100.00% | available |
| [sympy.physics.quantum.commutator](https://docs.sympy.org/latest/modules/physics/quantum/commutator.html#module-sympy.physics.quantum.commutator) | 1 / 1 | 0 | 100.00% | available |
| [sympy.physics.quantum.constants](https://docs.sympy.org/latest/modules/physics/quantum/constants.html#module-sympy.physics.quantum.constants) | 2 / 2 | 0 | 100.00% | available |
| [sympy.physics.quantum.dagger](https://docs.sympy.org/latest/modules/physics/quantum/dagger.html#module-sympy.physics.quantum.dagger) | 1 / 1 | 0 | 100.00% | available |
| [sympy.physics.quantum.gate](https://docs.sympy.org/latest/modules/physics/quantum/gate.html#module-sympy.physics.quantum.gate) | 29 / 29 | 0 | 100.00% | available |
| [sympy.physics.quantum.grover](https://docs.sympy.org/latest/modules/physics/quantum/grover.html#module-sympy.physics.quantum.grover) | 5 / 5 | 0 | 100.00% | available |
| [sympy.physics.quantum.hilbert](https://docs.sympy.org/latest/modules/physics/quantum/hilbert.html#module-sympy.physics.quantum.hilbert) | 8 / 8 | 0 | 100.00% | available |
| [sympy.physics.quantum.innerproduct](https://docs.sympy.org/latest/modules/physics/quantum/innerproduct.html#module-sympy.physics.quantum.innerproduct) | 1 / 1 | 0 | 100.00% | available |
| [sympy.physics.quantum.operator](https://docs.sympy.org/latest/modules/physics/quantum/operator.html#module-sympy.physics.quantum.operator) | 6 / 6 | 0 | 100.00% | available |
| [sympy.physics.quantum.operatorset](https://docs.sympy.org/latest/modules/physics/quantum/operatorset.html#module-sympy.physics.quantum.operatorset) | 2 / 2 | 0 | 100.00% | available |
| [sympy.physics.quantum.piab](https://docs.sympy.org/latest/modules/physics/quantum/piab.html#module-sympy.physics.quantum.piab) | 3 / 3 | 0 | 100.00% | available |
| [sympy.physics.quantum.qapply](https://docs.sympy.org/latest/modules/physics/quantum/qapply.html#module-sympy.physics.quantum.qapply) | 1 / 1 | 0 | 100.00% | available |
| [sympy.physics.quantum.qft](https://docs.sympy.org/latest/modules/physics/quantum/qft.html#module-sympy.physics.quantum.qft) | 4 / 4 | 0 | 100.00% | available |
| [sympy.physics.quantum.qubit](https://docs.sympy.org/latest/modules/physics/quantum/qubit.html#module-sympy.physics.quantum.qubit) | 11 / 11 | 0 | 100.00% | available |
| [sympy.physics.quantum.represent](https://docs.sympy.org/latest/modules/physics/quantum/represent.html#module-sympy.physics.quantum.represent) | 6 / 6 | 0 | 100.00% | available |
| [sympy.physics.quantum.shor](https://docs.sympy.org/latest/modules/physics/quantum/shor.html#module-sympy.physics.quantum.shor) | 6 / 6 | 0 | 100.00% | available |
| [sympy.physics.quantum.spin](https://docs.sympy.org/latest/modules/physics/quantum/spin.html#module-sympy.physics.quantum.spin) | 25 / 25 | 0 | 100.00% | available |
| [sympy.physics.quantum.state](https://docs.sympy.org/latest/modules/physics/quantum/state.html#module-sympy.physics.quantum.state) | 13 / 13 | 0 | 100.00% | available |
| [sympy.physics.quantum.tensorproduct](https://docs.sympy.org/latest/modules/physics/quantum/tensorproduct.html#module-sympy.physics.quantum.tensorproduct) | 2 / 2 | 0 | 100.00% | available |
| [sympy.physics.secondquant](https://docs.sympy.org/latest/modules/physics/secondquant.html#module-sympy.physics.secondquant) | 37 / 37 | 0 | 100.00% | available |
| [sympy.physics.sho](https://docs.sympy.org/latest/modules/physics/sho.html#module-sympy.physics.sho) | 2 / 2 | 0 | 100.00% | available |
| `sympy.physics.units` | 361 / 361 | 0 | 100.00% | available |
| [sympy.physics.units.dimensions](https://docs.sympy.org/latest/modules/physics/units/dimensions.html#module-sympy.physics.units.dimensions) | 2 / 2 | 0 | 100.00% | available |
| [sympy.physics.units.prefixes](https://docs.sympy.org/latest/modules/physics/units/prefixes.html#module-sympy.physics.units.prefixes) | 2 / 2 | 0 | 100.00% | available |
| [sympy.physics.units.quantities](https://docs.sympy.org/latest/modules/physics/units/quantities.html#module-sympy.physics.units.quantities) | 2 / 2 | 0 | 100.00% | available |
| [sympy.physics.units.unitsystem](https://docs.sympy.org/latest/modules/physics/units/unitsystem.html#module-sympy.physics.units.unitsystem) | 1 / 1 | 0 | 100.00% | available |
| [sympy.physics.units.util](https://docs.sympy.org/latest/modules/physics/units/quantities.html#module-sympy.physics.units.util) | 3 / 3 | 0 | 100.00% | available |
| [sympy.physics.vector](https://docs.sympy.org/latest/modules/physics/vector/api/index.html#module-sympy.physics.vector) | 27 / 27 | 0 | 100.00% | available |
| [sympy.physics.vector.functions](https://docs.sympy.org/latest/modules/physics/vector/api/kinematics.html#module-sympy.physics.vector.functions) | 14 / 14 | 0 | 100.00% | available |
| [sympy.physics.vector.point](https://docs.sympy.org/latest/modules/physics/vector/api/kinematics.html#module-sympy.physics.vector.point) | 1 / 1 | 0 | 100.00% | available |
| [sympy.physics.wigner](https://docs.sympy.org/latest/modules/physics/wigner.html#module-sympy.physics.wigner) | 11 / 11 | 0 | 100.00% | available |

</details>

<details>
<summary>plotting</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| `sympy.plotting` | 11 / 11 | 0 | 100.00% | available |
| [sympy.plotting.plot](https://docs.sympy.org/latest/modules/plotting.html#module-sympy.plotting.plot) | 12 / 12 | 0 | 100.00% | available |
| [sympy.plotting.pygletplot](https://docs.sympy.org/latest/modules/plotting.html#module-sympy.plotting.pygletplot) | 1 / 1 | 0 | 100.00% | available |

</details>

<details>
<summary>polys</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| [sympy.polys](https://docs.sympy.org/latest/modules/polys/reference.html#module-sympy.polys) | 176 / 176 | 0 | 100.00% | available |
| `sympy.polys.domains` | 32 / 32 | 0 | 100.00% | available |
| `sympy.polys.fields` | 4 / 4 | 0 | 100.00% | available |
| `sympy.polys.matrices` | 2 / 2 | 0 | 100.00% | available |
| [sympy.polys.matrices._dfm](https://docs.sympy.org/latest/modules/polys/domainmatrix.html#module-sympy.polys.matrices._dfm) | 1 / 1 | 0 | 100.00% | available |
| [sympy.polys.matrices._typing](https://docs.sympy.org/latest/modules/polys/domainmatrix.html#module-sympy.polys.matrices._typing) | 1 / 1 | 0 | 100.00% | available |
| [sympy.polys.matrices.ddm](https://docs.sympy.org/latest/modules/polys/domainmatrix.html#module-sympy.polys.matrices.ddm) | 1 / 1 | 0 | 100.00% | available |
| [sympy.polys.matrices.dense](https://docs.sympy.org/latest/modules/polys/domainmatrix.html#module-sympy.polys.matrices.dense) | 17 / 17 | 0 | 100.00% | available |
| [sympy.polys.matrices.domainmatrix](https://docs.sympy.org/latest/modules/polys/domainmatrix.html#module-sympy.polys.matrices.domainmatrix) | 2 / 2 | 0 | 100.00% | available |
| [sympy.polys.matrices.sdm](https://docs.sympy.org/latest/modules/polys/domainmatrix.html#module-sympy.polys.matrices.sdm) | 13 / 13 | 0 | 100.00% | available |
| `sympy.polys.numberfields` | 10 / 10 | 0 | 100.00% | available |
| [sympy.polys.numberfields.modules](https://docs.sympy.org/latest/modules/polys/numberfields.html#module-sympy.polys.numberfields.modules) | 13 / 13 | 0 | 100.00% | available |
| [sympy.polys.numberfields.subfield](https://docs.sympy.org/latest/modules/polys/numberfields.html#module-sympy.polys.numberfields.subfield) | 3 / 3 | 0 | 100.00% | available |
| [sympy.polys.polyconfig](https://docs.sympy.org/latest/modules/polys/internals.html#module-sympy.polys.polyconfig) | 4 / 4 | 0 | 100.00% | available |
| [sympy.polys.polyoptions](https://docs.sympy.org/latest/modules/polys/internals.html#module-sympy.polys.polyoptions) | 2 / 2 | 0 | 100.00% | available |
| `sympy.polys.rings` | 4 / 4 | 0 | 100.00% | available |
| [sympy.polys.solvers](https://docs.sympy.org/latest/modules/polys/solvers.html#module-sympy.polys.solvers) | 5 / 5 | 0 | 100.00% | available |

</details>

<details>
<summary>printing</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| `sympy.printing` | 41 / 41 | 0 | 100.00% | available |
| [sympy.printing.aesaracode](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.aesaracode) | 4 / 4 | 0 | 100.00% | available |
| [sympy.printing.c](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.c) | 10 / 10 | 0 | 100.00% | available |
| [sympy.printing.codeprinter](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.codeprinter) | 11 / 11 | 0 | 100.00% | available |
| [sympy.printing.conventions](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.conventions) | 2 / 2 | 0 | 100.00% | available |
| [sympy.printing.cxx](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.cxx) | 3 / 3 | 0 | 100.00% | available |
| [sympy.printing.fortran](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.fortran) | 3 / 3 | 0 | 100.00% | available |
| [sympy.printing.gtk](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.gtk) | 1 / 1 | 0 | 100.00% | available |
| [sympy.printing.jscode](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.jscode) | 4 / 4 | 0 | 100.00% | available |
| [sympy.printing.julia](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.julia) | 5 / 5 | 0 | 100.00% | available |
| [sympy.printing.lambdarepr](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.lambdarepr) | 6 / 6 | 0 | 100.00% | available |
| [sympy.printing.latex](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.latex) | 7 / 7 | 0 | 100.00% | available |
| [sympy.printing.maple](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.maple) | 3 / 3 | 0 | 100.00% | available |
| [sympy.printing.mathematica](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.mathematica) | 3 / 3 | 0 | 100.00% | available |
| [sympy.printing.mathml](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.mathml) | 6 / 6 | 0 | 100.00% | available |
| [sympy.printing.octave](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.octave) | 5 / 5 | 0 | 100.00% | available |
| [sympy.printing.precedence](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.precedence) | 12 / 12 | 0 | 100.00% | available |
| [sympy.printing.pretty](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.pretty) | 6 / 6 | 0 | 100.00% | available |
| [sympy.printing.pretty.pretty](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.pretty.pretty) | 19 / 19 | 0 | 100.00% | available |
| [sympy.printing.pretty.pretty_symbology](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.pretty.pretty_symbology) | 25 / 25 | 0 | 100.00% | available |
| [sympy.printing.pretty.stringpict](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.pretty.stringpict) | 2 / 2 | 0 | 100.00% | available |
| [sympy.printing.preview](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.preview) | 3 / 3 | 0 | 100.00% | available |
| [sympy.printing.printer](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.printer) | 3 / 3 | 0 | 100.00% | available |
| [sympy.printing.pycode](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.pycode) | 7 / 7 | 0 | 100.00% | available |
| [sympy.printing.python](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.python) | 3 / 3 | 0 | 100.00% | available |
| [sympy.printing.rcode](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.rcode) | 4 / 4 | 0 | 100.00% | available |
| [sympy.printing.repr](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.repr) | 2 / 2 | 0 | 100.00% | available |
| [sympy.printing.rust](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.rust) | 5 / 5 | 0 | 100.00% | available |
| [sympy.printing.smtlib](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.smtlib) | 2 / 2 | 0 | 100.00% | available |
| [sympy.printing.str](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.str) | 4 / 4 | 0 | 100.00% | available |
| [sympy.printing.tree](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.tree) | 4 / 4 | 0 | 100.00% | available |

</details>

<details>
<summary>series</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| `sympy.series` | 18 / 18 | 0 | 100.00% | available |
| `sympy.series.sequences` | 10 / 10 | 0 | 100.00% | available |

</details>

<details>
<summary>sets</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| `sympy.sets` | 26 / 26 | 0 | 100.00% | available |
| [sympy.sets.conditionset](https://docs.sympy.org/latest/modules/sets.html#module-sympy.sets.conditionset) | 3 / 3 | 0 | 100.00% | available |
| [sympy.sets.fancysets](https://docs.sympy.org/latest/modules/sets.html#module-sympy.sets.fancysets) | 12 / 12 | 0 | 100.00% | available |
| [sympy.sets.powerset](https://docs.sympy.org/latest/modules/sets.html#module-sympy.sets.powerset) | 1 / 1 | 0 | 100.00% | available |
| [sympy.sets.sets](https://docs.sympy.org/latest/modules/sets.html#module-sympy.sets.sets) | 22 / 22 | 0 | 100.00% | available |

</details>

<details>
<summary>simplify</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| `sympy.simplify` | 32 / 32 | 0 | 100.00% | available |
| [sympy.simplify.combsimp](https://docs.sympy.org/latest/modules/simplify/simplify.html#module-sympy.simplify.combsimp) | 1 / 1 | 0 | 100.00% | available |
| [sympy.simplify.cse_main](https://docs.sympy.org/latest/modules/simplify/simplify.html#module-sympy.simplify.cse_main) | 11 / 11 | 0 | 100.00% | available |
| [sympy.simplify.epathtools](https://docs.sympy.org/latest/modules/simplify/simplify.html#module-sympy.simplify.epathtools) | 2 / 2 | 0 | 100.00% | available |
| [sympy.simplify.fu](https://docs.sympy.org/latest/modules/simplify/fu.html#module-sympy.simplify.fu) | 31 / 31 | 0 | 100.00% | available |
| [sympy.simplify.hyperexpand](https://docs.sympy.org/latest/modules/simplify/simplify.html#module-sympy.simplify.hyperexpand) | 37 / 37 | 0 | 100.00% | available |
| [sympy.simplify.hyperexpand_doc](https://docs.sympy.org/latest/modules/simplify/hyperexpand.html#module-sympy.simplify.hyperexpand_doc) | 0 / 0 | 0 | N/A | no_exports |
| [sympy.simplify.powsimp](https://docs.sympy.org/latest/modules/simplify/simplify.html#module-sympy.simplify.powsimp) | 2 / 2 | 0 | 100.00% | available |
| [sympy.simplify.radsimp](https://docs.sympy.org/latest/modules/simplify/simplify.html#module-sympy.simplify.radsimp) | 17 / 17 | 0 | 100.00% | available |
| [sympy.simplify.ratsimp](https://docs.sympy.org/latest/modules/simplify/simplify.html#module-sympy.simplify.ratsimp) | 2 / 2 | 0 | 100.00% | available |
| [sympy.simplify.simplify](https://docs.sympy.org/latest/explanation/gotchas.html#module-sympy.simplify.simplify) | 21 / 21 | 0 | 100.00% | available |
| [sympy.simplify.sqrtdenest](https://docs.sympy.org/latest/modules/simplify/simplify.html#module-sympy.simplify.sqrtdenest) | 6 / 6 | 0 | 100.00% | available |
| [sympy.simplify.trigsimp](https://docs.sympy.org/latest/modules/simplify/simplify.html#module-sympy.simplify.trigsimp) | 5 / 5 | 0 | 100.00% | available |

</details>

<details>
<summary>solvers</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| [sympy.solvers](https://docs.sympy.org/latest/modules/solvers/solvers.html#module-sympy.solvers) | 46 / 46 | 0 | 100.00% | available |
| [sympy.solvers.inequalities](https://docs.sympy.org/latest/modules/solvers/inequalities.html#module-sympy.solvers.inequalities) | 8 / 8 | 0 | 100.00% | available |
| [sympy.solvers.ode](https://docs.sympy.org/latest/modules/solvers/ode.html#module-sympy.solvers.ode) | 11 / 11 | 0 | 100.00% | available |
| [sympy.solvers.ode.ode](https://docs.sympy.org/latest/modules/solvers/ode.html#module-sympy.solvers.ode.ode) | 23 / 23 | 0 | 100.00% | available |
| [sympy.solvers.pde](https://docs.sympy.org/latest/modules/solvers/pde.html#module-sympy.solvers.pde) | 9 / 9 | 0 | 100.00% | available |
| [sympy.solvers.recurr](https://docs.sympy.org/latest/modules/solvers/solvers.html#module-sympy.solvers.recurr) | 4 / 4 | 0 | 100.00% | available |
| [sympy.solvers.simplex](https://docs.sympy.org/latest/modules/solvers/solvers.html#module-sympy.solvers.simplex) | 6 / 6 | 0 | 100.00% | available |
| `sympy.solvers.solvers` | 15 / 15 | 0 | 100.00% | available |
| [sympy.solvers.solveset](https://docs.sympy.org/latest/modules/solvers/solveset.html#module-sympy.solvers.solveset) | 14 / 14 | 0 | 100.00% | available |

</details>

<details>
<summary>stats</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| [sympy.stats](https://docs.sympy.org/latest/modules/stats.html#module-sympy.stats) | 157 / 157 | 0 | 100.00% | available |
| [sympy.stats.crv](https://docs.sympy.org/latest/modules/stats.html#module-sympy.stats.crv) | 9 / 9 | 0 | 100.00% | available |
| [sympy.stats.crv_types](https://docs.sympy.org/latest/modules/stats.html#module-sympy.stats.crv_types) | 54 / 54 | 0 | 100.00% | available |
| [sympy.stats.frv](https://docs.sympy.org/latest/modules/stats.html#module-sympy.stats.frv) | 9 / 9 | 0 | 100.00% | available |
| [sympy.stats.frv_types](https://docs.sympy.org/latest/modules/stats.html#module-sympy.stats.frv_types) | 11 / 11 | 0 | 100.00% | available |
| [sympy.stats.rv](https://docs.sympy.org/latest/modules/stats.html#module-sympy.stats.rv) | 41 / 41 | 0 | 100.00% | available |

</details>

<details>
<summary>tensor</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| [sympy.tensor](https://docs.sympy.org/latest/modules/tensor/index.html#module-sympy.tensor) | 19 / 19 | 0 | 100.00% | available |
| [sympy.tensor.array](https://docs.sympy.org/latest/modules/tensor/array.html#module-sympy.tensor.array) | 16 / 16 | 0 | 100.00% | available |
| [sympy.tensor.array.expressions](https://docs.sympy.org/latest/modules/tensor/array_expressions.html#module-sympy.tensor.array.expressions) | 16 / 16 | 0 | 100.00% | available |
| [sympy.tensor.index_methods](https://docs.sympy.org/latest/modules/tensor/index_methods.html#module-sympy.tensor.index_methods) | 3 / 3 | 0 | 100.00% | available |
| [sympy.tensor.indexed](https://docs.sympy.org/latest/modules/tensor/indexed.html#module-sympy.tensor.indexed) | 4 / 4 | 0 | 100.00% | available |
| [sympy.tensor.tensor](https://docs.sympy.org/latest/modules/tensor/tensor.html#module-sympy.tensor.tensor) | 35 / 35 | 0 | 100.00% | available |
| [sympy.tensor.toperators](https://docs.sympy.org/latest/modules/tensor/toperators.html#module-sympy.tensor.toperators) | 1 / 1 | 0 | 100.00% | available |

</details>

<details>
<summary>testing</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| [sympy.testing](https://docs.sympy.org/latest/modules/testing/index.html#module-sympy.testing) | 2 / 2 | 0 | 100.00% | available |
| [sympy.testing.pytest](https://docs.sympy.org/latest/modules/testing/pytest.html#module-sympy.testing.pytest) | 16 / 16 | 0 | 100.00% | available |
| [sympy.testing.randtest](https://docs.sympy.org/latest/modules/testing/randtest.html#module-sympy.testing.randtest) | 0 / 0 | 0 | N/A | no_exports |
| [sympy.testing.runtests](https://docs.sympy.org/latest/modules/testing/runtests.html#module-sympy.testing.runtests) | 21 / 21 | 0 | 100.00% | available |

</details>

<details>
<summary>utilities</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| [sympy.utilities](https://docs.sympy.org/latest/modules/utilities/index.html#module-sympy.utilities) | 25 / 25 | 0 | 100.00% | available |
| [sympy.utilities.autowrap](https://docs.sympy.org/latest/modules/utilities/autowrap.html#module-sympy.utilities.autowrap) | 9 / 9 | 0 | 100.00% | available |
| [sympy.utilities.codegen](https://docs.sympy.org/latest/modules/utilities/codegen.html#module-sympy.utilities.codegen) | 16 / 16 | 0 | 100.00% | available |
| [sympy.utilities.decorator](https://docs.sympy.org/latest/modules/utilities/decorator.html#module-sympy.utilities.decorator) | 10 / 10 | 0 | 100.00% | available |
| [sympy.utilities.enumerative](https://docs.sympy.org/latest/modules/utilities/enumerative.html#module-sympy.utilities.enumerative) | 6 / 6 | 0 | 100.00% | available |
| [sympy.utilities.exceptions](https://docs.sympy.org/latest/modules/utilities/exceptions.html#module-sympy.utilities.exceptions) | 3 / 3 | 0 | 100.00% | available |
| [sympy.utilities.iterables](https://docs.sympy.org/latest/modules/utilities/iterables.html#module-sympy.utilities.iterables) | 59 / 59 | 0 | 100.00% | available |
| [sympy.utilities.lambdify](https://docs.sympy.org/latest/modules/utilities/lambdify.html#module-sympy.utilities.lambdify) | 3 / 3 | 0 | 100.00% | available |
| [sympy.utilities.memoization](https://docs.sympy.org/latest/modules/utilities/memoization.html#module-sympy.utilities.memoization) | 2 / 2 | 0 | 100.00% | available |
| [sympy.utilities.misc](https://docs.sympy.org/latest/modules/utilities/misc.html#module-sympy.utilities.misc) | 13 / 13 | 0 | 100.00% | available |
| [sympy.utilities.source](https://docs.sympy.org/latest/modules/utilities/source.html#module-sympy.utilities.source) | 2 / 2 | 0 | 100.00% | available |
| [sympy.utilities.timeutils](https://docs.sympy.org/latest/modules/utilities/timeutils.html#module-sympy.utilities.timeutils) | 2 / 2 | 0 | 100.00% | available |

</details>

<details>
<summary>vector</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| [sympy.vector](https://docs.sympy.org/latest/modules/vector/index.html#module-sympy.vector) | 43 / 43 | 0 | 100.00% | available |

</details>
