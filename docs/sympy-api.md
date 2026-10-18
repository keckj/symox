# SymPy module API inventory

[README](../README.md)

Measured on 2026-10-18 against SymPy 1.14.0.

Source: [official Sphinx module inventory](https://docs.sympy.org/latest/objects.inv) (282 documented modules), extended with parent packages and the README import surfaces.

All 311 resulting module paths are listed below. 174 native module paths are missing; 0 native imports failed; 1 reference imports could not be measured.

Names come from `__all__`, or public callables defined in the module when `__all__` is absent, plus documented functions, classes and data actually bound in the installed reference module. Imported helper names are excluded unless explicitly exported or documented. Constants absent from these sources and individual class methods are outside this module export measure.

The native column measures the corresponding `symox.*` module path. An API available only at the package root is still absent at its compatible module path. Existing README rows with an explicit root mapping retain their original meaning; they can differ from these module rows.

`missing_module` means the module or a parent package is absent. `native_import_error` means an import failed for another reason. `reference_import_error` means the reference requires an unavailable dependency or the upstream documentation names an unimportable module. `no_exports` has no countable names under the stated rule. N/A is never treated as 100%.

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
| [sympy.abc](https://docs.sympy.org/latest/modules/abc.html#module-sympy.abc) | 0 / 0 | 0 | N/A | missing_module |

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
| [sympy.assumptions](https://docs.sympy.org/latest/modules/assumptions/index.html#module-sympy.assumptions) | 8 / 12 | 4 | 66.67% | available |
| [sympy.assumptions.ask](https://docs.sympy.org/latest/modules/assumptions/ask.html#module-sympy.assumptions.ask) | 2 / 6 | 4 | 33.33% | available |
| [sympy.assumptions.assume](https://docs.sympy.org/latest/modules/assumptions/assume.html#module-sympy.assumptions.assume) | 4 / 6 | 2 | 66.67% | available |
| [sympy.assumptions.predicates](https://docs.sympy.org/latest/modules/assumptions/predicates.html#module-sympy.assumptions.predicates) | 0 / 0 | 0 | N/A | missing_module |
| [sympy.assumptions.refine](https://docs.sympy.org/latest/modules/assumptions/refine.html#module-sympy.assumptions.refine) | 1 / 9 | 8 | 11.11% | available |

</details>

<details>
<summary>calculus</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| [sympy.calculus](https://docs.sympy.org/latest/modules/calculus/index.html#module-sympy.calculus) | 15 / 17 | 2 | 88.24% | available |
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
| `sympy.codegen` | 2 / 13 | 11 | 15.38% | available |
| [sympy.codegen.algorithms](https://docs.sympy.org/latest/modules/codegen.html#module-sympy.codegen.algorithms) | 0 / 2 | 2 | 0.00% | missing_module |
| [sympy.codegen.approximations](https://docs.sympy.org/latest/modules/codegen.html#module-sympy.codegen.approximations) | 0 / 2 | 2 | 0.00% | missing_module |
| [sympy.codegen.ast](https://docs.sympy.org/latest/modules/codegen.html#module-sympy.codegen.ast) | 3 / 43 | 40 | 6.98% | available |
| [sympy.codegen.cfunctions](https://docs.sympy.org/latest/modules/codegen.html#module-sympy.codegen.cfunctions) | 0 / 11 | 11 | 0.00% | missing_module |
| [sympy.codegen.cnodes](https://docs.sympy.org/latest/modules/codegen.html#module-sympy.codegen.cnodes) | 0 / 11 | 11 | 0.00% | missing_module |
| [sympy.codegen.cutils](https://docs.sympy.org/latest/modules/codegen.html#module-sympy.codegen.cutils) | 0 / 1 | 1 | 0.00% | missing_module |
| [sympy.codegen.cxxnodes](https://docs.sympy.org/latest/modules/codegen.html#module-sympy.codegen.cxxnodes) | 0 / 1 | 1 | 0.00% | missing_module |
| [sympy.codegen.fnodes](https://docs.sympy.org/latest/modules/codegen.html#module-sympy.codegen.fnodes) | 0 / 32 | 32 | 0.00% | missing_module |
| [sympy.codegen.futils](https://docs.sympy.org/latest/modules/codegen.html#module-sympy.codegen.futils) | 0 / 1 | 1 | 0.00% | missing_module |
| [sympy.codegen.matrix_nodes](https://docs.sympy.org/latest/modules/codegen.html#module-sympy.codegen.matrix_nodes) | 0 / 1 | 1 | 0.00% | missing_module |
| [sympy.codegen.pyutils](https://docs.sympy.org/latest/modules/codegen.html#module-sympy.codegen.pyutils) | 0 / 1 | 1 | 0.00% | missing_module |
| [sympy.codegen.rewriting](https://docs.sympy.org/latest/modules/codegen.html#module-sympy.codegen.rewriting) | 0 / 18 | 18 | 0.00% | missing_module |

</details>

<details>
<summary>combinatorics</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| `sympy.combinatorics` | 33 / 33 | 0 | 100.00% | available |
| [sympy.combinatorics.galois](https://docs.sympy.org/latest/modules/combinatorics/galois.html#module-sympy.combinatorics.galois) | 4 / 21 | 17 | 19.05% | available |
| [sympy.combinatorics.generators](https://docs.sympy.org/latest/modules/combinatorics/permutations.html#module-sympy.combinatorics.generators) | 4 / 6 | 2 | 66.67% | available |
| [sympy.combinatorics.graycode](https://docs.sympy.org/latest/modules/combinatorics/graycode.html#module-sympy.combinatorics.graycode) | 1 / 6 | 5 | 16.67% | available |
| [sympy.combinatorics.group_constructs](https://docs.sympy.org/latest/modules/combinatorics/group_constructs.html#module-sympy.combinatorics.group_constructs) | 0 / 1 | 1 | 0.00% | missing_module |
| [sympy.combinatorics.group_numbers](https://docs.sympy.org/latest/modules/combinatorics/group_numbers.html#module-sympy.combinatorics.group_numbers) | 0 / 4 | 4 | 0.00% | missing_module |
| [sympy.combinatorics.named_groups](https://docs.sympy.org/latest/modules/combinatorics/named_groups.html#module-sympy.combinatorics.named_groups) | 6 / 8 | 2 | 75.00% | available |
| [sympy.combinatorics.partitions](https://docs.sympy.org/latest/modules/combinatorics/partitions.html#module-sympy.combinatorics.partitions) | 5 / 7 | 2 | 71.43% | available |
| [sympy.combinatorics.perm_groups](https://docs.sympy.org/latest/modules/combinatorics/perm_groups.html#module-sympy.combinatorics.perm_groups) | 1 / 4 | 3 | 25.00% | available |
| [sympy.combinatorics.permutations](https://docs.sympy.org/latest/modules/combinatorics/permutations.html#module-sympy.combinatorics.permutations) | 2 / 4 | 2 | 50.00% | available |
| [sympy.combinatorics.polyhedron](https://docs.sympy.org/latest/modules/combinatorics/polyhedron.html#module-sympy.combinatorics.polyhedron) | 1 / 1 | 0 | 100.00% | available |
| [sympy.combinatorics.prufer](https://docs.sympy.org/latest/modules/combinatorics/prufer.html#module-sympy.combinatorics.prufer) | 1 / 1 | 0 | 100.00% | available |
| [sympy.combinatorics.subsets](https://docs.sympy.org/latest/modules/combinatorics/subsets.html#module-sympy.combinatorics.subsets) | 1 / 2 | 1 | 50.00% | available |
| [sympy.combinatorics.tensor_can](https://docs.sympy.org/latest/modules/combinatorics/tensor_can.html#module-sympy.combinatorics.tensor_can) | 0 / 12 | 12 | 0.00% | missing_module |
| [sympy.combinatorics.testutil](https://docs.sympy.org/latest/modules/combinatorics/testutil.html#module-sympy.combinatorics.testutil) | 0 / 2 | 2 | 0.00% | missing_module |
| [sympy.combinatorics.util](https://docs.sympy.org/latest/modules/combinatorics/util.html#module-sympy.combinatorics.util) | 0 / 0 | 0 | N/A | missing_module |

</details>

<details>
<summary>core</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| `sympy.core` | 30 / 101 | 71 | 29.70% | available |
| [sympy.core.add](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.add) | 0 / 1 | 1 | 0.00% | missing_module |
| [sympy.core.assumptions](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.assumptions) | 0 / 8 | 8 | 0.00% | missing_module |
| [sympy.core.basic](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.basic) | 0 / 3 | 3 | 0.00% | missing_module |
| [sympy.core.cache](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.cache) | 0 / 5 | 5 | 0.00% | missing_module |
| [sympy.core.containers](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.containers) | 2 / 5 | 3 | 40.00% | available |
| [sympy.core.evalf](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.evalf) | 3 / 44 | 41 | 6.82% | available |
| [sympy.core.expr](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.expr) | 0 / 6 | 6 | 0.00% | missing_module |
| [sympy.core.exprtools](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.exprtools) | 0 / 7 | 7 | 0.00% | missing_module |
| [sympy.core.function](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.function) | 8 / 28 | 20 | 28.57% | available |
| [sympy.core.intfunc](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.intfunc) | 0 / 10 | 10 | 0.00% | missing_module |
| [sympy.core.kind](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.kind) | 0 / 3 | 3 | 0.00% | missing_module |
| [sympy.core.mod](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.mod) | 0 / 1 | 1 | 0.00% | missing_module |
| [sympy.core.mul](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.mul) | 0 / 4 | 4 | 0.00% | missing_module |
| [sympy.core.multidimensional](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.multidimensional) | 1 / 4 | 3 | 25.00% | available |
| [sympy.core.numbers](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.numbers) | 8 / 34 | 26 | 23.53% | available |
| [sympy.core.power](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.power) | 0 / 1 | 1 | 0.00% | missing_module |
| [sympy.core.random](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.random) | 0 / 4 | 4 | 0.00% | missing_module |
| [sympy.core.relational](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.relational) | 0 / 14 | 14 | 0.00% | missing_module |
| [sympy.core.singleton](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.singleton) | 0 / 3 | 3 | 0.00% | missing_module |
| [sympy.core.symbol](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.symbol) | 0 / 8 | 8 | 0.00% | missing_module |
| [sympy.core.sympify](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.sympify) | 2 / 4 | 2 | 50.00% | available |
| [sympy.core.traversal](https://docs.sympy.org/latest/modules/core.html#module-sympy.core.traversal) | 0 / 7 | 7 | 0.00% | missing_module |

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
| [sympy.functions](https://docs.sympy.org/latest/modules/functions/index.html#module-sympy.functions) | 0 / 178 | 178 | 0.00% | available |
| `sympy.functions.elementary` | 0 / 0 | 0 | N/A | no_exports |
| `sympy.functions.elementary.complexes` | 14 / 14 | 0 | 100.00% | available |
| `sympy.functions.special` | 0 / 0 | 0 | N/A | no_exports |
| [sympy.functions.special.bessel](https://docs.sympy.org/latest/modules/functions/special.html#module-sympy.functions.special.bessel) | 21 / 21 | 0 | 100.00% | available |
| [sympy.functions.special.beta_functions](https://docs.sympy.org/latest/modules/functions/special.html#module-sympy.functions.special.beta_functions) | 3 / 4 | 1 | 75.00% | available |
| [sympy.functions.special.elliptic_integrals](https://docs.sympy.org/latest/modules/functions/special.html#module-sympy.functions.special.elliptic_integrals) | 0 / 4 | 4 | 0.00% | missing_module |
| [sympy.functions.special.error_functions](https://docs.sympy.org/latest/modules/functions/special.html#module-sympy.functions.special.error_functions) | 0 / 21 | 21 | 0.00% | missing_module |
| [sympy.functions.special.gamma_functions](https://docs.sympy.org/latest/modules/functions/special.html#module-sympy.functions.special.gamma_functions) | 9 / 9 | 0 | 100.00% | available |
| `sympy.functions.special.hyper` | 17 / 17 | 0 | 100.00% | available |
| [sympy.functions.special.mathieu_functions](https://docs.sympy.org/latest/modules/functions/special.html#module-sympy.functions.special.mathieu_functions) | 5 / 5 | 0 | 100.00% | available |
| [sympy.functions.special.polynomials](https://docs.sympy.org/latest/modules/functions/special.html#module-sympy.functions.special.polynomials) | 0 / 14 | 14 | 0.00% | missing_module |
| [sympy.functions.special.singularity_functions](https://docs.sympy.org/latest/modules/functions/special.html#module-sympy.functions.special.singularity_functions) | 1 / 1 | 0 | 100.00% | available |
| [sympy.functions.special.zeta_functions](https://docs.sympy.org/latest/modules/functions/special.html#module-sympy.functions.special.zeta_functions) | 0 / 6 | 6 | 0.00% | missing_module |

</details>

<details>
<summary>geometry</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| `sympy.geometry` | 30 / 30 | 0 | 100.00% | available |
| [sympy.geometry.curve](https://docs.sympy.org/latest/modules/geometry/curves.html#module-sympy.geometry.curve) | 1 / 1 | 0 | 100.00% | available |
| [sympy.geometry.ellipse](https://docs.sympy.org/latest/modules/geometry/ellipses.html#module-sympy.geometry.ellipse) | 2 / 2 | 0 | 100.00% | available |
| [sympy.geometry.entity](https://docs.sympy.org/latest/modules/geometry/entities.html#module-sympy.geometry.entity) | 1 / 5 | 4 | 20.00% | available |
| [sympy.geometry.line](https://docs.sympy.org/latest/modules/geometry/lines.html#module-sympy.geometry.line) | 10 / 12 | 2 | 83.33% | available |
| [sympy.geometry.plane](https://docs.sympy.org/latest/modules/geometry/plane.html#module-sympy.geometry.plane) | 1 / 1 | 0 | 100.00% | available |
| [sympy.geometry.point](https://docs.sympy.org/latest/modules/geometry/points.html#module-sympy.geometry.point) | 3 / 3 | 0 | 100.00% | available |
| [sympy.geometry.polygon](https://docs.sympy.org/latest/modules/geometry/polygons.html#module-sympy.geometry.polygon) | 3 / 5 | 2 | 60.00% | available |
| [sympy.geometry.util](https://docs.sympy.org/latest/modules/geometry/utils.html#module-sympy.geometry.util) | 7 / 9 | 2 | 77.78% | available |

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
| [sympy.integrals](https://docs.sympy.org/latest/modules/integrals/integrals.html#module-sympy.integrals) | 0 / 30 | 30 | 0.00% | missing_module |
| [sympy.integrals.intpoly](https://docs.sympy.org/latest/modules/integrals/integrals.html#module-sympy.integrals.intpoly) | 0 / 22 | 22 | 0.00% | missing_module |
| [sympy.integrals.meijerint](https://docs.sympy.org/latest/modules/integrals/g-functions.html#module-sympy.integrals.meijerint) | 0 / 3 | 3 | 0.00% | missing_module |
| [sympy.integrals.meijerint_doc](https://docs.sympy.org/latest/modules/integrals/g-functions.html#module-sympy.integrals.meijerint_doc) | 0 / 0 | 0 | N/A | missing_module |
| [sympy.integrals.transforms](https://docs.sympy.org/latest/modules/integrals/integrals.html#module-sympy.integrals.transforms) | 0 / 32 | 32 | 0.00% | missing_module |

</details>

<details>
<summary>interactive</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| [sympy.interactive](https://docs.sympy.org/latest/modules/interactive.html#module-sympy.interactive) | 2 / 3 | 1 | 66.67% | available |
| [sympy.interactive.printing](https://docs.sympy.org/latest/modules/interactive.html#module-sympy.interactive.printing) | 0 / 1 | 1 | 0.00% | missing_module |
| [sympy.interactive.session](https://docs.sympy.org/latest/modules/interactive.html#module-sympy.interactive.session) | 0 / 6 | 6 | 0.00% | missing_module |

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
| [sympy.logic](https://docs.sympy.org/latest/modules/logic.html#module-sympy.logic) | 0 / 20 | 20 | 0.00% | missing_module |
| [sympy.logic.inference](https://docs.sympy.org/latest/modules/logic.html#module-sympy.logic.inference) | 0 / 7 | 7 | 0.00% | missing_module |

</details>

<details>
<summary>matrices</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| [sympy.matrices](https://docs.sympy.org/latest/modules/matrices/index.html#module-sympy.matrices) | 42 / 73 | 31 | 57.53% | available |
| [sympy.matrices.expressions](https://docs.sympy.org/latest/modules/matrices/expressions.html#module-sympy.matrices.expressions) | 30 / 40 | 10 | 75.00% | available |
| [sympy.matrices.expressions.blockmatrix](https://docs.sympy.org/latest/modules/matrices/expressions.html#module-sympy.matrices.expressions.blockmatrix) | 3 / 16 | 13 | 18.75% | available |
| [sympy.matrices.immutable](https://docs.sympy.org/latest/modules/matrices/immutablematrices.html#module-sympy.matrices.immutable) | 0 / 6 | 6 | 0.00% | missing_module |
| [sympy.matrices.kind](https://docs.sympy.org/latest/modules/matrices/kind.html#module-sympy.matrices.kind) | 0 / 1 | 1 | 0.00% | missing_module |
| [sympy.matrices.matrixbase](https://docs.sympy.org/latest/modules/matrices/matrices.html#module-sympy.matrices.matrixbase) | 0 / 7 | 7 | 0.00% | missing_module |
| [sympy.matrices.sparse](https://docs.sympy.org/latest/modules/matrices/sparse.html#module-sympy.matrices.sparse) | 0 / 3 | 3 | 0.00% | missing_module |
| [sympy.matrices.sparsetools](https://docs.sympy.org/latest/modules/matrices/sparsetools.html#module-sympy.matrices.sparsetools) | 0 / 1 | 1 | 0.00% | missing_module |

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
| [sympy.ntheory.factor_](https://docs.sympy.org/latest/modules/ntheory.html#module-sympy.ntheory.factor_) | 37 / 38 | 1 | 97.37% | available |
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
| [sympy.parsing.sym_expr](https://docs.sympy.org/latest/modules/parsing.html#module-sympy.parsing.sym_expr) | 0 / 1 | 1 | 0.00% | missing_module |
| `sympy.parsing.sympy_parser` | 21 / 21 | 0 | 100.00% | available |

</details>

<details>
<summary>physics</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| [sympy.physics](https://docs.sympy.org/latest/reference/public/physics/index.html#module-sympy.physics) | 0 / 5 | 5 | 0.00% | missing_module |
| [sympy.physics.biomechanics](https://docs.sympy.org/latest/modules/physics/biomechanics/api/index.html#module-sympy.physics.biomechanics) | 0 / 15 | 15 | 0.00% | missing_module |
| [sympy.physics.biomechanics.activation](https://docs.sympy.org/latest/modules/physics/biomechanics/api/activation.html#module-sympy.physics.biomechanics.activation) | 0 / 3 | 3 | 0.00% | missing_module |
| [sympy.physics.biomechanics.curve](https://docs.sympy.org/latest/modules/physics/biomechanics/api/curve.html#module-sympy.physics.biomechanics.curve) | 0 / 9 | 9 | 0.00% | missing_module |
| [sympy.physics.biomechanics.musculotendon](https://docs.sympy.org/latest/modules/physics/biomechanics/api/musculotendon.html#module-sympy.physics.biomechanics.musculotendon) | 0 / 3 | 3 | 0.00% | missing_module |
| `sympy.physics.continuum_mechanics` | 0 / 4 | 4 | 0.00% | missing_module |
| [sympy.physics.continuum_mechanics.Arch](https://docs.sympy.org/latest/modules/physics/continuum_mechanics/arches.html#module-sympy.physics.continuum_mechanics.Arch) | N/A | N/A | N/A | reference_import_error |
| `sympy.physics.continuum_mechanics.arch` | 0 / 1 | 1 | 0.00% | missing_module |
| [sympy.physics.continuum_mechanics.beam](https://docs.sympy.org/latest/modules/physics/continuum_mechanics/beam.html#module-sympy.physics.continuum_mechanics.beam) | 0 / 2 | 2 | 0.00% | missing_module |
| [sympy.physics.continuum_mechanics.cable](https://docs.sympy.org/latest/modules/physics/continuum_mechanics/cable.html#module-sympy.physics.continuum_mechanics.cable) | 0 / 1 | 1 | 0.00% | missing_module |
| [sympy.physics.continuum_mechanics.truss](https://docs.sympy.org/latest/modules/physics/continuum_mechanics/truss.html#module-sympy.physics.continuum_mechanics.truss) | 0 / 1 | 1 | 0.00% | missing_module |
| [sympy.physics.control](https://docs.sympy.org/latest/modules/physics/control/index.html#module-sympy.physics.control) | 0 / 33 | 33 | 0.00% | missing_module |
| [sympy.physics.control.lti](https://docs.sympy.org/latest/modules/physics/control/lti.html#module-sympy.physics.control.lti) | 0 / 16 | 16 | 0.00% | missing_module |
| `sympy.physics.hep` | 0 / 0 | 0 | N/A | missing_module |
| [sympy.physics.hep.gamma_matrices](https://docs.sympy.org/latest/modules/physics/hep/index.html#module-sympy.physics.hep.gamma_matrices) | 0 / 5 | 5 | 0.00% | missing_module |
| [sympy.physics.hydrogen](https://docs.sympy.org/latest/modules/physics/hydrogen.html#module-sympy.physics.hydrogen) | 0 / 4 | 4 | 0.00% | missing_module |
| [sympy.physics.matrices](https://docs.sympy.org/latest/modules/physics/matrices.html#module-sympy.physics.matrices) | 0 / 4 | 4 | 0.00% | missing_module |
| [sympy.physics.mechanics](https://docs.sympy.org/latest/modules/physics/mechanics/api/index.html#module-sympy.physics.mechanics) | 0 / 74 | 74 | 0.00% | missing_module |
| [sympy.physics.mechanics.actuator](https://docs.sympy.org/latest/modules/physics/mechanics/api/actuator.html#module-sympy.physics.mechanics.actuator) | 0 / 7 | 7 | 0.00% | missing_module |
| [sympy.physics.mechanics.joint](https://docs.sympy.org/latest/modules/physics/mechanics/api/joint.html#module-sympy.physics.mechanics.joint) | 0 / 7 | 7 | 0.00% | missing_module |
| [sympy.physics.mechanics.kane](https://docs.sympy.org/latest/modules/physics/mechanics/api/kane_lagrange.html#module-sympy.physics.mechanics.kane) | 0 / 1 | 1 | 0.00% | missing_module |
| [sympy.physics.mechanics.lagrange](https://docs.sympy.org/latest/modules/physics/mechanics/api/kane_lagrange.html#module-sympy.physics.mechanics.lagrange) | 0 / 1 | 1 | 0.00% | missing_module |
| [sympy.physics.mechanics.linearize](https://docs.sympy.org/latest/modules/physics/mechanics/api/linearize.html#module-sympy.physics.mechanics.linearize) | 0 / 1 | 1 | 0.00% | missing_module |
| [sympy.physics.mechanics.pathway](https://docs.sympy.org/latest/modules/physics/mechanics/api/pathway.html#module-sympy.physics.mechanics.pathway) | 0 / 4 | 4 | 0.00% | missing_module |
| [sympy.physics.mechanics.wrapping_geometry](https://docs.sympy.org/latest/modules/physics/mechanics/api/wrapping_geometry.html#module-sympy.physics.mechanics.wrapping_geometry) | 0 / 3 | 3 | 0.00% | missing_module |
| `sympy.physics.optics` | 0 / 39 | 39 | 0.00% | missing_module |
| [sympy.physics.optics.gaussopt](https://docs.sympy.org/latest/modules/physics/optics/gaussopt.html#module-sympy.physics.optics.gaussopt) | 0 / 16 | 16 | 0.00% | missing_module |
| [sympy.physics.optics.medium](https://docs.sympy.org/latest/modules/physics/optics/medium.html#module-sympy.physics.optics.medium) | 0 / 1 | 1 | 0.00% | missing_module |
| [sympy.physics.optics.polarization](https://docs.sympy.org/latest/modules/physics/optics/polarization.html#module-sympy.physics.optics.polarization) | 0 / 11 | 11 | 0.00% | missing_module |
| [sympy.physics.optics.utils](https://docs.sympy.org/latest/modules/physics/optics/utils.html#module-sympy.physics.optics.utils) | 0 / 10 | 10 | 0.00% | missing_module |
| [sympy.physics.optics.waves](https://docs.sympy.org/latest/modules/physics/optics/waves.html#module-sympy.physics.optics.waves) | 0 / 1 | 1 | 0.00% | missing_module |
| [sympy.physics.paulialgebra](https://docs.sympy.org/latest/modules/physics/paulialgebra.html#module-sympy.physics.paulialgebra) | 0 / 1 | 1 | 0.00% | missing_module |
| [sympy.physics.qho_1d](https://docs.sympy.org/latest/modules/physics/qho_1d.html#module-sympy.physics.qho_1d) | 0 / 3 | 3 | 0.00% | missing_module |
| `sympy.physics.quantum` | 0 / 44 | 44 | 0.00% | missing_module |
| [sympy.physics.quantum.anticommutator](https://docs.sympy.org/latest/modules/physics/quantum/anticommutator.html#module-sympy.physics.quantum.anticommutator) | 0 / 1 | 1 | 0.00% | missing_module |
| [sympy.physics.quantum.cartesian](https://docs.sympy.org/latest/modules/physics/quantum/cartesian.html#module-sympy.physics.quantum.cartesian) | 0 / 15 | 15 | 0.00% | missing_module |
| [sympy.physics.quantum.cg](https://docs.sympy.org/latest/modules/physics/quantum/cg.html#module-sympy.physics.quantum.cg) | 0 / 5 | 5 | 0.00% | missing_module |
| [sympy.physics.quantum.circuitplot](https://docs.sympy.org/latest/modules/physics/quantum/circuitplot.html#module-sympy.physics.quantum.circuitplot) | 0 / 7 | 7 | 0.00% | missing_module |
| [sympy.physics.quantum.commutator](https://docs.sympy.org/latest/modules/physics/quantum/commutator.html#module-sympy.physics.quantum.commutator) | 0 / 1 | 1 | 0.00% | missing_module |
| [sympy.physics.quantum.constants](https://docs.sympy.org/latest/modules/physics/quantum/constants.html#module-sympy.physics.quantum.constants) | 0 / 2 | 2 | 0.00% | missing_module |
| [sympy.physics.quantum.dagger](https://docs.sympy.org/latest/modules/physics/quantum/dagger.html#module-sympy.physics.quantum.dagger) | 0 / 1 | 1 | 0.00% | missing_module |
| [sympy.physics.quantum.gate](https://docs.sympy.org/latest/modules/physics/quantum/gate.html#module-sympy.physics.quantum.gate) | 0 / 29 | 29 | 0.00% | missing_module |
| [sympy.physics.quantum.grover](https://docs.sympy.org/latest/modules/physics/quantum/grover.html#module-sympy.physics.quantum.grover) | 0 / 5 | 5 | 0.00% | missing_module |
| [sympy.physics.quantum.hilbert](https://docs.sympy.org/latest/modules/physics/quantum/hilbert.html#module-sympy.physics.quantum.hilbert) | 0 / 8 | 8 | 0.00% | missing_module |
| [sympy.physics.quantum.innerproduct](https://docs.sympy.org/latest/modules/physics/quantum/innerproduct.html#module-sympy.physics.quantum.innerproduct) | 0 / 1 | 1 | 0.00% | missing_module |
| [sympy.physics.quantum.operator](https://docs.sympy.org/latest/modules/physics/quantum/operator.html#module-sympy.physics.quantum.operator) | 0 / 6 | 6 | 0.00% | missing_module |
| [sympy.physics.quantum.operatorset](https://docs.sympy.org/latest/modules/physics/quantum/operatorset.html#module-sympy.physics.quantum.operatorset) | 0 / 2 | 2 | 0.00% | missing_module |
| [sympy.physics.quantum.piab](https://docs.sympy.org/latest/modules/physics/quantum/piab.html#module-sympy.physics.quantum.piab) | 0 / 3 | 3 | 0.00% | missing_module |
| [sympy.physics.quantum.qapply](https://docs.sympy.org/latest/modules/physics/quantum/qapply.html#module-sympy.physics.quantum.qapply) | 0 / 1 | 1 | 0.00% | missing_module |
| [sympy.physics.quantum.qft](https://docs.sympy.org/latest/modules/physics/quantum/qft.html#module-sympy.physics.quantum.qft) | 0 / 4 | 4 | 0.00% | missing_module |
| [sympy.physics.quantum.qubit](https://docs.sympy.org/latest/modules/physics/quantum/qubit.html#module-sympy.physics.quantum.qubit) | 0 / 11 | 11 | 0.00% | missing_module |
| [sympy.physics.quantum.represent](https://docs.sympy.org/latest/modules/physics/quantum/represent.html#module-sympy.physics.quantum.represent) | 0 / 6 | 6 | 0.00% | missing_module |
| [sympy.physics.quantum.shor](https://docs.sympy.org/latest/modules/physics/quantum/shor.html#module-sympy.physics.quantum.shor) | 0 / 6 | 6 | 0.00% | missing_module |
| [sympy.physics.quantum.spin](https://docs.sympy.org/latest/modules/physics/quantum/spin.html#module-sympy.physics.quantum.spin) | 0 / 25 | 25 | 0.00% | missing_module |
| [sympy.physics.quantum.state](https://docs.sympy.org/latest/modules/physics/quantum/state.html#module-sympy.physics.quantum.state) | 0 / 13 | 13 | 0.00% | missing_module |
| [sympy.physics.quantum.tensorproduct](https://docs.sympy.org/latest/modules/physics/quantum/tensorproduct.html#module-sympy.physics.quantum.tensorproduct) | 0 / 2 | 2 | 0.00% | missing_module |
| [sympy.physics.secondquant](https://docs.sympy.org/latest/modules/physics/secondquant.html#module-sympy.physics.secondquant) | 0 / 37 | 37 | 0.00% | missing_module |
| [sympy.physics.sho](https://docs.sympy.org/latest/modules/physics/sho.html#module-sympy.physics.sho) | 0 / 2 | 2 | 0.00% | missing_module |
| `sympy.physics.units` | 0 / 361 | 361 | 0.00% | missing_module |
| [sympy.physics.units.dimensions](https://docs.sympy.org/latest/modules/physics/units/dimensions.html#module-sympy.physics.units.dimensions) | 0 / 2 | 2 | 0.00% | missing_module |
| [sympy.physics.units.prefixes](https://docs.sympy.org/latest/modules/physics/units/prefixes.html#module-sympy.physics.units.prefixes) | 0 / 2 | 2 | 0.00% | missing_module |
| [sympy.physics.units.quantities](https://docs.sympy.org/latest/modules/physics/units/quantities.html#module-sympy.physics.units.quantities) | 0 / 2 | 2 | 0.00% | missing_module |
| [sympy.physics.units.unitsystem](https://docs.sympy.org/latest/modules/physics/units/unitsystem.html#module-sympy.physics.units.unitsystem) | 0 / 1 | 1 | 0.00% | missing_module |
| [sympy.physics.units.util](https://docs.sympy.org/latest/modules/physics/units/quantities.html#module-sympy.physics.units.util) | 0 / 3 | 3 | 0.00% | missing_module |
| [sympy.physics.vector](https://docs.sympy.org/latest/modules/physics/vector/api/index.html#module-sympy.physics.vector) | 0 / 27 | 27 | 0.00% | missing_module |
| [sympy.physics.vector.functions](https://docs.sympy.org/latest/modules/physics/vector/api/kinematics.html#module-sympy.physics.vector.functions) | 0 / 14 | 14 | 0.00% | missing_module |
| [sympy.physics.vector.point](https://docs.sympy.org/latest/modules/physics/vector/api/kinematics.html#module-sympy.physics.vector.point) | 0 / 1 | 1 | 0.00% | missing_module |
| [sympy.physics.wigner](https://docs.sympy.org/latest/modules/physics/wigner.html#module-sympy.physics.wigner) | 0 / 11 | 11 | 0.00% | missing_module |

</details>

<details>
<summary>plotting</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| `sympy.plotting` | 11 / 11 | 0 | 100.00% | available |
| [sympy.plotting.plot](https://docs.sympy.org/latest/modules/plotting.html#module-sympy.plotting.plot) | 10 / 12 | 2 | 83.33% | available |
| [sympy.plotting.pygletplot](https://docs.sympy.org/latest/modules/plotting.html#module-sympy.plotting.pygletplot) | 0 / 1 | 1 | 0.00% | missing_module |

</details>

<details>
<summary>polys</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| [sympy.polys](https://docs.sympy.org/latest/modules/polys/reference.html#module-sympy.polys) | 90 / 176 | 86 | 51.14% | available |
| `sympy.polys.domains` | 32 / 32 | 0 | 100.00% | available |
| `sympy.polys.fields` | 4 / 4 | 0 | 100.00% | available |
| `sympy.polys.matrices` | 2 / 2 | 0 | 100.00% | available |
| [sympy.polys.matrices._dfm](https://docs.sympy.org/latest/modules/polys/domainmatrix.html#module-sympy.polys.matrices._dfm) | 0 / 1 | 1 | 0.00% | missing_module |
| [sympy.polys.matrices._typing](https://docs.sympy.org/latest/modules/polys/domainmatrix.html#module-sympy.polys.matrices._typing) | 0 / 1 | 1 | 0.00% | missing_module |
| [sympy.polys.matrices.ddm](https://docs.sympy.org/latest/modules/polys/domainmatrix.html#module-sympy.polys.matrices.ddm) | 0 / 1 | 1 | 0.00% | missing_module |
| [sympy.polys.matrices.dense](https://docs.sympy.org/latest/modules/polys/domainmatrix.html#module-sympy.polys.matrices.dense) | 0 / 17 | 17 | 0.00% | missing_module |
| [sympy.polys.matrices.domainmatrix](https://docs.sympy.org/latest/modules/polys/domainmatrix.html#module-sympy.polys.matrices.domainmatrix) | 1 / 2 | 1 | 50.00% | available |
| [sympy.polys.matrices.sdm](https://docs.sympy.org/latest/modules/polys/domainmatrix.html#module-sympy.polys.matrices.sdm) | 0 / 13 | 13 | 0.00% | missing_module |
| `sympy.polys.numberfields` | 10 / 10 | 0 | 100.00% | available |
| [sympy.polys.numberfields.modules](https://docs.sympy.org/latest/modules/polys/numberfields.html#module-sympy.polys.numberfields.modules) | 3 / 13 | 10 | 23.08% | available |
| [sympy.polys.numberfields.subfield](https://docs.sympy.org/latest/modules/polys/numberfields.html#module-sympy.polys.numberfields.subfield) | 3 / 3 | 0 | 100.00% | available |
| [sympy.polys.polyconfig](https://docs.sympy.org/latest/modules/polys/internals.html#module-sympy.polys.polyconfig) | 0 / 4 | 4 | 0.00% | missing_module |
| [sympy.polys.polyoptions](https://docs.sympy.org/latest/modules/polys/internals.html#module-sympy.polys.polyoptions) | 1 / 2 | 1 | 50.00% | available |
| `sympy.polys.rings` | 4 / 4 | 0 | 100.00% | available |
| [sympy.polys.solvers](https://docs.sympy.org/latest/modules/polys/solvers.html#module-sympy.polys.solvers) | 0 / 5 | 5 | 0.00% | missing_module |

</details>

<details>
<summary>printing</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| `sympy.printing` | 41 / 41 | 0 | 100.00% | available |
| [sympy.printing.aesaracode](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.aesaracode) | 0 / 4 | 4 | 0.00% | missing_module |
| [sympy.printing.c](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.c) | 1 / 10 | 9 | 10.00% | available |
| [sympy.printing.codeprinter](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.codeprinter) | 0 / 11 | 11 | 0.00% | missing_module |
| [sympy.printing.conventions](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.conventions) | 0 / 2 | 2 | 0.00% | missing_module |
| [sympy.printing.cxx](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.cxx) | 3 / 3 | 0 | 100.00% | available |
| [sympy.printing.fortran](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.fortran) | 1 / 3 | 2 | 33.33% | available |
| [sympy.printing.gtk](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.gtk) | 1 / 1 | 0 | 100.00% | available |
| [sympy.printing.jscode](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.jscode) | 2 / 4 | 2 | 50.00% | available |
| [sympy.printing.julia](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.julia) | 2 / 5 | 3 | 40.00% | available |
| [sympy.printing.lambdarepr](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.lambdarepr) | 0 / 6 | 6 | 0.00% | missing_module |
| [sympy.printing.latex](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.latex) | 4 / 7 | 3 | 57.14% | available |
| [sympy.printing.maple](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.maple) | 2 / 3 | 1 | 66.67% | available |
| [sympy.printing.mathematica](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.mathematica) | 2 / 3 | 1 | 66.67% | available |
| [sympy.printing.mathml](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.mathml) | 5 / 6 | 1 | 83.33% | available |
| [sympy.printing.octave](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.octave) | 2 / 5 | 3 | 40.00% | available |
| [sympy.printing.precedence](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.precedence) | 2 / 12 | 10 | 16.67% | available |
| [sympy.printing.pretty](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.pretty) | 3 / 6 | 3 | 50.00% | available |
| [sympy.printing.pretty.pretty](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.pretty.pretty) | 6 / 19 | 13 | 31.58% | available |
| [sympy.printing.pretty.pretty_symbology](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.pretty.pretty_symbology) | 0 / 25 | 25 | 0.00% | missing_module |
| [sympy.printing.pretty.stringpict](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.pretty.stringpict) | 0 / 2 | 2 | 0.00% | missing_module |
| [sympy.printing.preview](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.preview) | 3 / 3 | 0 | 100.00% | available |
| [sympy.printing.printer](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.printer) | 1 / 3 | 2 | 33.33% | available |
| [sympy.printing.pycode](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.pycode) | 1 / 7 | 6 | 14.29% | available |
| [sympy.printing.python](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.python) | 2 / 3 | 1 | 66.67% | available |
| [sympy.printing.rcode](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.rcode) | 2 / 4 | 2 | 50.00% | available |
| [sympy.printing.repr](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.repr) | 2 / 2 | 0 | 100.00% | available |
| [sympy.printing.rust](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.rust) | 0 / 5 | 5 | 0.00% | available |
| [sympy.printing.smtlib](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.smtlib) | 2 / 2 | 0 | 100.00% | available |
| [sympy.printing.str](https://docs.sympy.org/latest/modules/printing.html#module-sympy.printing.str) | 2 / 4 | 2 | 50.00% | available |
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
| `sympy.sets` | 21 / 26 | 5 | 80.77% | available |
| [sympy.sets.conditionset](https://docs.sympy.org/latest/modules/sets.html#module-sympy.sets.conditionset) | 0 / 3 | 3 | 0.00% | missing_module |
| [sympy.sets.fancysets](https://docs.sympy.org/latest/modules/sets.html#module-sympy.sets.fancysets) | 0 / 12 | 12 | 0.00% | missing_module |
| [sympy.sets.powerset](https://docs.sympy.org/latest/modules/sets.html#module-sympy.sets.powerset) | 0 / 1 | 1 | 0.00% | missing_module |
| [sympy.sets.sets](https://docs.sympy.org/latest/modules/sets.html#module-sympy.sets.sets) | 0 / 22 | 22 | 0.00% | missing_module |

</details>

<details>
<summary>simplify</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| `sympy.simplify` | 0 / 32 | 32 | 0.00% | missing_module |
| [sympy.simplify.combsimp](https://docs.sympy.org/latest/modules/simplify/simplify.html#module-sympy.simplify.combsimp) | 0 / 1 | 1 | 0.00% | missing_module |
| [sympy.simplify.cse_main](https://docs.sympy.org/latest/modules/simplify/simplify.html#module-sympy.simplify.cse_main) | 0 / 11 | 11 | 0.00% | missing_module |
| [sympy.simplify.epathtools](https://docs.sympy.org/latest/modules/simplify/simplify.html#module-sympy.simplify.epathtools) | 0 / 2 | 2 | 0.00% | missing_module |
| [sympy.simplify.fu](https://docs.sympy.org/latest/modules/simplify/fu.html#module-sympy.simplify.fu) | 0 / 31 | 31 | 0.00% | missing_module |
| [sympy.simplify.hyperexpand](https://docs.sympy.org/latest/modules/simplify/simplify.html#module-sympy.simplify.hyperexpand) | 0 / 37 | 37 | 0.00% | missing_module |
| [sympy.simplify.hyperexpand_doc](https://docs.sympy.org/latest/modules/simplify/hyperexpand.html#module-sympy.simplify.hyperexpand_doc) | 0 / 0 | 0 | N/A | missing_module |
| [sympy.simplify.powsimp](https://docs.sympy.org/latest/modules/simplify/simplify.html#module-sympy.simplify.powsimp) | 0 / 2 | 2 | 0.00% | missing_module |
| [sympy.simplify.radsimp](https://docs.sympy.org/latest/modules/simplify/simplify.html#module-sympy.simplify.radsimp) | 0 / 17 | 17 | 0.00% | missing_module |
| [sympy.simplify.ratsimp](https://docs.sympy.org/latest/modules/simplify/simplify.html#module-sympy.simplify.ratsimp) | 0 / 2 | 2 | 0.00% | missing_module |
| [sympy.simplify.simplify](https://docs.sympy.org/latest/explanation/gotchas.html#module-sympy.simplify.simplify) | 0 / 21 | 21 | 0.00% | missing_module |
| [sympy.simplify.sqrtdenest](https://docs.sympy.org/latest/modules/simplify/simplify.html#module-sympy.simplify.sqrtdenest) | 0 / 6 | 6 | 0.00% | missing_module |
| [sympy.simplify.trigsimp](https://docs.sympy.org/latest/modules/simplify/simplify.html#module-sympy.simplify.trigsimp) | 0 / 5 | 5 | 0.00% | missing_module |

</details>

<details>
<summary>solvers</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| [sympy.solvers](https://docs.sympy.org/latest/modules/solvers/solvers.html#module-sympy.solvers) | 29 / 46 | 17 | 63.04% | available |
| [sympy.solvers.inequalities](https://docs.sympy.org/latest/modules/solvers/inequalities.html#module-sympy.solvers.inequalities) | 8 / 8 | 0 | 100.00% | available |
| [sympy.solvers.ode](https://docs.sympy.org/latest/modules/solvers/ode.html#module-sympy.solvers.ode) | 11 / 11 | 0 | 100.00% | available |
| [sympy.solvers.ode.ode](https://docs.sympy.org/latest/modules/solvers/ode.html#module-sympy.solvers.ode.ode) | 4 / 23 | 19 | 17.39% | available |
| [sympy.solvers.pde](https://docs.sympy.org/latest/modules/solvers/pde.html#module-sympy.solvers.pde) | 9 / 9 | 0 | 100.00% | available |
| [sympy.solvers.recurr](https://docs.sympy.org/latest/modules/solvers/solvers.html#module-sympy.solvers.recurr) | 0 / 4 | 4 | 0.00% | missing_module |
| [sympy.solvers.simplex](https://docs.sympy.org/latest/modules/solvers/solvers.html#module-sympy.solvers.simplex) | 0 / 6 | 6 | 0.00% | missing_module |
| `sympy.solvers.solvers` | 15 / 15 | 0 | 100.00% | available |
| [sympy.solvers.solveset](https://docs.sympy.org/latest/modules/solvers/solveset.html#module-sympy.solvers.solveset) | 5 / 14 | 9 | 35.71% | available |

</details>

<details>
<summary>stats</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| [sympy.stats](https://docs.sympy.org/latest/modules/stats.html#module-sympy.stats) | 0 / 157 | 157 | 0.00% | missing_module |
| [sympy.stats.crv](https://docs.sympy.org/latest/modules/stats.html#module-sympy.stats.crv) | 0 / 9 | 9 | 0.00% | missing_module |
| [sympy.stats.crv_types](https://docs.sympy.org/latest/modules/stats.html#module-sympy.stats.crv_types) | 0 / 54 | 54 | 0.00% | missing_module |
| [sympy.stats.frv](https://docs.sympy.org/latest/modules/stats.html#module-sympy.stats.frv) | 0 / 9 | 9 | 0.00% | missing_module |
| [sympy.stats.frv_types](https://docs.sympy.org/latest/modules/stats.html#module-sympy.stats.frv_types) | 0 / 11 | 11 | 0.00% | missing_module |
| [sympy.stats.rv](https://docs.sympy.org/latest/modules/stats.html#module-sympy.stats.rv) | 0 / 41 | 41 | 0.00% | missing_module |

</details>

<details>
<summary>tensor</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| [sympy.tensor](https://docs.sympy.org/latest/modules/tensor/index.html#module-sympy.tensor) | 13 / 19 | 6 | 68.42% | available |
| [sympy.tensor.array](https://docs.sympy.org/latest/modules/tensor/array.html#module-sympy.tensor.array) | 16 / 16 | 0 | 100.00% | available |
| [sympy.tensor.array.expressions](https://docs.sympy.org/latest/modules/tensor/array_expressions.html#module-sympy.tensor.array.expressions) | 0 / 16 | 16 | 0.00% | missing_module |
| [sympy.tensor.index_methods](https://docs.sympy.org/latest/modules/tensor/index_methods.html#module-sympy.tensor.index_methods) | 0 / 3 | 3 | 0.00% | missing_module |
| [sympy.tensor.indexed](https://docs.sympy.org/latest/modules/tensor/indexed.html#module-sympy.tensor.indexed) | 0 / 4 | 4 | 0.00% | missing_module |
| [sympy.tensor.tensor](https://docs.sympy.org/latest/modules/tensor/tensor.html#module-sympy.tensor.tensor) | 35 / 35 | 0 | 100.00% | available |
| [sympy.tensor.toperators](https://docs.sympy.org/latest/modules/tensor/toperators.html#module-sympy.tensor.toperators) | 0 / 1 | 1 | 0.00% | missing_module |

</details>

<details>
<summary>testing</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| [sympy.testing](https://docs.sympy.org/latest/modules/testing/index.html#module-sympy.testing) | 0 / 2 | 2 | 0.00% | missing_module |
| [sympy.testing.pytest](https://docs.sympy.org/latest/modules/testing/pytest.html#module-sympy.testing.pytest) | 0 / 16 | 16 | 0.00% | missing_module |
| [sympy.testing.randtest](https://docs.sympy.org/latest/modules/testing/randtest.html#module-sympy.testing.randtest) | 0 / 0 | 0 | N/A | missing_module |
| [sympy.testing.runtests](https://docs.sympy.org/latest/modules/testing/runtests.html#module-sympy.testing.runtests) | 0 / 21 | 21 | 0.00% | missing_module |

</details>

<details>
<summary>utilities</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| [sympy.utilities](https://docs.sympy.org/latest/modules/utilities/index.html#module-sympy.utilities) | 19 / 25 | 6 | 76.00% | available |
| [sympy.utilities.autowrap](https://docs.sympy.org/latest/modules/utilities/autowrap.html#module-sympy.utilities.autowrap) | 0 / 9 | 9 | 0.00% | missing_module |
| [sympy.utilities.codegen](https://docs.sympy.org/latest/modules/utilities/codegen.html#module-sympy.utilities.codegen) | 0 / 16 | 16 | 0.00% | missing_module |
| [sympy.utilities.decorator](https://docs.sympy.org/latest/modules/utilities/decorator.html#module-sympy.utilities.decorator) | 3 / 10 | 7 | 30.00% | available |
| [sympy.utilities.enumerative](https://docs.sympy.org/latest/modules/utilities/enumerative.html#module-sympy.utilities.enumerative) | 0 / 6 | 6 | 0.00% | missing_module |
| [sympy.utilities.exceptions](https://docs.sympy.org/latest/modules/utilities/exceptions.html#module-sympy.utilities.exceptions) | 0 / 3 | 3 | 0.00% | missing_module |
| [sympy.utilities.iterables](https://docs.sympy.org/latest/modules/utilities/iterables.html#module-sympy.utilities.iterables) | 59 / 59 | 0 | 100.00% | available |
| [sympy.utilities.lambdify](https://docs.sympy.org/latest/modules/utilities/lambdify.html#module-sympy.utilities.lambdify) | 0 / 3 | 3 | 0.00% | missing_module |
| [sympy.utilities.memoization](https://docs.sympy.org/latest/modules/utilities/memoization.html#module-sympy.utilities.memoization) | 0 / 2 | 2 | 0.00% | missing_module |
| [sympy.utilities.misc](https://docs.sympy.org/latest/modules/utilities/misc.html#module-sympy.utilities.misc) | 1 / 13 | 12 | 7.69% | available |
| [sympy.utilities.source](https://docs.sympy.org/latest/modules/utilities/source.html#module-sympy.utilities.source) | 0 / 2 | 2 | 0.00% | missing_module |
| [sympy.utilities.timeutils](https://docs.sympy.org/latest/modules/utilities/timeutils.html#module-sympy.utilities.timeutils) | 1 / 2 | 1 | 50.00% | available |

</details>

<details>
<summary>vector</summary>

| Reference module | Present | Missing | API presence | Status |
|---|---:|---:|---:|---|
| [sympy.vector](https://docs.sympy.org/latest/modules/vector/index.html#module-sympy.vector) | 43 / 43 | 0 | 100.00% | available |

</details>


## Import errors

- `sympy.physics.continuum_mechanics.Arch`: `ModuleNotFoundError: No module named 'sympy.physics.continuum_mechanics.Arch'`
