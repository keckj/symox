# SymPy API implementation by section

[Measured inventory](sympy-api.md) · [Module behavior and limitations](api-implementation.md) · [README](../README.md)

Measured on 2026-10-19 against SymPy 1.14.0. The initial audit remains available in commit `e2a1816`.

310 module paths: 301 importable with measurable exports, 0 missing, 9 without measurable exports and 0 reference import errors. No measured rows have missing names. The incorrect upstream `Arch` module entry is recorded as a documentation correction and measured as the `Arch` class in `continuum_mechanics.arch`.

These counts measure name presence. They do not verify signatures or behavior. Rows overlap and cannot be added into global coverage. Paths without exports and genuine reference errors remain explicit. See the [algorithmic parity contract](algorithmic-parity.md) for the stronger acceptance criteria.

| Section | Paths | Missing paths | Rows with missing names |
|---|---:|---:|---:|
| [root](#root) | 1 | 0 | 0 |
| [abc](#abc) | 1 | 0 | 0 |
| [algebras](#algebras) | 1 | 0 | 0 |
| [assumptions](#assumptions) | 5 | 0 | 0 |
| [calculus](#calculus) | 5 | 0 | 0 |
| [categories](#categories) | 2 | 0 | 0 |
| [codegen](#codegen) | 13 | 0 | 0 |
| [combinatorics](#combinatorics) | 16 | 0 | 0 |
| [core](#core) | 23 | 0 | 0 |
| [crypto](#crypto) | 2 | 0 | 0 |
| [diffgeom](#diffgeom) | 1 | 0 | 0 |
| [discrete](#discrete) | 3 | 0 | 0 |
| [functions](#functions) | 14 | 0 | 0 |
| [geometry](#geometry) | 9 | 0 | 0 |
| [holonomic](#holonomic) | 1 | 0 | 0 |
| [integrals](#integrals) | 5 | 0 | 0 |
| [interactive](#interactive) | 3 | 0 | 0 |
| [liealgebras](#liealgebras) | 1 | 0 | 0 |
| [logic](#logic) | 2 | 0 | 0 |
| [matrices](#matrices) | 8 | 0 | 0 |
| [ntheory](#ntheory) | 14 | 0 | 0 |
| [parsing](#parsing) | 3 | 0 | 0 |
| [physics](#physics) | 67 | 0 | 0 |
| [plotting](#plotting) | 3 | 0 | 0 |
| [polys](#polys) | 17 | 0 | 0 |
| [printing](#printing) | 31 | 0 | 0 |
| [series](#series) | 2 | 0 | 0 |
| [sets](#sets) | 5 | 0 | 0 |
| [simplify](#simplify) | 13 | 0 | 0 |
| [solvers](#solvers) | 9 | 0 | 0 |
| [stats](#stats) | 6 | 0 | 0 |
| [tensor](#tensor) | 7 | 0 | 0 |
| [testing](#testing) | 4 | 0 | 0 |
| [utilities](#utilities) | 12 | 0 | 0 |
| [vector](#vector) | 1 | 0 | 0 |

<a id="root"></a>
<details>
<summary>root</summary>

All 921 top-level exports are present. This does not guarantee compatible signatures or behavior.

No missing names in the measured paths.

</details>

<a id="abc"></a>
<details>
<summary>abc</summary>

The module provides 76 predefined Latin/Greek symbols and the `_clash` dictionaries. The 0/0 count does not measure these constants: the path exists, but N/A remains N/A.

No missing names in the measured paths.

</details>

<a id="algebras"></a>
<details>
<summary>algebras</summary>

No missing connections in the inventoried surface.

No missing names in the measured paths.

</details>

<a id="assumptions"></a>
<details>
<summary>assumptions</summary>

Extensible predicates and binary relations, CNF facts and refine helpers connect to the native engine. Available rules remain those of the engine and registered handlers.

No missing names in the measured paths.

</details>

<a id="calculus"></a>
<details>
<summary>calculus</summary>

All 17 measured exports are connected, including `AccumBounds` and `euler_equations`. Measured submodules also have complete name presence.

No missing names in the measured paths.

</details>

<a id="categories"></a>
<details>
<summary>categories</summary>

No missing connections in the two measured paths.

No missing names in the measured paths.

</details>

<a id="codegen"></a>
<details>
<summary>codegen</summary>

Program ASTs, types, declarations, loops, functions and C/Fortran nodes are available. Numerical rewrites, Newton methods and C, Fortran, Python, Rust, Julia and Octave generators are connected. See the behavior and limitations guide. Custom PyTorch, JAX, Slang, SYCL/AdaptiveCpp and CUDA generators live in `symox.extra.codegen`; `symox.extra` exports their entry points.

No missing names in the measured paths.

</details>

<a id="combinatorics"></a>
<details>
<summary>combinatorics</summary>

Transitive group constructors of degrees 5 and 6, group counts, Gray/RGS and stabilizer helpers are available. Tensor canonicalization uses bounded exact exploration of signed orbits.

No missing names in the measured paths.

</details>

<a id="core"></a>
<details>
<summary>core</summary>

Types and metadata, assumptions, singletons, factors, traversal and numerical helpers are connected. Noncommutative symbols and undefined functions (`commutative=False`) are supported. `match`, `matches` and `replace` use a port of SymPy's `matches` protocol, including Wild `properties`. The evalf precision tuples adapt the native engine; they are not interval certificates.

No missing names in the measured paths.

</details>

<a id="crypto"></a>
<details>
<summary>crypto</summary>

No missing connections in the two measured paths.

No missing names in the measured paths.

</details>

<a id="diffgeom"></a>
<details>
<summary>diffgeom</summary>

No missing connections in the inventoried surface.

No missing names in the measured paths.

</details>

<a id="discrete"></a>
<details>
<summary>discrete</summary>

No missing connections in the three measured paths.

No missing names in the measured paths.

</details>

<a id="functions"></a>
<details>
<summary>functions</summary>

Function-family classes connect to native expressions, with real decomposition for supported functions and numerical incomplete-beta corrections.

No missing names in the measured paths.

</details>

<a id="geometry"></a>
<details>
<summary>geometry</summary>

`GeometrySet` and the intermediate `LinearEntity2D`/`LinearEntity3D` classes are integrated into the hierarchy. Rotation, translation and scaling matrices, symbol lookup and coplanarity are implemented.

No missing names in the measured paths.

</details>

<a id="holonomic"></a>
<details>
<summary>holonomic</summary>

No missing connections in the inventoried surface.

No missing names in the measured paths.

</details>

<a id="integrals"></a>
<details>
<summary>integrals</summary>

Exact polynomial integration over polygons and polyhedra, moments and facet helpers are available. Meijer adapters use native integration and preserve unresolved failures or conditions.

No missing names in the measured paths.

</details>

<a id="interactive"></a>
<details>
<summary>interactive</summary>

Integer transformation and automatic symbol creation connect to Python/IPython sessions. IPython hooks require its optional installation.

No missing names in the measured paths.

</details>

<a id="liealgebras"></a>
<details>
<summary>liealgebras</summary>

No missing connections in the inventoried surface.

No missing names in the measured paths.

</details>

<a id="logic"></a>
<details>
<summary>logic</summary>

Propositional knowledge bases use CNF conversion and the local SAT solver. `gateinputcount` and partial logical evaluation are available.

No missing names in the measured paths.

</details>

<a id="matrices"></a>
<details>
<summary>matrices</summary>

`MatrixKind`, `MatrixSet`, `CompanionMatrix` and block matrix rewrites are connected. Representation classes adapt native matrices; the new objects do not replace the native engine.

No missing names in the measured paths.

</details>

<a id="ntheory"></a>
<details>
<summary>ntheory</summary>

`digits` connects to `ntheory.factor_`. No names are missing in measured paths.

No missing names in the measured paths.

</details>

<a id="parsing"></a>
<details>
<summary>parsing</summary>

`SymPyExpression` parses an explicit subset of numerical C and Fortran programs and generates an executable AST. Unsupported grammar constructs are rejected.

No missing names in the measured paths.

</details>

<a id="physics"></a>
<details>
<summary>physics</summary>

No names are missing in the measured importable paths. The upstream inventory incorrectly labels the `Arch` class as a module. Its correct path is `symox.physics.continuum_mechanics.arch.Arch`, which is available. Name presence does not remove the limitations described in the physics guide.

No missing names in the measured paths.

</details>

<a id="plotting"></a>
<details>
<summary>plotting</summary>

Argument validation and plot construction connect to existing Plot objects. Graphics dependencies remain optional.

No missing names in the measured paths.

</details>

<a id="polys"></a>
<details>
<summary>polys</summary>

Exact Appell polynomials, context-local configuration, DDM/SDM/DFM storage, dense/sparse computation kernels, polynomial linear systems and module morphisms are connected. DFM uses the Symox engine; no FLINT backend is claimed.

No missing names in the measured paths.

</details>

<a id="printing"></a>
<details>
<summary>printing</summary>

Python/NumPy/mpmath, C89/C99/C11, Fortran and Rust printers include program-node rendering. Pretty protocols, conventions, precedence and other language mappings connect to the printers. Aesara requires its optional package. Custom PyTorch, JAX, Slang and SYCL/AdaptiveCpp printers live in `symox.extra.printing`, and the native `cudacode` printer is exported from `symox.extra`; compatible SymPy printing paths remain available.

No missing names in the measured paths.

</details>

<a id="series"></a>
<details>
<summary>series</summary>

No missing connections in the two inventoried paths. The report does not cover every internal submodule.

No missing names in the measured paths.

</details>

<a id="sets"></a>
<details>
<summary>sets</summary>

`SetKind`, Cartesian/polar complex regions and angle normalization are available. Elementwise arithmetic uses interval rules, finite enumeration or symbolic images.

No missing names in the measured paths.

</details>

<a id="simplify"></a>
<details>
<summary>simplify</summary>

CSE helpers, sum/product combination, radicals and Gröbner trigonometric reduction are connected. Hypergeometric and Meijer operators use actual differential polynomials; the formula catalog remains limited.

No missing names in the measured paths.

</details>

<a id="solvers"></a>
<details>
<summary>solvers</summary>

An exact two-phase linear programming solver supports equality constraints, bounds and free variables in `simplex`. Solveset helpers, constants and initial conditions, ODE series recurrences and triangular nonlinear systems connect to the APIs. Frobenius series and coupled-system limitations are documented.

No missing names in the measured paths.

</details>

<a id="stats"></a>
<details>
<summary>stats</summary>

Finite/continuous domains and probability spaces connect to `pspace` and expectation/probability engines. Finite joint distributions retain dependencies during calculations and sampling. The statistics guide limitations still apply.

No missing names in the measured paths.

</details>

<a id="tensor"></a>
<details>
<summary>tensor</summary>

Abstract array expressions, matrix/indexed conversions, axis operations and component derivatives are connected. Symbolic sizes remain lazy; materialization is bounded. General indexed conversions remain limited.

No missing names in the measured paths.

</details>

<a id="testing"></a>
<details>
<summary>testing</summary>

Executable pytest, doctest and randomized-test adapters use the available Python runners. They do not reproduce SymPy's entire internal testing infrastructure.

No missing names in the measured paths.

</details>

<a id="utilities"></a>
<details>
<summary>utilities</summary>

Memoization, warnings, source and string utilities are connected. `lambdastr` and `implemented_function` connect to numerical calculations. `codegen`, `autowrap`/f2py and NumPy `ufuncify` compile real routines; Cython remains optional.

No missing names in the measured paths.

</details>

<a id="vector"></a>
<details>
<summary>vector</summary>

No missing connections in the inventoried surface.

No missing names in the measured paths.

</details>
