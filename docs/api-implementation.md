# Internal APIs: behavior and limitations

[Full inventory](sympy-api.md) · [Status by section](api-connections.md) · [README](../README.md)

These implementations use Symox's native expressions. They do not import SymPy to perform calculations. A name in the inventory does not establish full compatibility with its SymPy protocol. The [algorithmic parity contract](algorithmic-parity.md) and [versioned audit](sympy-parity-audit.md) track behavioral evidence separately.

Native expressions do not imply native algorithms: the current Python package
also contains mathematical algorithms, including solving, integration,
domain-matrix operations and statistics. Porting these computations into Rust
while preserving their capabilities is outstanding architectural work, tracked
in the [active TODO](../TODO.md) and inventoried in the
[Python runtime audit](../audits/architecture/python-runtime-2026-10-28-r9/README.md). API coverage and Python test counts do not
establish compliance with the Rust-engine architecture stated in the README.

## Linear programming and logic

`solvers.simplex` implements an exact two-phase simplex method with Bland pivoting. `linprog` accepts inequalities, equalities and bounds; `lpmin` and `lpmax` convert symbolic relations and handle free variables. Infeasible and unbounded problems raise the corresponding exceptions. Rational pivoting runs in the native `simplex` module; other coefficient domains keep the Python path. Coefficients must support the exact comparisons required for pivoting; this is not a constraint solver for undecidable symbolic parameters.

`PropKB` stores CNF clauses and queries the local SAT solver. `pl_true` also handles partial assignments. These operations concern propositional logic.

## Polynomials and domain matrices

All five Appell polynomial families are constructed exactly, with an optional `Poly` result. `DomainMatrix`, DDM, SDM and the `ddm_*`/`sdm_*` functions follow the SymPy 1.14 implementation on native domain elements. Selected kernels run in Rust behind admission modules (`polys/matrices/_*_native.py`): fraction-free LU over ZZ/QQ, LU and LU solve over QQ, RREF over ZZ/QQ, sparse QQ Gauss-Jordan, Hermite normal form (including modulo D), Smith normal form and LLL. Inputs these modules decline run the retained Python algorithm.

**DFM requires FLINT, as in SymPy without it.** Symox reports `GROUND_TYPES='python'`, so `DFM` is the reference placeholder and `to_dfm()` raises `NotImplementedError("DFM requires GROUND_TYPES=flint.")`.

`polyconfig` provides context-local configuration, defaults and temporary overrides. Algorithm selectors do not control every computation path in the Rust engine.

Module morphisms use exact matrices for images, kernels and composition; minimal polynomials use Krylov powers. Finite searches modulo an integer are bounded. These interfaces do not replace all of SymPy's ideal algorithms or arbitrary-rank module algorithms.

## Geometry, matrices and sets

The relevant geometric loci inherit from `GeometrySet` and `Set`. The 2D/3D line classes have dedicated intermediate classes. Geometric operations use existing exact methods; inheriting from `Set` does not add a general algebra for every geometric set to the native engine.

`CompanionMatrix` constructs the companion matrix of a monic univariate polynomial. `MatrixSet` checks dimensions and membership of each entry. `MatrixKind` and `SetKind` provide interned metadata; the native matrix graph retains its existing internal discriminants.

Block matrix helpers reuse native simplifications and Schur complements. Set arithmetic distinguishes elementwise operations from union and intersection. It enumerates finite sets, applies supported interval rules and otherwise preserves a symbolic image.

## Code generation and compilation

The AST includes types, variables, declarations, assignments, loops, calls, functions and returns, plus C- and Fortran-specific nodes. Printers and `utilities.codegen` emit C, Fortran, Python, Rust, Julia and Octave for supported constructs. Tests compile and execute C, Fortran and Rust routines, and execute generated Python.

Numerical optimizations include `expm1`, `log1p`, powers and transformations of exponential sums. Symbolic intrinsics retain their intent for the printer; their derivatives and rewrites do not reproduce every internal SymPy protocol. Polynomial approximations remain subject to their methods' bounds and checks.

`autowrap(..., backend='f2py')` compiles an extension with NumPy f2py and a Fortran compiler. Default builds use Meson when installed with a compatible NumPy; explicit compiler/backend flags remain delegated to f2py. The Cython backend requires an optional Cython installation. `ufuncify` compiles a C extension providing a real `numpy.ufunc`, with NumPy broadcasting and multiple outputs in double precision. Supported scalar signatures do not provide a general adaptation of indexed arrays and complex types. Compilation errors include compiler diagnostics.

`implemented_function` combines Python functions with native numerical expressions, including nested calls. `lambdastr` produces executable lambdas with nested arguments. Custom functions must provide an implementation compatible with the selected numerical backend.

Symox-specific PyTorch, JAX, Slang, SYCL/AdaptiveCpp and CUDA extensions live in `symox.extra`. Their source printers do not require execution libraries. Compatible SymPy entry points remain available; see [code generation](codegen/overview.md) and [test environments](testing.md).

## Abstract arrays and tensor derivatives

`ArraySymbol`, tensor products, contractions, diagonals, permutations, sums and reshapes retain dimensions without materializing entries. Explicit operations reuse native arrays and their 100,000-entry limit. Contractions of symbolic size remain exact sums.

Abstract entries use symbolic `ArrayElement` applications to retain index dependencies in native calculations. Their internal representation differs from SymPy's `Indexed` objects. Canonicalization, `convert_array_to_matrix`, `convert_indexed_to_array`, `convert_array_to_indexed` and `array_derive` are Rust ports of the SymPy 1.14 algorithms. Inputs the native code declines run the vendored, unmodified SymPy 1.14 implementation (`symox._arrayexpr`) on a structural translation.

`PartialDerivative` computes free and contracted indices and evaluates derivatives from supplied components. Differentiation components must be differentiable coordinates; symbolic product algebra involving these new objects is not fully adapted.

## Probability spaces

Domains are separate from their measures. Simple, product and conditional finite or continuous domains provide predicates and exact integration or enumeration. `pspace` returns specialized spaces for supported scalar random variables.

General finite spaces retain joint weights for probabilities, expectations and sampling. General continuous spaces use their joint density and can retain unevaluated integrals. Conditional sampling by inversion can be expensive; it requires interval supports and numerical parameters. The general restrictions in the [statistics guide](statistics.md) still apply.

## Core, assumptions and numerical evaluation

Number classes, function types, `Kind` metadata, singletons and containers identify or construct the corresponding native objects. `StdFactKB` derives assumptions through the native engine; custom binary relations can register `ask` handlers. Factor, traversal, symbol and cache helpers operate on local expressions.

The `core.evalf` helpers expose mpmath's tuple protocol with precision in bits. Evaluation, sums and quadrature use the available numerical engines. Precision fields are working estimates, not general interval error certificates. Exact `Float` equality accounts for binary values; `all_close` provides a separate tolerance comparison.

Real and imaginary projections retain unknown finiteness for unconstrained arguments. `Abs(x+oo)` therefore remains symbolic until assumptions or substitution justify evaluation. Principal half-integer powers of strictly negative real arguments have a provable imaginary part.

## Combinatorics, groups and tensors

Standard transitive groups of degrees 5 and 6, Rubik's Cube generators, Gray codes and RGS ranks have concrete representations. Permutation and glyph tables are static data; using them does not load SymPy. `find_transitive_subgroups_of_S6` checks catalog groups without reproducing an exhaustive randomized search.

`combinatorics.tensor_can` is a port of SymPy's Butler-Portugal implementation; `canonicalize_naive` enumerates groups as in the reference. `graph_certificate` uses `canonicalize` and additionally renumbers vertex labels that are not `0..n-1`. The tensor-level canonicalization in `symox.tensor` still explores signed orbits of slot symmetries and dummy indices, limited to 100,000 states.

Arithmetic criteria for counting cyclic, abelian or nilpotent groups are exact. General group counts remain limited to orders covered by available formulas and tables; other cases are rejected.

## ODE series and simplification

`dsolve` and `classify_ode` are adapted from SymPy 1.14's `ode.py`: hint machinery, hint order and solution forms follow the reference, and `classify_ode` matching runs in Rust. The earlier native solvers remain a fallback. Power-series hints (`2nd_power_series_ordinary`, `2nd_power_series_regular`) return the reference output; for a regular singular point with a double or resonant pair of indicial roots, Symox additionally appends the logarithmic second solution, which is an extension beyond SymPy.

Nonlinear systems use the reference `sysode_nonlinear_2eq_order1` (types 1 to 5) and `3eq_order1` solvers. The second-order nonlinear inspectors return `None`, like their reference counterparts. See the [differential equations guide](differential-equations.md) for current gaps.

CSE helpers order dependencies, track common arguments and release temporaries after their last use. Sum and product combination handles identical or adjacent limits. `nthroot` accepts a candidate algebraic simplification after exactly verifying its power. `inversecombine` is a formal transformation that deliberately ignores branch restrictions.

`trigsimp_groebner` is a port of SymPy's algorithm: an ideal of trigonometric relations, a Groebner basis over ZZ and `ratsimpmodprime`, with the `hints` and `polynomial` options. When the reference raises `ValueError` for rational multiples of an angle, Symox falls back to its earlier common-angle reduction (at most eight distinct arguments, integer multipliers up to 16); `trigsimp(..., method='matching')` also uses that reduction.

`hyperexpand` is a Rust port (`special/hyperexpand`) of the SymPy 1.14 formula tables, contiguous shift operators and `devise_plan` bucket-by-bucket planning; `devise_plan_meijer` follows the same reference. Remaining bounds of the Python layer: order reduction up to 256 factors, terminating series up to 10,000 terms, shifted sums up to 256 terms, `Formula` bases of order 1 to 32, and at most eight parameters per group for matching. Inputs outside the ported catalog stay unevaluated.

## Polynomial integration and program parsing

`polytope_integrate` computes moments through boundary integrals and facet triangulation, using exact arithmetic in two and three dimensions. It accepts oriented polygons, faceted polyhedra and bounded planar half-space intersections. Maximum degree is 64. Orientation determines the sign; inputs must describe a consistent domain and facets.

`meijerint_indefinite`, `meijerint_definite` and `meijerint_inversion` are Rust ports of SymPy's Meijer G integration (`integrate/meijerint`: table, conditions, recursive and inversion steps). They return `None` when the method does not apply, keep convergence conditions explicit and raise the reference errors.

`SymPyExpression` parses a numerical subset of C and Fortran: recognized declarations, assignments, expressions, functions and loops. The AST can emit C, Fortran or Python. This parser does not implement the full grammars, preprocessor, linker, or all typing and scoping rules. Programs outside the subset must be adapted or processed by a dedicated compiler.

## Printing, sessions and test runners

Naming conventions, precedence tables, baseline-aware pretty layouts and language mappings connect to the printers. Symbolic expressions retain each backend's rendering limitations. The Rust printer handles integer scalar declarations and calculations, loop ranges and integer-to-float conversions in tested functions; it is not a general Rust type checker.

Python sessions can transform integer literals and create unknown loaded symbols automatically. IPython sessions use transformation hooks when IPython is installed. `AesaraPrinter` builds Aesara graphs when a compatible package is available; it does not install that package and otherwise fails explicitly.

`testing` reuses pytest and Python doctests. These interfaces execute real tests and handle supported assertions, warnings and subprocesses without embedding SymPy's own collection and rewriting engine. Unadapted internal selectors and events do not expose every feature of the SymPy runner.

## Validation

The `tests/test_remaining_*` files check exact results, dependencies and dimensions, errors, independent numerical values of special operators and generated code execution. Code tests compile C, Fortran and Rust programs; some require optional tools.

The curated 273-expression differential corpus remains separate from API name presence. The generated inventory measures neither every signature nor every SymPy behavior. Current test results and immutable audit checkpoints are recorded in the [work log](../audits/sympy-1.14.0/WORKLOG.md); previous counts describe their own snapshots and are not current parity certificates.
