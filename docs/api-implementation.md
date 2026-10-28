# Internal APIs: behavior and limitations

[Full inventory](sympy-api.md) · [Status by section](api-connections.md) · [README](../README.md)

These implementations use Symox's native expressions. They do not import SymPy to perform calculations. A name in the inventory does not establish full compatibility with its SymPy protocol. The [algorithmic parity contract](algorithmic-parity.md) and [versioned audit](sympy-parity-audit.md) track behavioral evidence separately.

Native expressions do not imply native algorithms: the current Python package
also contains mathematical algorithms, including solving, integration,
domain-matrix operations and statistics. Porting these computations into Rust
while preserving their capabilities is outstanding architectural work, tracked
in the [active TODO](../TODO.md). API coverage and Python test counts do not
establish compliance with the Rust-engine architecture stated in the README.

## Linear programming and logic

`solvers.simplex` implements an exact two-phase simplex method with Bland pivoting. `linprog` accepts inequalities, equalities and bounds; `lpmin` and `lpmax` convert symbolic relations and handle free variables. Infeasible and unbounded problems raise the corresponding exceptions. Coefficients must support the exact comparisons required for pivoting; this is not a constraint solver for undecidable symbolic parameters.

`PropKB` stores CNF clauses and queries the local SAT solver. `pl_true` also handles partial assignments. These operations concern propositional logic.

## Polynomials and domain matrices

All five Appell polynomial families are constructed exactly, with an optional `Poly` result. DDM, SDM and the `ddm_*`/`sdm_*` functions operate on actual dense or sparse representations and reuse `DomainMatrix` for exact algorithms.

**DFM uses the Symox engine without FLINT.** This container does not provide all optimizations or protocols of SymPy's FLINT backend. Field and ring algorithms remain distinct; not every storage method outside the inventory has been adapted.

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

Abstract entries use symbolic `ArrayElement` applications to retain index dependencies in native calculations. Their internal representation differs from SymPy's `Indexed` objects. Conversions cover recognized matrix sums and products, transposes, Kronecker products and direct indexed accesses; they do not recognize every arbitrary Einstein expression.

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

`tensor_can` explores signed orbits of slot symmetries and dummy indices. It detects cancellation from opposite signs and selects a canonical representative. Exploration is limited to 100,000 states; this is not the optimized Butler–Portugal algorithm. `canonicalize_naive` uses the same orbit engine. The graph certificate uses bounded exhaustive search and differs in format from SymPy's tensor certificate.

Arithmetic criteria for counting cyclic, abelian or nilpotent groups are exact. General group counts remain limited to orders covered by available formulas and tables; other cases are rejected.

## ODE series and simplification

`dsolve` accepts first-order Taylor series, ordinary second-order linear series and homogeneous Frobenius series whose two branches can be computed by recurrence. Calculations are limited to 64 coefficients. Repeated indicial roots and resonances requiring logarithms raise `NotImplementedError`. Supported initial conditions are substituted into solutions and solved for constants. `classify_ode(..., n=...)` explicitly requests series classification; default classification retains exact solution families.

Triangular nonlinear systems of two or three functions reduce to scalar ODEs, with solution substitution and distinct constants. General coupled nonlinear systems remain outside this solver. Internal second-order inspectors recognize no nonlinear families, like their reference counterparts.

CSE helpers order dependencies, track common arguments and release temporaries after their last use. Sum and product combination handles identical or adjacent limits. `nthroot` accepts a candidate algebraic simplification after exactly verifying its power. `inversecombine` is a formal transformation that deliberately ignores branch restrictions.

The Gröbner trigonometric reduction uses polynomial relations for sin/cos or sinh/cosh and compares two generator orders, with at most eight distinct arguments. Commensurable angles reduce to a common variable when their integer multipliers do not exceed 16. Advanced `hints` and `polynomial` options do not yet select a strategy equivalent to SymPy's.

Hypergeometric and Meijer contiguous operators are actual polynomials in the Euler operator. Inverse shifts come from division by the affine operator in the differential equation; singular cancellations are rejected. Hypergeometric plans search for a monotone path avoiding these singularities, bounded by 256 steps and 10,000 states. Each parameter group has at most eight parameters for matching.

`Formula` constructs a differential basis and its companion matrix. Formula classes are Python adapters without SymPy's full `Expr` hierarchy. The elementary catalog and Lerch recognition are limited; search also uses expansions already known to the native engine. Terminating series and some sum shifts are computed exactly. These classes do not imply that SymPy's full hypergeometric catalog has been ported.

## Polynomial integration and program parsing

`polytope_integrate` computes moments through boundary integrals and facet triangulation, using exact arithmetic in two and three dimensions. It accepts oriented polygons, faceted polyhedra and bounded planar half-space intersections. Maximum degree is 64. Orientation determines the sign; inputs must describe a consistent domain and facets.

The `meijerint` adapters call native integration. An unresolved result returns `None`; known conditions remain explicit. They are not a second general Meijer G integration engine.

`SymPyExpression` parses a numerical subset of C and Fortran: recognized declarations, assignments, expressions, functions and loops. The AST can emit C, Fortran or Python. This parser does not implement the full grammars, preprocessor, linker, or all typing and scoping rules. Programs outside the subset must be adapted or processed by a dedicated compiler.

## Printing, sessions and test runners

Naming conventions, precedence tables, baseline-aware pretty layouts and language mappings connect to the printers. Symbolic expressions retain each backend's rendering limitations. The Rust printer handles integer scalar declarations and calculations, loop ranges and integer-to-float conversions in tested functions; it is not a general Rust type checker.

Python sessions can transform integer literals and create unknown loaded symbols automatically. IPython sessions use transformation hooks when IPython is installed. `AesaraPrinter` builds Aesara graphs when a compatible package is available; it does not install that package and otherwise fails explicitly.

`testing` reuses pytest and Python doctests. These interfaces execute real tests and handle supported assertions, warnings and subprocesses without embedding SymPy's own collection and rewriting engine. Unadapted internal selectors and events do not expose every feature of the SymPy runner.

## Validation

The `tests/test_remaining_*` files check exact results, dependencies and dimensions, errors, independent numerical values of special operators and generated code execution. Code tests compile C, Fortran and Rust programs; some require optional tools.

The curated 273-expression differential corpus remains separate from API name presence. The generated inventory measures neither every signature nor every SymPy behavior. Current test results and immutable audit checkpoints are recorded in the [work log](../audits/sympy-1.14.0/WORKLOG.md); previous counts describe their own snapshots and are not current parity certificates.
