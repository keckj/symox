# Public core, utilities and sessions

[README](../README.md)

All 921 names in SymPy 1.14.0's root `__all__` are exposed by symox. This is an
import-surface measurement. It does not mean every signature, class hierarchy,
method, mathematical domain or algorithm is compatible. The
[full inventory](sympy-api.md) continues to show incomplete and absent modules.

## Core interfaces

`Dict` is an immutable Python mapping of native symbolic keys and values, with
symbolic substitution, replacement and free-symbol queries. `RealNumber` is the
native-backed `Float` constructor. `FunctionClass` supports checks of public
function-head factories. `WildFunction` filters native wildcard matching to
applied functions with optional allowed argument counts; these constructors do
not recreate SymPy's Python inheritance hierarchy.

`evaluate(False)` changes Python scalar arithmetic, `Add`, `Mul`, `Pow`, ordinary
function wrappers, relations and `Piecewise` construction to retain held nodes.
The context is nested and uses `ContextVar`, so state restores on exceptions and
is isolated between independent contexts. An explicit `evaluate=True` overrides
the context on supported constructors. Internal Rust algorithms still perform
canonical arithmetic. `seterr(divide=True)` makes Python native-scalar `0/0`
raise; it does not alter all internal Rust numerical algorithms.

`SympifyError` wraps public sympification failures. Low-level `evalf(expr,prec,
options)` accepts binary precision and returns mpmath-compatible mpf tuples and
precision fields; `PrecisionExhausted` reports a requested precision above
`maxprec`. This adapter does not implement SymPy's complete adaptive error
tracking or quadrature options. High-level `N` and `Expr.evalf` retain their
usual decimal-precision interfaces.

`Options(gens,args,flags=None,strict=False)` validates supported polynomial
options and normalizes generators, simple coefficient-domain names, modulus
and monomial order. Legacy `ZZ_python`, `ZZ_gmpy`, `QQ_python`, `QQ_gmpy`,
`FF_python`, `FF_gmpy` and their domain class names construct native domains.
The backend suffixes do not select GMPY computations. `PythonRational` is an
exact Python rational interoperating with native coefficient domains.

## Matrices and indexed expressions

`DeferredVector` creates lazily named scalar components and supports numerical
argument packing in `lambdify`. `Permanent` holds a matrix permanent; `per`
evaluates concrete matrices with the native algorithm. `blockcut` constructs
validated matrix blocks, `diagonalize_vector` constructs diagonal matrices, and
`matrix_symbols` finds abstract matrix symbols recursively.

`Idx` stores an integer label and optional bounds. `get_indices` checks Einstein
free-index conformance; `get_contraction_structure` identifies contractions and
nested contractions. The current scalar code-generation architecture stores
indexed references as registered native symbol atoms. Shape-less references
retain their displayed indices, while shaped references follow the existing
flattened C/C++ index representation. This is not a new native indexed DAG:
registration metadata is local to the Python process, and using the same index
label with conflicting bounds is outside the supported model.

## Utilities, assumptions and dispatch

`arity`, `comp`, `capture`, `filldedent`, `public`, `threaded`, `xthreaded`,
`vectorize` and `timed` provide callable arity inspection, numeric comparison,
stdout capture, formatting, export registration, elementwise operations and
adaptive timing. `vectorize` maps selected arguments through nested sequences;
selecting several sequence arguments produces their Cartesian application.
`threaded` also distributes over sums; `xthreaded` preserves the outer sum.

`check_assumptions` and `failing_assumptions` retain three-valued native query
results. `register_handler` and `remove_handler` add callable or class-based
handlers for unresolved predicates. Contradictory handler answers raise.
`bool_map` searches for a variable bijection and checks Boolean equivalence
with the native SAT solver; its permutation search is intended for small
Boolean expressions. `expand_multinomial` selects the native multinomial
expansion mode.

`multipledispatch.dispatch` and `Dispatcher` resolve type signatures by subclass
specificity, cache resolutions and reject ambiguities. The `strategies` package
provides deterministic strategy composition, fixed points, conditional rewrites,
identity removal and grouping. These packages expose useful subsets rather than
the entire corresponding SymPy helper-module trees.

## Interactive and optional modules

`init_printing` configures Python's display hook and optionally IPython's plain
and LaTeX formatters. `init_session` explicitly launches a console with common
symbols and functions; automatic source transformations remain unsupported.
Importing symox never launches a console. `plot_backends` exposes the existing
plot backend registry.

`test` invokes pytest on supplied paths, or the source checkout's test directory;
`doctest` executes supplied modules or text files and reports actual failures.
pytest is needed only when calling the test runner. `concrete` exposes native
sums/products; `external.import_module` loads optional dependencies, with
optional version checks; `release` provides symox's own release version.
