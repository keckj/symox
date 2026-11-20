# Architecture

[README](../README.md)

Modules of `crates/symox-core/src` (crate name `symox`):

| module | content |
|---|---|
| `number.rs` | `Number`: `i64` / big integer / rational / `f64` / arbitrary-precision float |
| `bigfloat`, `evalf_mp` | arbitrary-precision floats and multiprecision evaluation (mpmath replacement) |
| `expr`, `symbol`, `host` | `Expr` (Arc'd node with cached hash and symbol bloom mask), `Kind`, symbols with assumptions, opaque handles for Python subclasses |
| `ops`, `functions`, `relational` | canonical constructors `add`, `mul`, `pow`, operators, function heads with automatic evaluation, relations and `Piecewise` |
| `assumptions`, `assumption_facts`, `ask` | three-valued assumption queries, fact tables, `Q`/`ask`/`refine` |
| `diff`, `expand`, `subs`, `evalf`, `traverse`, `pattern`, `sortkey` | core algorithms, matching and canonical ordering |
| `printer.rs`, `parser.rs`, `printing/` | `str`/`srepr`, `sympify`, LaTeX/pretty printing and code generation (C, CUDA, Rust, Python, JAX, PyTorch, Slang, SYCL, ...) |
| `lambdify`, `lambdify_vec`, `cephes_gamma` | numerical compilation and array evaluation |
| `poly/`, `ntheory` | exact polynomial arithmetic, factoring, Groebner bases, rational functions, Hermite and Smith normal forms, number theory |
| `matrix/`, `linear`, `simplex` | symbolic and exact matrix algorithms, domain-matrix kernels (LU, RREF, LLL), linear systems, exact rational simplex |
| `matrixexpr/`, `matrixcontraction`, `arrayexpr*` | typed matrix expressions and SymPy 1.14 array-expression ports (canonicalization, array to matrix, indexed conversions, derivatives) |
| `simp`, `simplify`, `simptools/` | simplification and structural helpers |
| `integrate/` | indefinite and definite integration, including Risch, rational integration and `meijerint` |
| `solve/`, `nsolve_mp` | `solve`, transcendental and radical solving, solution checks, `nsolve` |
| `series/`, `limits/` | series expansion and Gruntz limits |
| `summation`, `products`, `convergence`, `intlimits`, `recurr` | sums, products, convergence tests and recurrences |
| `ode/` | `dsolve` solvers and `classify_ode` matching |
| `logic/`, `sets/` | Boolean algebra, SAT, set operations and `solveset` |
| `special/` | special functions, numerical evaluation and `hyperexpand` |
| `stats`, `transforms/` | exact probability kernels and integral transforms |

Rust integration tests live in `crates/symox-core/tests` (for example
`lll_parity.rs` against SymPy 1.14 fixtures).

| path | content |
|---|---|
| `crates/symox-py` | PyO3 bindings, built as `symox._core` |
| `python/symox` | the Python package |
| `python/symox/extra/` | optional tensor/GPU printers, source generation and execution adapters |

The Python package exposes native nodes under SymPy's class names and method
resolution order (`_node_types.py`), so `isinstance` checks against `Float`,
`Number`, `Function` or `Basic` behave as in SymPy. User subclasses of
`Symbol`, `Function`, `Expr`, `Atom` or `Basic` are carried through the native
engine as opaque host handles (`host.rs`, `_host.py`).

Many Rust ports sit behind Python admission modules (for example
`polys/matrices/_*_native.py`): the module checks that the input is in the
ported domain, calls `symox._core`, and otherwise runs the retained Python
algorithm. Mathematical algorithms still implemented in Python are inventoried
in the [Python runtime audit](../audits/architecture/python-runtime-2026-10-28-r9/README.md)
and tracked in the [TODO](../TODO.md).

Invariants: ordinary `Expr` nodes are canonical. Only `ops`/`functions`
constructors create canonical nodes (`Expr::from_kind` is reserved for code
that guarantees the canonical form itself). The Python `evaluate=False`
constructors explicitly mark held nodes with `Expr::from_unevaluated_kind`.
Branch-aware evaluators can also retain a specific public expression form with
that marker, such as `im(sqrt(n)) = -I*sqrt(n)` for strictly negative real `n`.
`Add` stores `(term, coefficient)` pairs, `Mul`
stores `(base, exponent)` pairs, both sorted by `Expr`'s `Ord` (`Mul`
factors by base then exponent).
