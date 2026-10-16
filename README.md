# symox

**Fast symbolic mathematics in Rust, with a SymPy compatible Python API.**

symox is a computer algebra system written from scratch in Rust. It follows
SymPy's semantics (canonical forms, automatic simplifications, printing) so
that existing SymPy code can switch with minimal changes, while being one to
three orders of magnitude faster.

```python
from symox import symbols, sin, cos, exp, expand, diff, lambdify

x, y = symbols("x y")
expand((x + y)**3)            # x**3 + 3*x**2*y + 3*x*y**2 + y**3
diff(sin(x)*exp(x), x)        # exp(x)*sin(x) + exp(x)*cos(x)
(x**2 + y).subs(x, 2)         # y + 4
f = lambdify((x, y), sin(x)*cos(y))   # compiled, vectorised over numpy arrays
```

It is also a plain Rust library:

```rust
use symox::prelude::*;

let x = Expr::symbol("x");
let e = expand(&(&x + 1).powi(10));
println!("{}", diff(&e, &x));
```

## Status

Implemented (milestone 1, the core):

* expression DAG with cached hashes, canonical `Add`/`Mul`/`Pow` following
  SymPy's automatic simplification rules (like-term collection, power
  merging, exact radicals `sqrt(8) = 2*sqrt(2)`, `I**2 = -1`, numeric
  coefficient distribution, infinities and `nan`, ...)
* exact integers (machine word fast path, transparent big integers),
  rationals and floats
* elementary and special functions with automatic evaluation (special values
  of trig functions, parity, periodicity, inverse compositions, ...)
* `diff`, `expand` (sparse polynomial fast path), `subs`, `xreplace`,
  `evalf`, `free_symbols`, `has`, relationals, `Piecewise`
* `str`/`srepr` printers matching SymPy, a parser (`sympify`) accepting
  SymPy's `str` and `srepr` output
* `lambdify`: compiles expressions into a register machine with common
  subexpression sharing, vectorised over NumPy arrays (GIL released)
* Python bindings (PyO3, abi3 wheels for CPython >= 3.9), conversion to and
  from SymPy objects

Planned next: polynomials (`Poly`, `gcd`, `factor`, `cancel`, `apart`),
`solve`, `series`/`limit`, `simplify`, integration, LaTeX
printing, code generation (C/CUDA), arbitrary precision `evalf`,
assumptions.

Special functions include combinatorial numbers, gamma/zeta families, error
and exponential integrals, Bessel/Airy functions, elliptic integrals, Lambert W,
hypergeometric functions and orthogonal polynomials, with symbolic evaluation,
derivatives and double precision real/complex numerics. `expand_func`,
`besselsimp`, `hyperexpand`, Stirling numbers and multinomial coefficients are
also available. Unsupported symbolic cases remain unevaluated.

`special_nstr(expr, dps)` returns an arbitrary precision decimal string for
numeric expressions supported by the multiprecision evaluator, including
factorial, binomial, rising/falling factorials, beta, Catalan, Fibonacci and
Lucas via exact gamma/elementary identities. Other special function families
return `None` when multiprecision evaluation is unsupported; no double
precision fallback is presented as extra significant digits. The regular
`N`/`evalf` expression API still uses double precision.

## Compatibility with SymPy

`tests/compare_sympy.py` runs every case through SymPy and symox and
compares the printed canonical forms and, when they differ, numerical
equivalence. Known deliberate deviations:

* `diff(Abs(x), x)` returns `sign(x)` (the real-variable formula).
* `evalf` is double precision for now (`n` is accepted but ignored).

## Building

```bash
cargo test                      # Rust core
pip install maturin
maturin develop --release       # Python package into the active venv
python tests/compare_sympy.py   # differential tests (needs sympy)
python benches/bench_vs_sympy.py
```

## Symbolic matrices

`Matrix` provides mutable dense matrices; `ImmutableMatrix` and
`ImmutableDenseMatrix` provide immutable, hashable matrices. Constructors
accept nested sequences, NumPy arrays, SymPy matrices and nested matrix
blocks. Arithmetic, slicing, entry transformations, symbolic derivatives,
Jacobian/Hessian construction, exact determinant/inverse, row reduction,
linear solves, characteristic polynomials and supported eigenproblems run
in Rust. Cholesky, LDL, QR and diagonalization are also available.

```python
from symox import Matrix, ImmutableMatrix, symbols

x = symbols("x")
a = Matrix([[x, 1], [1, x]])
a.det()                          # x**2 - 1
a.inv()
a.gauss_jordan_solve(Matrix([1, 2]), freevar=True)
ImmutableMatrix(a)               # rejects mutation; usable as a dictionary key
```

`dot` supports Hermitian conjugation conventions; `diagonalize` supports
`reals_only`, `sort` and `normalize`. Matrix pickles preserve entry
assumptions and unevaluated arithmetic. Immutable results preserve their
immutability through arithmetic, slicing and decompositions; `as_mutable()`
returns an independent editable copy.

This API does not cover every SymPy matrix method. Eigenproblems beyond the
supported exact algebraic cases, custom zero/simplification callbacks and
unsupported constructor/algorithm options raise `NotImplementedError`.
Matrix `evalf` currently supports its default 15-digit double precision;
other precision requests raise `NotImplementedError`. Vector norms and
matrix Frobenius/1/infinity norms are supported; matrix spectral norms
are not implemented.

```bash
python tests/compare_matrices.py
python -m pytest tests/test_matrix_compat.py tests/compare_matrices.py
```

## Symbolic code generation compatibility

For the scalar polynomial API used by TriSplyne, replace imports consistently:

```python
import symox as sm
from symox.printing.cxx import CXX17CodePrinter
from symox.printing.precedence import precedence, PRECEDENCE
from symox.codegen.ast import Assignment, CodeBlock, Comment
```

`Add`, `Mul`, `Pow`, `Symbol`, `Integer`, `Rational`, `Float`, and
`UnevaluatedExpr` provide constructors and `isinstance` checks for native Rust
expressions. `S.Half`, tree reconstruction through `expr.func(*expr.args)`,
coefficient extraction, operation counts, deterministic `default_sort_key`,
NumPy object arrays, and expression pickling are supported. Free symbol
discovery is cached on immutable expression nodes; callers receive independent
Python sets.

`UnevaluatedExpr` and `Mul(..., evaluate=False)` preserve coefficients and
grouping for Horner code generation. `subs` supports algebraic sum and product
patterns; `xreplace` performs exact replacement. `reduced(f, divisors, *gens)`
returns polynomial quotients and a remainder, with deterministic generator
ordering when generators are inferred. Reduction currently supports `order="lex"`;
explicit `domain`/`modulus`, `polys=True`, and `auto=False` raise
`NotImplementedError`.

The extensible Python `CXX17CodePrinter` supports `_print_<Class>` overrides,
`doprint`, and `parenthesize`. The existing top-level `cxxcode` uses the native
Rust printer. `StrPrinter` and `PrettyPrinter` are available for debug output.
This compatibility surface targets scalar polynomial workflows; it is not a
replacement for every SymPy module.

Consumer integration tests load the real TriSplyne modules with each backend
and compare generated text. They need a sibling checkout or an explicit path:

```bash
TRISPLYNE_PATH=/path/to/trisplyne python -m pytest tests/test_trisplyne_integration.py
```

## Differential fuzzing

The fuzzer builds the same seeded expression trees through both APIs. Each
native worker is isolated in a subprocess and restarted after a crash or case
timeout. A wall-clock budget also cancels shrinking; `--jobs` controls the
number of independent workers. Findings are shrunk and deduplicated into
`findings.jsonl`, with run counts in `stats.json`.

```bash
python tests/fuzz/fuzz_core.py --mode core --seeds 0:100 --jobs 2 --time 60
python tests/fuzz/fuzz_core.py --replay core 7 srepr --timeout 5
python tests/fuzz/triage.py tests/fuzz/out/findings.jsonl wrong,crash,died
```

Modes also cover integration, solving, series, limits, summation,
simplification, lambdify, and parsing (`--mode all` runs all modes). Findings
distinguish numerical errors, canonical-form differences, unsupported results,
crashes, timeouts, and slow calls. Runs report findings without failing by
default; use `--fail-on wrong,crash,died` for an explicit CI failure policy.
Replay uses the same process isolation and timeout as a normal run.

## Architecture

| module | content |
|---|---|
| `number.rs` | `Number`: `i64` / big integer / rational / `f64` |
| `expr.rs` | `Expr` (Arc'd node with cached hash and symbol bloom mask), `Kind` |
| `ops.rs` | canonical constructors `add`, `mul`, `pow`, operators |
| `functions.rs` | function heads and their automatic evaluation |
| `diff.rs`, `expand.rs`, `subs.rs`, `evalf.rs` | core algorithms |
| `printer.rs`, `parser.rs` | `str`/`srepr` and `sympify` |
| `lambdify.rs` | numerical compilation |
| `crates/symox-py` | Python bindings, `python/symox` the Python package |

Invariants: ordinary `Expr` nodes are canonical. Only `ops`/`functions`
constructors create canonical nodes (`Expr::from_kind` is reserved for code
that guarantees the canonical form itself). The Python `evaluate=False`
constructors explicitly mark held nodes with `Expr::from_unevaluated_kind`.
`Add` stores `(term, coefficient)` pairs, `Mul`
stores `(base, exponent)` pairs, both sorted by `Expr`'s `Ord` (`Mul`
factors by base then exponent).

## License

Dual licensed under MIT or Apache-2.0, at your option.
