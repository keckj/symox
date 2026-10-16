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
`solve`, `series`/`limit`, matrices, `simplify`, integration, LaTeX
printing, code generation (C/CUDA), arbitrary precision `evalf`,
assumptions.

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

Invariants: every `Expr` is canonical. Only `ops`/`functions` constructors
create nodes (`Expr::from_kind` is reserved for code that guarantees the
canonical form itself). `Add` stores `(term, coefficient)` pairs, `Mul`
stores `(base, exponent)` pairs, both sorted by `Expr`'s `Ord` (`Mul`
factors by base then exponent).

## License

Dual licensed under MIT or Apache-2.0, at your option.
