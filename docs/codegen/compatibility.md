# Symbolic code generation compatibility

[README](../../README.md) · [Code generation](overview.md)

For SymPy-compatible scalar polynomial code generation, replace imports
consistently:

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
ordering when generators are inferred. Numeric and supported symbolic
coefficient fields support `lex`, `grlex` and `grevlex`.
`domain="ZZ"` / `"QQ"`, prime `modulus`, `polys=True` and `auto` select
the coefficient domain, return representation and ring/field division.

`printing.printer`, `str`, `precedence`, `codeprinter`, `c` and `cxx` are
ports of the SymPy 1.14 modules, so subclasses of `C99CodePrinter`,
`CXX11CodePrinter` or `CXX17CodePrinter` override `_print_<Class>`,
`type_mappings` and precedence as in SymPy, and use `doprint` and
`parenthesize`. The top-level `ccode` and `cxxcode` use these ported printers;
`ccode_function` and the [CUDA printers](cuda.md) remain native Rust
printers. `StrPrinter` is also a port; `PrettyPrinter` drives the native
layout engine. This compatibility surface targets scalar polynomial
workflows; it is not a replacement for every SymPy module.
