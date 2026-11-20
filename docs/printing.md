# Printing and code generation

[README](../README.md) · [Native code generation](codegen/overview.md)

`symox.printing` exports all 41 names in the measured SymPy printing surface.
These functions operate on native expressions. The measured percentage counts
entry points, rather than complete support for every printer setting or every
SymPy object type.

## Code targets

The code printers for C (`ccode`), C++ (`cxxcode`), Fortran (`fcode`), Rust
(`rust_code`), JavaScript (`jscode`), Julia (`julia_code`), Octave
(`octave_code`), R (`rcode`), GLSL (`glsl_code`), Maple (`maple_code`) and
Mathematica (`mathematica_code`) are Python ports of the SymPy 1.14 printers
in `printing.c`, `cxx`, `fortran`, `rust`, `jscode`, `julia`, `octave`,
`rcode`, `glsl`, `maple` and `mathematica`, built on the ported
`printing.codeprinter`. Their classes, such as `C99CodePrinter`,
`CXX17CodePrinter` or `MCodePrinter`, accept `_print_<Class>` overrides as
in SymPy. `pycode` (also `symox.printing.pycode`) keeps the native printer.
Optional array and GPU backends are described in the code generation guides.

```python
from symox import symbols, sin, jscode

x = symbols('x')
assert jscode(x**2 + sin(x)) == 'Math.pow(x, 2) + Math.sin(x)'
```

`assign_to` emits an assignment (as in SymPy, `mathematica_code` has no
`assign_to`). `user_functions` supplies target names for
functions beyond each target's builtins; it can select names with predicates.
Unknown functions raise `PrintMethodNotImplementedError` unless
`allow_unknown_functions=True`; Maple and Mathematica print them by name
(`f[x]` in Mathematica). `strict=False` emits partial code with a
comment listing the unsupported functions. The target runtime must provide
any user-mapped functions and optional special-function libraries. The C, C++,
Rust and R printers suffix reserved identifiers. `human=False` returns the constants,
unsupported set and code tuple instead of a string. Piecewise expressions need a final True
branch (`ValueError` otherwise), except in R, which emits `NA`, and Python,
which emits `None`. As in SymPy, `mathematica_code` has no Piecewise printer
and raises `PrintMethodNotImplementedError` for every Piecewise.

With `assign_to` and ranged `Idx` indices, the C, Fortran, Rust, JavaScript,
R, GLSL, Julia and Octave printers generate the reference loops, including
the summation loops of a contraction. Matrix `assign_to` prints one
assignment per entry. Julia, Octave, GLSL, Maple and Mathematica also print
inline matrices; JavaScript and R need a matrix `assign_to`. These printers
do not cover the full expression vocabularies of all target systems. Code
fragments require a suitable surrounding program and runtime. Scalar complex
arithmetic in JavaScript and GLSL is unsupported. Known differences from
the reference printers are listed in
[the parity contract](algorithmic-parity.md).

## MathML, expression graphs and symbolic programs

`mathml(expr, printer='content')` returns XML directly from the native tree.
Content output covers arithmetic, elementary and named functions, constants,
relations, Piecewise, explicit matrices, integrals, sums, products and
formal derivatives. `printer='presentation'` provides presentation markup for
scalar expressions and matrices, with generic function notation for less
common heads. Symbols and text are XML-escaped. `symbol_names` customizes
symbol labels; the full SymPy MathML settings vocabulary is not implemented.

`dotprint` emits Graphviz DOT, with node styles, label callbacks, repeated-node
control and a depth limit. Graphviz is needed only to render the emitted text.
`print_tree` writes an indented tree, optionally including native assumption
facts.

`python` returns a program importing Symox and assigning a reconstructed
symbolic object to `e`. `srepr` provides its constructor representation.
`pycode` instead emits a numerical Python expression. `sstr` and `sstrrepr`
provide ordinary text and quoted string output.

## Layout and console output

`TableForm` formats rectangular or ragged rows, headings, column alignment,
padding and zero suppression, and can return a native matrix or LaTeX table.
Column formats customize LaTeX entries. Input containers are copied.
`multiline_latex` splits additive expressions across `align*`, `eqnarray` or
`IEEEeqnarray` lines, with configurable term counts and continuation dots.

`pprint`, `pretty_print` and `pretty` use the native layout engine.
`pprint_use_unicode` gets or changes the default Unicode preference, returning
the previous setting. `pprint_try_use_unicode` checks the output encoding.
`pager_print` sends pretty output to Python's pager. The `print_*` entry points
write the corresponding formatted result to standard output.

`preview` produces PNG, SVG, PDF, PS or DVI through a TeX installation and the
appropriate converters. It supports files, BytesIO buffers, a custom preamble,
packages, converter options and viewer callbacks. The default Euler font needs
its TeX package; `euler=False` selects the standard font. For headless output,
use `viewer='file'` with `filename`, or `viewer='BytesIO'` with `outputbuffer`.
`print_gtk` emits presentation MathML to GtkMathView's optional `mathmlviewer`
executable; `start_viewer=False` validates the document without opening it.
Missing external renderers produce an explicit dependency error.

## SMT-LIB

`smtlib_code` emits typed variable declarations, arithmetic, comparisons,
Boolean expressions and conditional `ite` expressions. Integer assumptions
infer Int; other numeric symbols default to Real. `symbol_table` selects
explicit types, and `known_types`, `known_constants` and `known_functions`
customize mappings. Undefined functions require an argument/result type tuple
in `symbol_table`, or an explicit function mapping and caller-provided
prologue. `auto_declare`, `auto_assert`, `prefix_expressions` and
`suffix_expressions` control the surrounding script.

The default transcendental functions and `pow` follow dReal-style extensions;
other solvers may require mappings. Ordinary linear integer/real constraints
are tested with Z3. General quantifiers, complete type inference across arbitrary
expressions, all SymPy assumption predicates and full solver portability are
outside this implementation. `precision` is accepted for compatibility;
float literals retain their native textual precision.
