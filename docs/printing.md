# Printing and code generation

[README](../README.md) · [Native code generation](codegen/overview.md)

`symox.printing` exports all 41 names in the measured SymPy printing surface.
These functions operate on native expressions. The measured percentage counts
entry points, rather than complete support for every printer setting or every
SymPy object type.

## Code targets

The existing native printers are exposed as `ccode`, `cxxcode`, `fcode`,
`pycode` and `rust_code`. Compatible modules include `printing.c`, `cxx`,
`fortran`, `pycode` and `rust`. Their scalar capabilities and optional array
backends are described in the code generation guides.

Additional scalar printers support JavaScript (`jscode`), Julia (`julia_code`),
Octave (`octave_code`), R (`rcode`), Maple (`maple_code`), Mathematica
(`mathematica_code`) and GLSL (`glsl_code`). They support scalar arithmetic,
real constants, elementary mathematical functions, relations, Boolean
connectives and Piecewise expressions with a final True branch.

```python
from symox import symbols, sin, jscode

x = symbols('x')
assert jscode(x**2 + sin(x)) == 'Math.pow(x, 2) + Math.sin(x)'
```

`assign_to` emits an assignment. `user_functions` supplies target names for
functions beyond each target's builtins; it can select names with predicates.
Unknown functions raise an error unless explicitly permitted. The target
runtime must provide any user-mapped functions and optional special-function
libraries. Reserved identifiers can be suffixed or rejected. Ambiguous symbol
names are rejected. `human=False` returns the constants, unsupported set and
code tuple for these additional scalar printers.

Julia and Octave use elementwise arithmetic operators. JavaScript, Julia,
Octave, R and Mathematica also accept concrete matrices, emitting their target
literal syntax. These printers do not generate indexed contraction loops,
perform matrix shape inference or cover the full expression vocabularies of
all target systems. Code fragments require a suitable surrounding program and
runtime. Scalar complex arithmetic in JavaScript and GLSL is unsupported.

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
