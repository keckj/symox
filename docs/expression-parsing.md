# Expression parsing

[README](../README.md)

`parse_expr` is available at the package root and through
`symox.parsing.sympy_parser`. Token transformations construct native Symox
objects, without importing SymPy during parsing.

```python
from symox import symbols, Rational
from symox.parsing.sympy_parser import parse_expr

x, y, z = symbols('x y z')
assert parse_expr('2x + 3xyz', transformations='all') == 2*x + 3*x*y*z
assert parse_expr('0.[3]') == Rational(1, 3)
```

`parse_expr` evaluates transformed Python code and accepts trusted input.
The token transformations do not provide a sandbox.

## Transformations and namespaces

The standard transformations are `lambda_notation`, `auto_symbol`,
`repeated_decimals`, `auto_number` and `factorial_notation`. They support native
symbols and undefined functions, exact integer literals, decimal literals,
repeating decimals, symbolic lambdas, factorials and double factorials.

Optional transformations add:

- `convert_xor`: interpret `^` as exponentiation.
- `implicit_multiplication`: accept expressions such as `2x` or `2(x+1)`.
- `implicit_application`: accept function calls such as `sin x`.
- `implicit_multiplication_application`: combine splitting, multiplication,
  application and function powers, including `sin^2 x` and `3xyz`.
- `split_symbols` and `split_symbols_custom`: choose which names split into
  individual factors. Greek names and names with underscores remain intact by
  default.
- `convert_equals_signs`: turn eligible `=` tokens into symbolic equalities.
- `function_exponentiation`: accept notation such as `sin**2(x)`.
- `rationalize`: turn decimal literals into exact rationals.

Pass a tuple of transformations, `transformations='implicit'` or
`transformations='all'`. `T` supports indexing and slicing; `T[:5]` selects the
standard transformations. Transformation order affects interpretation.

`local_dict` supplies symbols, function heads or numerical callables, taking
precedence over the default namespace. `global_dict` replaces that namespace.
By default it contains public Symox exports and Python builtins. Names such as
`E`, `I`, `Q`, `sin` and `beta` therefore retain their existing meanings.
The parser's mathematical capabilities follow the native constructors and
functions available in the chosen namespace.

## Held expressions and lower-level helpers

`evaluate=False` constructs held arithmetic and elementary function nodes and
preserves input argument order. For example, `parse_expr('1+x', evaluate=False)`
and `parse_expr('x+1', evaluate=False)` have different argument tuples.
`sin(0)`, `sqrt(4)` and `Eq(1,1)` remain explicit under this option. Comparisons
can be chained in this mode through a held conjunction.

`stringify_expr` returns transformed Python source. `eval_expr` evaluates that
source with explicit dictionaries. `evaluateFalse` returns the transformed
Python AST, and `EvaluateFalseTransformer` exposes the AST visitor.
`AppliedFunction` and `ParenthesisGroup` expose the token groups used by the
implicit transformations. These helpers have the same trusted-input scope as
`parse_expr`.

Held construction controls parsing; subsequent native transformations,
substitutions and conversion to another engine can evaluate nodes. Not every
SymPy constructor or constructor option is implemented by the native engine.

## Numerical C and Fortran programs

`symox.parsing.sym_expr.SymPyExpression(source, mode='c')` and `mode='f'`
parse the supported numerical subset of C or Fortran into code generation nodes.
The object can convert this AST to Python, C or Fortran. Declarations,
assignments, functions and recognized loops are supported. This is a bounded
parser for the [documented subset](api-implementation.md#intégration-polynomiale-et-analyse-de-programmes),
with explicit errors for unsupported constructs; it does not implement a full
compiler frontend or all source-language type rules.
