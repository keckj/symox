# Code generation AST

[README](../../README.md) · [Code generation](overview.md)

`symox.codegen` describes programs as trees of statement nodes that the code
printers turn into C, C++, Fortran, R or Python source. The package mirrors the
module layout and the names of SymPy 1.14's `sympy.codegen`, but most of it is a
compact reimplementation, not a copy of SymPy's sources. Only
`codegen.abstract_nodes` (`List`) is carried over from SymPy. The printers that
render the nodes (`ccode`, `cxxcode`, `fcode` and the other code printers) are
ports of SymPy's printing modules, as described in [printing](../printing.md).

## Nodes

Every node derives from `CodegenAST` and declares its constructor fields in
`_fields`, with optional `defaults` and a per-field `_construct_<field>`
coercion. Nodes are immutable and hashable and support `args`, `func`, `subs`,
`xreplace`, `replace`, `has` and `kwargs()`, so a tree can be rebuilt with
`node.func(*node.args)`. Unlike in SymPy, codegen nodes are not `Basic`
instances: `isinstance(Assignment(x, 1), Basic)` is false, and generic
expression traversals do not descend into them.

`symox.codegen.ast` provides:

- assignments: `Assignment` and the augmented forms `AddAugmentedAssignment`,
  `SubAugmentedAssignment`, `MulAugmentedAssignment`, `DivAugmentedAssignment`,
  `ModAugmentedAssignment`, also built by `aug_assign(lhs, op, rhs)`;
- blocks: `CodeBlock` with `left_hand_sides`, `topological_sort` and `cse`,
  `Scope`, `Comment`, `NoneToken`/`none`, `BreakToken`/`break_`,
  `ContinueToken`/`continue_`;
- loops and control flow: `For`, `While`, `Return`, `Print`, `Stream`
  (`stdout`, `stderr`), `QuotedString`, `String`;
- declarations: `Variable`, `Pointer`, `Element`, `Declaration` and the
  attributes `value_const`, `pointer_const`, `Attribute`, `Node`, `Token`;
- functions: `FunctionPrototype`, `FunctionDefinition`, `FunctionCall`;
- types: `Type`, `IntBaseType`, `SignedIntType`, `UnsignedIntType`,
  `FloatBaseType`, `FloatType`, `ComplexType` and the instances `untyped`,
  `real`, `integer`, `bool_`, `complex_`, `intc`, `intp`, `int8` ...
  `int64`, `uint8` ... `uint64`, `float16`, `float32`, `float64`, `float80`,
  `float128`, `float256`, `complex64`, `complex128`.

Numeric types expose `max`, `min`, `eps`, `dig` and `cast_check`:
`int8.max` is 127, `float32.eps` is `1/8388608`,
`float64.cast_check(Rational(1, 3))` gives `0.33333333333333333`, and
`int8.cast_check(300)` raises `ValueError: Value is too big: 300 > 127`.
`Type.from_expr` infers `integer`, `real` or `complex_` from numbers and
symbol assumptions, and raises `ValueError` for a symbol without them.

## A complete C function

```python
from symox import symbols, ccode, Range
from symox.codegen.ast import (AddAugmentedAssignment, Declaration, For,
    FunctionDefinition, Return, Variable, integer, real)

x, s, i, n = symbols('x s i n')
body = [Declaration(Variable(s, real, 0)),
        For(i, Range(n), [AddAugmentedAssignment(s, x**i)]),
        Return(s)]
f = FunctionDefinition(real, 'geom', [Variable(x, real), Variable(n, integer)], body)
print(ccode(f))
```

prints, as in SymPy 1.14:

```c
double geom(double x, int n){
   double s = 0;
   for (i = 0; i < n; i += 1) {
      s += pow(x, i);
   };
   return s;
}
```

`fcode(f, source_format='free', standard=95)` prints the same Fortran function
as SymPy (`real*8 function geom(x, n)`, `real*8 :: s = 0`,
`do i = 0, n - 1, 1`, `s = s + x**i`, `geom = s`). `str(node)` gives the
SymPy-style form `Declaration(Variable(s, type=real, value=0))`.

`CodeBlock.cse()` extracts common subexpressions into new assignments, and
`CodeBlock.topological_sort` orders assignments by their dependencies:

```python
from symox import sin, cos
from symox.codegen.ast import Assignment, CodeBlock

y = symbols('y')
CodeBlock(Assignment(y, sin(x)**2 + cos(x)), Assignment(s, sin(x)**2*y)).cse()
# CodeBlock(Assignment(x0, sin(x)**2), Assignment(y, x0 + cos(x)),
#           Assignment(s, x0*y))
CodeBlock.topological_sort([Assignment(s, y + 1), Assignment(y, x**2)])
# CodeBlock(Assignment(y, x**2), Assignment(s, y + 1))
```

## Language-specific nodes

`codegen.cnodes` adds `PreIncrement`, `PostIncrement`, `PreDecrement`,
`PostDecrement`, `CommaOperator`, `Label`, `goto`, `sizeof`, `alignof`,
`struct`, `union`, `void`, `restrict`, `volatile` and `static`. `ccode(PreIncrement(i))`
is `++(i)` and `ccode(sizeof('double'))` is `sizeof(double)`. Pointers render
as `double * x`, or `const double * x` with `value_const`.

`codegen.fnodes` adds `Program`, `Module`, `Subroutine`, `SubroutineCall`,
`Do`, `ImpliedDoLoop`, `ArrayConstructor`, `GoTo`, `FortranReturn`, `use`,
`use_rename`, `sum_`, `product_`, the intents `intent_in`, `intent_out`,
`intent_inout`, `array`, `dimension`, `bind_C` and the intrinsics `allocated`,
`lbound`, `ubound`, `shape`, `size`, `reshape`, `merge`, `kind`, `cmplx`,
`dsign`, `isign`, `literal_sp` and `literal_dp`. A `Subroutine` with
`array(x, n, 'in')` arguments and a `Do` loop, wrapped in
`Module('mod', ['implicit none'], [sub])` and printed with `standard=2003`,
gives the same text as SymPy.

`codegen.cxxnodes` defines `using(type, alias=None)` and `constexpr`:
`cxxcode(using(real, 'R'))` is `using R = double`, as in SymPy.
`codegen.matrix_nodes` defines `MatrixSolve`; `codegen.abstract_nodes`
defines `List`, and `codegen.pynodes` its Python subclass.
`codegen.numpy_nodes` defines `logaddexp`, `logaddexp2`, `amin`, `amax`,
`minimum` and `maximum`, and `codegen.scipy_nodes` defines `cosm1` and
`powm1`.

## C99 functions and rewriting

`codegen.cfunctions` defines `Sqrt`, `Cbrt`, `exp2`, `expm1`, `log2`, `log10`,
`log1p`, `isnan`, `isinf`, `fma` and `hypot` as symbolic functions with
derivatives, numerical evaluation and C printing (`expm1(x).diff(x)` is
`exp(x)`; `N(expm1(Rational(1, 10**10)), 20)` is
`1.0000000000500000000e-10`).

`codegen.rewriting.optimize(expr, optimizations)` applies the rules in
decreasing priority; a rule with a cost function keeps the cheaper of the
original and the rewritten expression:

```python
from symox import exp, log
from symox.codegen.rewriting import optimize, optims_c99, create_expand_pow_optimization

y = symbols('y')
optimize(log(x + 1) + 2*exp(x) - 2 + log(y)/log(2), optims_c99)
# 2*expm1(x) + log1p(x) + log2(y)
ccode(create_expand_pow_optimization(3)(x**3 + y**5))
# 'pow(y, 5) + x*x*x'
```

The module also provides `expm1_opt`, `log1p_opt`, `exp2_opt`, `log2_opt`,
`log2const_opt`, `sinc_opts`, `logsumexp_2terms_opt`, `FuncMinusOneOptim`,
`matinv_opt` and the groups `optims_numpy` and `optims_scipy`.
`codegen.approximations` provides `SumApprox` and `SeriesApprox`, and
`codegen.algorithms` provides `newtons_method` and `newtons_method_function`.

`codegen.cutils.render_as_source_file` prints C99 with `#include` lines,
`codegen.futils.render_as_module` wraps Fortran in a module, and
`codegen.pyutils.render_as_module` prints Python with `import math`.

## Printer coverage

`ccode` and `cxxcode` render the shared statement nodes (checked for
`While`, `Print`, `FunctionCall`, `Scope`, pointer declarations, `Comment` and
`break_`). `fcode` renders the Fortran nodes and most shared ones; like SymPy,
it raises for `Scope` and for pointer declarations. `rcode` prints `For`
(`for(i in seq(from=0, to=n - 1, by=1){`). `rust_code`, `jscode`, `julia_code`
and `octave_code` print `Assignment` but raise
`PrintMethodNotImplementedError` for `For`, as in SymPy. `pycode` prints the
statement nodes SymPy's `PythonCodePrinter` prints (`Assignment`,
`CodeBlock`, `FunctionDefinition`, `While`, `Return`, `Print`, `Raise`, ...)
and raises `PrintMethodNotImplementedError` for the others, such as `For`;
with `strict=False` it prints them commented as unsupported, as SymPy does.

## Differences from SymPy 1.14

- `codegen.numpy_nodes`, `codegen.scipy_nodes` and `codegen.pynodes`
  exist and `logaddexp_opt`, `cosm1_opt` and `powm1_opt` rewrite to their
  nodes, as in SymPy. `lambdify` of these nodes and of the `cfunctions`
  classes follows SymPy's printers: NumPy and SciPy fold `minimum`/`maximum`
  with `functools.reduce` (the inputs are never used as `out`), and a name a
  printer does not define (`cosm1` with NumPy, math or mpmath, `fma`,
  `Cbrt`) raises `NameError` when the generated function is called.
  `N(logaddexp2(x, y), 30)` keeps a 15-digit coefficient, and `.series()`
  of these functions raises `NotImplementedError`.
- With the default Fortran 77 standard, a declaration with an initial value
  raises `NotImplementedError` (`F77 init./parameter statem. req. multiple
  lines.`), as in SymPy, including a standalone `Declaration` with value 0,
  which SymPy prints as `real*8 s` without the initializer. A declaration
  without a value prints `real*8 s` where SymPy raises.
- `Type` and `Attribute` names are `str` (SymPy: `String`); reprs match.
- `Do(..., concurrent=True)` prints the Fortran 2008 form, and `GoTo` without
  an expression prints. About a dozen node reprs differ in detail from SymPy's
  doctests, as does `str(none)`.
- Nodes are not `Basic`: code that tests `isinstance(node, Basic)` or
  dispatches on `Basic` does not see them as SymPy code does.

Routine-level generation (`codegen`, `make_routine`) and compiled wrappers
are covered in [routines and autowrap](autowrap.md).
