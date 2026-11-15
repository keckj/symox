# Routine generation and autowrap

[README](../../README.md) · [Code generation](overview.md)

`symox.utilities.codegen` writes complete source files (functions with typed
arguments, headers and modules) from expressions, and
`symox.utilities.autowrap` compiles such routines into callables. For
expression-to-callable translation without a compiler, use `lambdify`, described
in [compatibility](../compatibility.md) and the
[API implementation notes](../api-implementation.md). For statement-level
program trees, see the [code generation AST](ast.md).

## Routines and source files

`utilities.codegen` is SymPy 1.14's module with its imports redirected to symox
and two adaptations: matrix outputs are held in a small tuple container, and
argument type detection reads `is_integer` and `is_real` through `getattr`.
It provides `codegen`, `make_routine`, `get_code_generator`, `Routine`, the
argument classes `InputArgument`, `OutputArgument`, `InOutArgument` and
`Result`, `DataType`, `CodeGenError`, `CodeGenArgumentListError` and the
generators `CCodeGen` (`C89CodeGen`, `C99CodeGen`), `FCodeGen`,
`JuliaCodeGen`, `OctaveCodeGen` and `RustCodeGen`.

```python
from symox import symbols, Eq
from symox.utilities.codegen import codegen, make_routine

x, y, z = symbols('x y z')
[(c_name, c_code), (h_name, h_code)] = codegen(
    ('f', x + y*z), 'C99', 'test', header=False, empty=False)
print(c_code)
```

prints, as SymPy does:

```c
#include "test.h"
#include <math.h>
double f(double x, double y, double z) {
   double f_result;
   f_result = x + y*z;
   return f_result;
}
```

and `h_code` declares the prototype inside a `PROJECT__TEST__H` include guard.
An equation names an output argument: `codegen(('g', Eq(z, x*y)), 'F95', 'g')`
generates `subroutine g(x, y, z)` with `REAL*8, intent(out) :: z`. A routine
can return several values:

```python
r = make_routine('fcn', [x*y, Eq(z, x + y)])
[type(a).__name__ for a in r.arguments]   # ['InputArgument', 'InputArgument', 'OutputArgument']
[type(v).__name__ for v in r.results]     # ['Result']
```

Rust, Julia and Octave generators produce the same text as SymPy for these
scalar cases, for example
`fn f(x: f64, y: f64) -> f64 { let out1 = y*x.sin(); out1 }` for
`codegen(('f', y*sin(x)), 'Rust', 'f')`. Without `header=False` the generated
files start with SymPy's banner, naming symox and its release
(`Code generated with Symox 0.1.1` instead of `... SymPy 1.14.0`).

Matrix outputs and inputs (`MatrixSymbol`, including `InOutArgument`s) and
`Indexed`/`Idx` loops generate the same code as SymPy:
`codegen(('rot', Eq(MatrixSymbol('A', 2, 1), Matrix([x*y, x - y]))), 'C99',
'rot')` gives `void rot(double x, double y, double *A)`. C row-major strides
of 2-D arrays print as `A[i*n + j]` where SymPy prints `A[n*i + j]`.

## Compiled wrappers

`autowrap(expr, language=None, backend='f2py', tempdir=None, args=None,
flags=None, verbose=False, helpers=None)` generates a routine named
`autofunc`, compiles it and returns the compiled function. Arguments are
sorted by name unless `args` gives their order. Backends:

- `'f2py'` (default) writes Fortran 95 and runs
  `python -m numpy.f2py -c -m <module> <file>.f90`. With NumPy 1.26 or later
  and Meson installed, and no explicit `flags`, it adds `--backend meson`.
  It needs NumPy, a Fortran compiler and, on that path, Meson and Ninja.
- `'cython'` writes C, a `.pyx` wrapper and a `setup.py`, and builds with
  Cython and setuptools. Only scalar `double` arguments are supported; array
  arguments raise `NotImplementedError`. SymPy's Cython options
  (`include_dirs`, `libraries`, `extra_compile_args`, ...) are accepted but
  not used. Without Cython it raises `ImportError`.
- `'dummy'` compiles nothing: the result has `args` and `returns` attributes
  and returns the expression text, as in SymPy.

Any other backend raises `ValueError`. Generated modules are named
`symox_wrapped_<n>`. Without `tempdir`, sources and binaries go to a new
`symox-autowrap-*` directory under the system temporary directory, which is
removed once the module is loaded, as SymPy does; a given `tempdir` is kept.

```python
from symox import Matrix
from symox.utilities.autowrap import autowrap, binary_function, ufuncify

f = autowrap(x**2 + y)               # f2py: compiles Fortran
f(2.0, 0.5)                          # 4.5
g = autowrap(x**2 + y, args=(y, x))  # explicit argument order
g(2.0, 3.0)                          # 11.0
autowrap(Matrix([x*y, x - y]))(2.0, 3.0)   # array([[ 6.], [-1.]])
```

`binary_function('h', expr, **kwargs)` returns
`implemented_function('h', autowrap(expr, **kwargs))`: `lambdify` calls the
compiled code, but `h(2, 3).evalf()` stays `h(2, 3)` (SymPy evaluates it
through the attached implementation; see the differences below).

## NumPy ufuncs

`ufuncify(args, expr, backend='numpy')` writes a C source file for a genuine
NumPy ufunc and compiles it with `cc -shared -fPIC -O3 -std=c99` against
NumPy's and Python's headers. All inputs and outputs are `double`; a list of
expressions gives several outputs, and NumPy's limit of 32 inputs plus
outputs applies. `helpers` is a list of `(name, expr, args)` triples compiled
as static C functions that `expr` can call. Any other backend returns
`numpy.vectorize` over the scalar `autowrap` result.

```python
import numpy as np
from symox import sin, exp, Function

u = ufuncify((x, y), x**2 + y)
u(np.array([1.0, 2.0, 3.0]), 0.5)        # array([1.5, 4.5, 9.5])
(u.nin, u.nout)                          # (2, 1)
ufuncify(x, sin(x)*exp(-x))(np.linspace(0, 1, 3))
# array([0.        , 0.29078629, 0.30955988])
h = Function('h')
ufuncify(x, h(x) + 1, helpers=[('h', x**2, [x])])(np.array([1.0, 2.0]))
# array([2., 5.])
```

## What was exercised

On Linux with gcc, gfortran, NumPy 2.2, Meson and Ninja, the examples above
were run: the `numpy` ufunc path, the `f2py` path (scalar, reordered and
matrix outputs) and the `dummy` backend. Each compilation took between
under a second (ufunc) and about two seconds (f2py). In the same environment
SymPy 1.14's f2py backend fails with `CodeWrapError`, because it calls
`numpy.f2py` without the Meson backend. The Cython backend was not exercised
(Cython was not installed). Windows, macOS, MSVC and non-GCC compilers were
not tested; the ufunc path calls `cc` directly and has no compiler selection.

## Differences from SymPy 1.14

- `autowrap.py` is a short reimplementation, not SymPy's module: wrapper
  classes do not expose SymPy's setup-file templates, the Cython backend
  ignores its compiler options.
- `implemented_function` and `binary_function` results evaluate under `evalf`
  and `N` as in SymPy (`implemented_function('h', lambda a, b: a*b + 1)(2,
  3).evalf()` is `7.00000000000000`). Unlike SymPy, a plain `Function('h')`
  of the same name also evaluates once `h` has been implemented.
