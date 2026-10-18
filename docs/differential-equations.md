# Ordinary and partial differential equations

[README](../README.md) · [Solution checks](solver-checks.md)

Differential interfaces use native scalar differentiation, integration,
algebra and matrix algorithms. They do not import SymPy at runtime. The
classifiers recognize the families described here; their hint lists are not
the complete lists of SymPy's methods.

## Scalar ordinary differential equations

`dsolve(eq, func=None, hint='default', ics=None)` retains the native Rust
solver for linear constant-coefficient equations, homogeneous Cauchy–Euler
equations, first-order linear, separable and Bernoulli equations. Initial
conditions use function values and supported derivative substitutions.

`classify_ode(eq, func=None)` recognizes separable, exact, first-order linear
and Bernoulli equations, homogeneous constant-coefficient equations, finite
undetermined-coefficient forcing, second-order variation of parameters and
homogeneous Cauchy–Euler equations. `dict=True` returns the differential order,
default hint, ordered matching hints and native coefficient dictionaries.
These dictionaries are native data, not SymPy's internal solver instances.

The hints in `symox.solvers.ode.allhints` can be passed explicitly to `dsolve`.
The exact method constructs a potential function and returns an explicit
solution when isolation succeeds, or an implicit equation otherwise. The
default solver also tries this exact method and verified Lie point symmetries
when its original methods fail.
`checkodesol` supports explicit solutions, lists of solutions and implicit
solutions whose total derivatives give a provably zero residual.

```python
from symox import Function, Symbol, cos, classify_ode, dsolve, checkodesol

x = Symbol('x')
f = Function('f')(x)
equation = (x + cos(f))*f.diff(x) + f + 2*x
assert '1st_exact' in classify_ode(equation, f)
solution = dsolve(equation, f)
assert checkodesol(equation, solution, f) == (True, 0)
```

Unknown or ambiguous dependent functions require an explicit `func` argument.
Generated integration constants avoid existing free symbol names. Structural
inspection reads held derivative counts without enumerating their variables.
Differential order is limited to 256 and inspection to 10,000 graph nodes.
Power-series ODE solvers, implicit exact initial conditions and the full SymPy
classifier/hint machinery remain unimplemented.
Classification records structural matches; integration or root isolation can
still exceed the native algorithms for a matched equation.

### Arbitrary constants and Lie point symmetries

`symox.solvers.ode.constantsimp(expr, constants)` absorbs redundant arbitrary
constants while retaining parameter symbols and correlations between repeated
occurrences. It handles scalar expressions, equations and nested sequences.
Only symbols listed in `constants` can be renamed or absorbed.

`infinitesimals(eq, func, hint='default')` searches first-order point symmetry
generators. The named searches are `abaco1_simple`, `abaco1_product`,
`abaco2_similar`, `abaco2_unique_unknown`, `abaco2_unique_general`, `linear`,
`function_sum`, `bivariate` and `chi`; `hint='all'` combines their results.
Every returned generator satisfies its exact determining equation.
`checkinfsol(eq, generators, func)` independently substitutes generators and
returns `(is_valid, residual)` pairs, including nonzero residuals for invalid
generators. Polynomial searches use exact nullspaces with degree at most six;
the separated-polynomial fallback searches degrees up to three.

`dsolve(..., hint='lie_group')` constructs a first integral by integrating the
closed form `(dy - h*dx)/(eta - xi*h)` for `y' = h(x, y)`. It verifies both
partial derivatives of that integral before returning a solution. Pass both
`xi=` and `eta=` to supply a custom generator; invalid generators raise
`ValueError`. Function-value initial conditions are supported. Solutions can
be implicit and retain unevaluated integrals when native integration or root
isolation cannot produce an explicit expression.

```python
from symox.solvers.ode import infinitesimals, checkinfsol, constantsimp

equation = f.diff(x) - x**2*f
generators = infinitesimals(equation, f)
assert all(valid and residual == 0
           for valid, residual in checkinfsol(equation, generators, f))
solution = dsolve(equation, f, hint='lie_group')
assert checkodesol(equation, solution, f) == (True, 0)
```

These are bounded first-order searches, not a complete classification of all
Lie symmetries. A failed search raises `NotImplementedError` (or `ValueError`
for an explicitly selected heuristic). Higher-order symmetry solving is not
implemented. Generated constants also avoid names present in initial data.

## Linear systems

The compatible module `symox.solvers.ode.systems` provides:

- `linear_ode_to_matrix(eqs, funcs, t, order)`: exact coefficient extraction
  with the convention `A_n X^(n) = sum(A_k X^(k)) + b`. Coefficient denominators
  are preserved. Nonlinear systems and insufficient requested orders raise
  `ODENonlinearError` and `ODEOrderError`.
- `canonical_odes(eqs, funcs, t)`: isolates each function's highest derivative
  using native algebraic solving and returns a list of supported branches.
- `linodesolve(A, t, b=None, B=None, type='auto', doit=False)`: solves
  `X' = A X + b` by a matrix exponential when `A` is constant or commutes with
  its antiderivative `B`. An explicit `B` is validated by differentiation.
  Forcing integrals remain held unless `doit=True`.

Public `dsolve` accepts a square first-order linear system and a list of
dependent functions, with optional function-value initial conditions.

```python
from symox import symbols, Function, dsolve

t = symbols('t')
f, g = Function('f')(t), Function('g')(t)
solutions = dsolve([f.diff(t)-f-g, g.diff(t)-g], [f, g],
                  ics={f.subs(t, 0): 1, g.subs(t, 0): 2})
```

Matrix extraction is limited to 64 equations/functions. Matrix-exponential
solving is limited to eight functions and the native exact spectral models.
Noncommuting time-dependent systems, singular differential-algebraic systems,
general nonlinear system solving and alternate time-coordinate methods are
not implemented by this interface.

## Partial differential equations

`classify_pde`, `pdsolve` and `checkpdesol` support first-order linear equations
in two independent variables. Constant coefficients use characteristic
coordinates and an integrating factor, including inhomogeneous forcing.
`hint='1st_linear_constant_coeff_Integral'` retains a held characteristic
integral. `solvefun=Function('H')` chooses the arbitrary function; the default
name avoids collisions with function heads already in the equation. Fresh
dummy characteristic coordinates cannot capture user symbols.

```python
from symox import symbols, Function, pdsolve, checkpdesol

x, y = symbols('x y')
u = Function('u')(x, y)
equation = 2*u.diff(x) + 3*u.diff(y) + 4*u
solution = pdsolve(equation, u)
assert checkpdesol(equation, solution, u) == (True, 0)
```

The compatible method entry points in `symox.solvers.pde` are
`pde_1st_linear_constant_coeff_homogeneous`,
`pde_1st_linear_constant_coeff` and `pde_1st_linear_variable_coeff`, with
arguments `(eq, func, order, match, solvefun)`. `pdsolve` dispatches to these
methods. They accept native coefficient dictionaries or SymPy-style Wild
match dictionaries and reject coefficients inconsistent with the equation.
The inhomogeneous constant-coefficient method retains its integral; ordinary
`pdsolve` evaluates it unless the `_Integral` hint was requested. Isotropic
complex directions, including `u_x + I*u_y`, use nondegenerate coordinates.

Variable coefficients cover equations with only one derivative direction,
diagonal dilation fields, independently separable transport speeds and
supported shear fields. Other cases reduce to a native characteristic ODE,
invert its first integral and integrate the forcing along the characteristic.
Both coordinate directions are tried. This includes mixed fields such as
`(x+y)*u_x + u_y` and variable zeroth-order terms such as
`x*u_x - y*u_y + y**2*u - y**2`. Unknown forcing functions can retain held
characteristic integrals without capturing their parameters.

Classification can recognize first-order linear equations whose characteristic
ODE or coordinate inversion is beyond the native solvers; `pdsolve` then
raises `NotImplementedError`. General characteristic systems, second-order
PDE solvers and PDE boundary conditions remain unsupported.

`pde_separate`, `pde_separate_add` and `pde_separate_mul` apply an additive or
multiplicative separated ansatz and isolate the first group of variables.
Variable groups must be disjoint and cover every original argument. They
support higher derivatives, grouped variables and separable nonlinear
exponential equations. A failed separation returns `None`; it does not claim
that no separated solution exists.

```python
from symox import pde_separate_mul

X, Y = Function('X')(x), Function('Y')(y)
assert pde_separate_mul(u.diff(x, 2)-u.diff(y, 2), u, [X, Y]) is not None
```

Solution checks substitute native derivatives and simplify the residual.
Integrals with identical limits combine by linearity only when coefficients
are independent of their bound variables. Solutions and sign divisions are
local to domains where the chosen characteristic coordinates and leading
coefficients are defined; this interface does not infer global boundary data.
