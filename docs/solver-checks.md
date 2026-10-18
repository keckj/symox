# Solution checks and linear isolation

[README](../README.md) · [Polynomials and systems](polynomials.md)

`checksol`, `solve_linear` and `solve_linear_system_LU` are exported from
`symox` and `symox.solvers.solvers`. They run on native scalar expressions
and matrices, without importing SymPy.

## Checking solutions

```python
from symox import symbols, sqrt, checksol, Eq

x, y = symbols('x y')
assert checksol(x*x-2, x, sqrt(2)) is True
assert checksol(1/(x-1), x, 1) is False
assert checksol([Eq(x+y, 3), Eq(x-y, 1)], {x: 2, y: 1}) is True
```

The function accepts an expression equal to zero, an equality, a collection of
equations or matrix entries. Values can be supplied as `(symbol, value)` or a
substitution mapping. Mappings are copied. Results are `True`, `False` or
`None` when available native simplification cannot decide. A false equation
short-circuits a collection; an undecided equation makes the collection
undecided unless another equation is false.

Exact rational residuals are checked exactly. The default `numerical=True`
also accepts numerical residuals below an absolute tolerance of `1e-9`, evaluated
with 30 decimal digits. This accepts approximate roots and is not an exact
identity proof. `numerical=False` disables this approximate check. NaN and
infinite substitutions or residuals are rejected.

`minimal=True` limits transformations, `simplify=False` omits optional global
simplification, and `warn=True` emits a warning for undecided results. `force`
is accepted as a compatibility hint; it does not introduce positivity assumptions
or override principal branches. Unsupported options raise `TypeError`.

## Isolating one linear variable

```python
from symox import solve_linear

assert solve_linear(2*x+y-1) == (x, (1-y)/2)
assert solve_linear(x+y, symbols=[y]) == (y, -x)
assert solve_linear(1/(1/x-2)) == (0, 0)
```

`solve_linear(lhs, rhs=0, symbols=(), exclude=())` analyzes a rational numerator
and tries ordinary scalar symbols in canonical order. Selected symbols and
excluded symbols can be supplied explicitly. Relevant held derivatives are
evaluated before isolation. Original negative-power bases are retained for
pole checks, including poles hidden in nested rational expressions.

The return conventions are:

| Result | Meaning |
|---|---|
| `(symbol, value)` | A linear variable was isolated without a known denominator violation |
| `(0, 1)` | The numerator is independent of the selected variables, including an identically zero numerator |
| `(0, 0)` | A single-symbol numerator has only a candidate excluded by a pole |
| `(numerator, denominator)` | No supported linear isolation was found |

Common denominators can differ structurally from SymPy while representing the
same rational expression. Symbolic coefficient denominators are treated
generically when their nonzero status cannot be decided. This is scalar linear
isolation rather than a general nonlinear solver.

## Square augmented systems

`solve_linear_system_LU(augmented, symbols)` accepts a square coefficient
matrix with one appended right hand side column and returns a symbol/value
dictionary. It uses native `LUsolve`, preserves the input matrix and rejects
invalid shapes, missing or duplicate symbols and noninvertible systems.

Structural rational isolation is limited to 256 nested arithmetic levels and
the pole/derivative traversal to 10,000 distinct nodes.

## Denominators, matrix shortcuts and sparse solutions

`symox.solvers.solvers.denoms(expr, *symbols)` recursively collects denominator
factors, optionally restricted to factors containing the selected symbols.
It accepts equations and expression collections. `recast_to_symbols(eqs,
symbols)` replaces nonsymbol unknowns with fresh dummy symbols and returns
the changed equations, unknowns in the original order and a restoration map.
Overlapping replacements are simultaneous.

`det_perm`, `det_minor` and `det_quick` compute exact determinants;
`inv_quick` computes exact inverses. Small symbolic matrices use Leibniz
expansion, memoized Laplace minors or cofactors; larger matrices use native
matrix kernels. Explicit `det_quick(..., method=...)` selects the native
determinant method. All require square inputs and preserve them. `det_quick`
and `inv_quick` are also exported at the package root. Permutation expansion
is limited to order eight, minors to twelve and cofactor inversion to six.

`minsolve_linear_system(augmented, *symbols, quick=False)` finds a particular
solution with minimal nonzero support. Exact search enumerates smaller supports
and checks consistency with native linear solving. `quick=True` uses greedy
support reduction and does not promise the global minimum. Nontrivial
homogeneous solutions follow the linear solver's free-variable convention.
Inconsistent systems return `None`; fully unconstrained systems return `{}`.
The input is preserved and coefficient parameters are retained. Exact search
is limited to 20 unknowns and greedy solving to 64.

## Removing radicals

`symox.solvers.solvers.unrad(eq, *symbols, take=None)` eliminates algebraic
radicals using native polynomial resultants. Nested radicals are replaced
with fresh coordinates and eliminated from the outside inward. Multiple
rational powers of one base can use a shared coordinate instead, returning
its defining relation alongside the transformed equation.

```python
from symox.solvers.solvers import unrad
from symox import root

polynomial, cov = unrad(sqrt(x)+root(x, 3)-2)
p, relation = cov
assert polynomial == p**3+p**2-2
assert relation == p**6-x
```

The result is `(equation, covariance)`; `None` means there were no selected
radicals. The transformed equation can have extraneous roots: candidates
must be checked in the original equation and against its denominator and
branch conditions. Elimination is limited to eight radical generators,
individual degrees up to 64 and combined resultant degree up to 256.
Shared-coordinate degree is limited to 64. Symbolic exponents and
transcendental functions of radical coordinates are outside this polynomial
elimination interface; unsupported transformations raise `NotImplementedError`.

## Undetermined coefficients

`solve_undetermined_coeffs(equation, coefficients, *variables, **flags)` matches
independent polynomial or functional generators, then solves their coefficient
equations with the native solver. Coefficients can enter nonlinearly, and bases
can include exponentials or trigonometric functions. Explicit variables can be
given individually or in a list. Without them, the common variables in the
coefficient-dependent and coefficient-independent parts are used first, keeping
other symbols as parameters when possible.

```python
from symox import solve_undetermined_coeffs

a, b = symbols('a b')
assert solve_undetermined_coeffs(a*x+b-2*x-3, [a, b], x) == {a: 2, b: 3}
assert solve_undetermined_coeffs(a*a*x+b-x, [a, b], x) == [{a: -1, b: 0}, {a: 1, b: 0}]
```

A unique solution returns a dictionary; multiple solutions return a list of
dictionaries. `dict=True` always returns a list. `set=True` returns solved
symbols and a set of solution tuples. Empty coefficient matches return `[]`.
Generators containing an unknown coefficient, such as `a**x` when solving for
`a`, return `None`. Arbitrary identities between dependent functional generators
are not resolved by this coefficient-matching method. Solver options are passed
to the native system solver; candidate denominator checks are omitted only
when explicitly requested with `check=False`.

The underlying scalar `Expr.as_coefficients_dict(*generators)` also supports
generator-dependent collection directly. Missing keys have coefficient zero;
held sums retain their public terms while separating numeric coefficients.

## Functional composition

`decompogen(f, x)` returns a chain with the outer function first.
`symox.solvers.decompogen.compogen(chain, x)` reconstructs its composition.

```python
from symox import decompogen, sin, cos
from symox.solvers.decompogen import compogen

chain = decompogen(sin(cos(x*x+1)), x)
assert chain == [sin(x), cos(x), x*x+1]
assert compogen(chain, x) == sin(cos(x*x+1))
```

Supported decompositions include unary function nesting, powers with one varying
argument, a common polynomial generator, univariate polynomial composition and
`Min`/`Max` arguments with a common inner chain. Variable exponents retain their
base: `2**(x*x+1)` decomposes into `[2**x, x*x+1]`. When no supported chain is
found, the result is `[f]`. Expressions must be arithmetic scalars and the
variable must be a symbol. Chains and recursion are limited to 256 functions or
levels; expression traversal is limited to 10,000 distinct nodes.

## Differential order

`ode_order(equation, f(x))` returns the highest total derivative order involving
that dependent function. It accepts expressions and equations, distinguishes
other dependent functions, and adds outer derivative orders for held composite
derivatives. A power such as `Derivative(f(x), x)**2` still has order one.
Mixed derivatives use the sum of their variable counts. This is structural
analysis of the held equation; it does not evaluate the derivatives.

The implementation visits shared subexpressions once, with a 100,000-node
limit, and reads derivative counts without enumerating their variables.

`homogeneous_order(expression, *coordinates)` proves a degree under simultaneous
positive scaling of the supplied coordinates. It handles arithmetic combinations,
constant powers, degree-zero function arguments and declared applied functions
whose arguments are also supplied. For example,
`homogeneous_order(x*x/sqrt(x*x+y*y), x, y)` is one.

Absolute values and complex projections follow positive scaling directly:
`Abs(x)` has degree one and `Abs(x**I)` has degree zero. Symbolic exponents can
cancel, as in `x**(n+1)*y**(-n)`, whose degree is one. These proofs include cases
where SymPy 1.14 returns `None`. A fractional power of a base with a complex
scaling degree remains unresolved to preserve its principal branch. Held
derivatives, order terms and other unproved degrees return `None`. The input
graph is limited to 10,000 distinct scalar nodes.

The native ODE solver and candidate checker limit differential order to 256
before allocating coefficient or derivative tables. Unsupported solver requests
raise `NotImplementedError`; the current candidate checker returns `False`
when it cannot validate an excessive order. Derivative placeholder conversion
retains repeated subgraphs rather than unfolding them.

User parameters and the independent variable are isolated from temporary solver
symbols. Generated integration constants skip names present in the equation or
initial conditions. Thus `f' - C1*f = 0` with `f(0)=3` produces
`f(x)=3*exp(C1*x)`, retaining the parameter `C1`. Initial conditions must refer
to the requested function; unrelated functions and derivative orders above 256
are rejected.
