# Equation solvers

[README](../README.md) · [Solution checks](solver-checks.md) · [Polynomials and systems](polynomials.md) · [Real inequalities](inequalities.md) · [Symbolic sets](symbolic-sets.md)

`solve`, `solveset`, `solveset_real`, `solveset_complex`, `linsolve`,
`nonlinsolve` and `nsolve` are exported from `symox` and do not import SymPy
at runtime. Candidate checks, `unrad` and `solve_undetermined_coeffs` are in
[solution checks](solver-checks.md), `roots` and polynomial systems in
[polynomials and systems](polynomials.md), `rsolve` in
[compatibility](compatibility.md), `dsolve` and `pdsolve` in
[differential equations](differential-equations.md).

## `solve`: results and shapes

```python
from symox import *

x, y, z = symbols('x y z')
xr = Symbol('xr', real=True)

assert solve(x**2 - 4, x) == [-2, 2]
assert solve(x**2 - 4, x, dict=True) == [{x: -2}, {x: 2}]
assert solve(x**2 - 4, x, set=True) == ([x], {(-2,), (2,)})
assert solve([x + y - 3, x - y - 1], [x, y]) == {x: 2, y: 1}
assert solve([x**2 + y**2 - 1, x - y], [x, y]) == \
    [(-sqrt(2)/2, -sqrt(2)/2), (sqrt(2)/2, sqrt(2)/2)]
assert solve(x + y - 3, [x, y]) == [(3 - y, y)]
assert solve(x + y - 1) == [{x: 1 - y}]
assert solve(xr**2 + 1, xr) == []
assert str(solve(x**2 > 4, x)) == '((-oo < x) & (x < -2)) | ((2 < x) & (x < oo))'
```

As in SymPy 1.14, one equation in one unknown gives a list of values; a
system with a unique linear solution gives a dictionary; several unknowns
give a list of tuples in the order of the requested symbols, with unsolved
symbols kept as parameters; omitted unknowns are inferred from the free
symbols. `dict=True` always gives a list of dictionaries and `set=True` a
pair `(symbols, set of tuples)`. A relational input is passed to
`reduce_inequalities` and returns a boolean condition.

An empty list does not distinguish its causes. `solve` returns `[]` when the
candidates are all rejected (assumptions on the unknown, poles,
principal-branch radicals), when a system is inconsistent, and also, as in
SymPy, for an identity: `solve(x - x, x) == []` although every `x` solves
it. The native binding documents this convention
(`crates/symox-py/src/solve.rs`). When no method applies, `solve` raises
`NotImplementedError` (`could not solve ... for x`, or
`could not solve the system`); it never returns `[]` in that case. Use
`solveset` to distinguish an empty solution set from the whole domain.

High-degree polynomials give exact `CRootOf` values when no radical form is
used (`cubics=False`, `quartics=False`, or an irreducible quintic such as
`x**5 - x + 1`). Floats are rationalized before solving and the solutions
converted back: `solve(x - 0.5, x) == [0.5]`.

## `solve` flags

| Flag | Behavior |
|---|---|
| `dict`, `set` | Output shapes as above |
| `check` | `False` skips the denominator and assumption checks of candidates |
| `cubics`, `quartics` | `False` gives `CRootOf` values instead of radicals |
| `rational` | `False` solves with the floats; otherwise they are rationalized |
| `implicit` | Allows solutions that still contain the unknown (`x = f(x)`) |
| `particular`, `quick` | Sparse particular solution of an underdetermined linear system; `quick=True` requires `particular=True` (`ValueError` otherwise) |
| `exclude` | Symbols that are not inferred as unknowns |
| `manual` | Accepted; for systems the result is computed by the same native solver |
| `force` | `False` keeps the reference `_tsolve` steps without the posify retry and without forced `logcombine`/`expand_log`; the checks of candidates do not posify |
| `bivariate` | `False` skips the Lambert W step, unless `force` (default `True`) re-enables it in the posify retry, as in SymPy |
| `minimal` | `True` makes the checks of candidates stop before simplification and keep undecided candidates (`solve(sqrt(x) + y, x, minimal=True) == [y**2]`) |
| `simplify`, `numerical`, `quintics`, `incomplete`, `tsolve`, `warn` | Accepted; they do not select a different method (`warn=True` emits no warning) |

Unknown keywords are ignored, as in SymPy 1.14
(`solve(x - sqrt(y), y, positive=True) == [x**2]`). `domain` is ignored by
`solve`. Non-symbol targets (`x**2`, `f(x)`) are replaced by dummies, solved,
then restored. Matrix equations are split into their entries.

Deviations from SymPy:

- `warn=True` emits no warning for undecided candidates; SymPy warns.
- `tsolve=True` on a polynomial equation raises `NotImplementedError`;
  SymPy's direct `_tsolve` returns the incomplete `[1]` for
  `solve(x**2 - 1, x, tsolve=True)`.
- With `bivariate=False, force=False`, an equation `f(x)**g(x) = c` whose
  logarithmic form is unsolved and whose integer guesses all fail raises
  `NotImplementedError`; SymPy returns `[]`, which is wrong
  (`x**x - 2` has the root `exp(LambertW(log(2)))`). The same holds
  without flags for `(x + 1)**x - 8`.
- Solvable quintics (`x**5 - 5*x + 12`) give `CRootOf` values; SymPy's
  `roots_quintic` gives radicals (not ported).
- `simplify=False` keeps the quadratic formula over one denominator:
  `solve(a*x**2 + b*x + 1, x, simplify=False)` gives
  `(-b - sqrt(-4*a + b**2))/(2*a)`, SymPy `-b/(2*a) - sqrt(-4*a + b**2)/(2*a)`.
- Several equations that SymPy cannot solve with `force=False` get their
  Lambert W solutions (`solve(x + exp(x), x, force=False) == [-LambertW(1)]`;
  SymPy raises `NotImplementedError`).
- `particular=True` without `quick` searches a solution with the fewest
  nonzero values, including a single nonzero value that SymPy's search
  skips (see `minsolve_linear_system` in
  [solution checks](solver-checks.md)). For
  `[x + 2*y + 3*z - 6, 2*x + y - 3]`, symox returns `{x: 0, y: 3, z: 0}`
  and SymPy `{x: 3/2, y: 0, z: 3/2}`; both are valid.
- `Abs(x)` with an unknown that is neither real nor imaginary raises
  `NotImplementedError`, as in SymPy; real unknowns are split into sign
  branches.

## How `solve` proceeds

The Python layer (`python/symox/__init__.py`, `_solve_relations.py`) admits
flags, rewrites the input and shapes the output. The solving itself is the
native kernel `crates/symox-core/src/solve.rs`. For one equation in one
unknown it tries, in order:

1. clearing denominators and removing integer content;
2. polynomial equations: closed forms through degree four, factoring, and
   `CRootOf` for the remaining exact factors;
3. products, factor by factor;
4. equations in several trigonometric functions, rewritten in `tan(x/2)`;
5. powers with a symbolic exponent (`symbolic_power.rs`);
6. polynomials in a single generator such as `exp(x)` or `sin(x)`;
7. isolating the unknown and inverting the outer functions;
8. Lambert W forms (`lambert.rs`);
9. combining logarithms, expanding, clearing exponential denominators and
   combining powers;
10. rewriting trigonometric functions as exponentials.

Each candidate is checked against the denominators, the assumptions of the
unknown and, for non-polynomial equations, a numeric test (`checksol`);
infinite and NaN candidates are dropped. Radicals use the native `unrad`
port (`crates/symox-core/src/solve/unrad.rs`), and the steps of SymPy's
`_tsolve` (`solve/tsolve.rs`) are tried when the kernel fails. Python
fallbacks remain for `Piecewise`, held derivatives and integrals, functions
with a user `inverse()`, and the resultant elimination of radicals described
in [solution checks](solver-checks.md).

```python
assert solve(exp(2*x) - 3*exp(x) + 2, x) == [0, log(2)]
assert solve(x*exp(x) - 1, x) == [LambertW(1)]
assert solve(x**x - 2, x) == [exp(LambertW(log(2)))]
assert solve(sin(x) - Rational(1, 2), x) == [pi/6, 5*pi/6]
assert solve(log(x) + log(x - 1) - log(2), x) == [2]
assert solve(sqrt(x) - x + 2, x) == [4]
assert solve((x - 1)/(x**2 - 1), x) == []
assert solve((x - 1)/(x**2 - 1), x, check=False) == [1]
assert solve(x - y*exp(-x), x, implicit=True) == [y*exp(-x)]
```

`solve` returns a finite list of representatives for periodic equations, as
SymPy does; `solveset` returns the complete families.

Systems are solved by `solve_system_with_check`: parameter-only equations
are dropped unless they are proven nonzero; single-unknown equations are
isolated; linear systems use exact Gauss-Jordan elimination; a factored
equation splits into branches; other systems are solved by successive
substitution, choosing the lowest-degree pivot whose leading coefficient is
proven nonzero. Unlike SymPy, `solve` does not call `solve_poly_system` for
polynomial systems.

```python
assert solve([exp(x) - 2, y - x], [x, y]) == [(log(2), log(2))]
assert solve([x - 1, x - 2], [x, y]) == []
```

## `solveset`

`solveset(f, symbol, domain=S.Complexes)` returns a set. The result kinds
are:

| Result | Meaning |
|---|---|
| `FiniteSet` | Exactly these solutions |
| `ImageSet` or a `Union` of them | Complete periodic or logarithmic families |
| `Interval`, `Union`, `Complement` | Solution regions, poles removed |
| `Intersection(FiniteSet, Reals)` | Parametric solutions whose reality is unknown |
| `EmptySet` | Proven absence of solutions in the domain |
| the domain | Proven identity |
| `ConditionSet` | Not solved; no claim about the solution set |

```python
assert solveset(x**2 + 1, x) == FiniteSet(-I, I)
assert solveset(x**2 + 1, x, S.Reals) == S.EmptySet
assert str(solveset(exp(x) - 1, x)) == 'ImageSet(Lambda(_n, 2*_n*I*pi), Integers)'
assert solveset(x - x, x) == S.Complexes
assert solveset((x**2 - 1)/(x - 1), x) == FiniteSet(-1)
assert solveset(cos(x) + x, x, S.Reals) == \
    ConditionSet(x, Eq(x + cos(x), 0), S.Reals)
assert solveset(Abs(x) - 1, x, S.Reals) == FiniteSet(-1, 1)
assert str(solveset(Mod(x, 3) - 1, x, S.Integers)) == \
    'ImageSet(Lambda(_n, 3*_n + 1), Integers)'
assert solveset(x**2 - 1 > 0, x, S.Reals) == \
    Union(Interval.open(-oo, -1), Interval.open(1, oo))
```

The `ImageSet` parameter is a dummy `_n`, so the examples compare printed
forms. `solveset_real(f, x)` and `solveset_complex(f, x)` fix the domain; in SymPy they are imported from
`sympy.solvers.solveset`. The assumptions of `symbol` are replaced by the
domain: `solveset(xr**2 + 1, xr) == FiniteSet(-I, I)`. Several free symbols
without a given symbol raise `ValueError`. `Abs` in the complex domain
raises SymPy's `ValueError`. Inequalities are solved over real domains only;
over `Complexes` they give a `ConditionSet`.

The native core (`crates/symox-core/src/sets/solveset.rs`) proceeds as
follows: modular equations; expressions free of the symbol (domain,
`EmptySet` or `ConditionSet` according to the zero test); identities;
products as unions of their factors; polynomials in one periodic or
exponential generator, giving `ImageSet` families; mixed trigonometric
equations; complex exponentials (logarithm branches); real absolute values;
and finally the `solve` kernel, whose candidates are filtered by the domain.
Python steps then handle polynomial `Abs`, `Min`, `Max` and `Piecewise`
splits, and a native radical step with Python callbacks handles remaining
radicals. A `ConditionSet` is returned whenever no step completes. In the
complex domain, or with parameters, Lambert W equations return a
`ConditionSet` because the finite set of real branches found by `solve` is
not complete there.

Deviations from SymPy 1.14 (each is listed in
[algorithmic parity](algorithmic-parity.md) unless stated):

- Extra Lambert solutions over the reals: `solveset(x*exp(x) - 1, x,
  S.Reals)` gives `{LambertW(1)}` and `solveset(x**x - 2, x, S.Reals)`
  `{exp(LambertW(log(2)))}`; SymPy returns `ConditionSet`.
- `solveset(2**x - 4, x)` gives the complete family
  `ImageSet(Lambda(_n, 2*_n*I*pi/log(2) + 2), Integers)`; SymPy returns a
  `ConditionSet`.
- `solveset(log(x - 3) + log(x + 3), x)` gives `{sqrt(10)}`; SymPy also
  returns `-sqrt(10)`, which does not satisfy the principal-branch equation.
- `solveset(x**Rational(1, 3) + 2, x, S.Reals)` gives `EmptySet` (principal
  branch); SymPy gives `{-8}`. This gap is open.
- When the expression does not contain the symbol and its zero test is
  undecided, symox returns `ConditionSet(x, Eq(y, 0), Complexes)` for
  `solveset(y, x)`; SymPy returns `EmptySet`. Not yet listed in the parity
  document.

## Linear and nonlinear systems

```python
assert linsolve([x + y + z - 1, x + y + 2*z - 3], (x, y, z)) == \
    FiniteSet((-y - 1, y, 2))
assert linsolve([x + y - 1, x + y - 2], (x, y)) == S.EmptySet
assert linsolve(Matrix([[1, 2, 3], [4, 5, 6]]), (x, y)) == FiniteSet((-1, 2))
assert nonlinsolve([x**2 + y**2 - 1, x - y], [x, y]) == \
    FiniteSet((-sqrt(2)/2, -sqrt(2)/2), (sqrt(2)/2, sqrt(2)/2))
assert str(nonlinsolve([exp(x) - 1, y - 2], [x, y])) == \
    '{(ImageSet(Lambda(_n, 2*_n*I*pi), Integers), 2)}'
```

`linsolve` (`python/symox/_linear.py`) accepts equations and symbols, an
augmented matrix or a pair `(A, b)`; the exact reduction is the native
`_linsolve_matrix`. Free variables are kept as parameters, or named
`tau0, tau1, ...` when no symbols are given. Inconsistent systems give
`EmptySet` and nonlinear terms raise `NonlinearError`.

`nonlinsolve` (`python/symox/_polysolve.py`) first tries the native
parametric decomposition (`crates/symox-core/src/poly/parametric.rs`) for
factorable and triangular rational systems, with explicit zero/nonzero
branches. Polynomial systems then use a Groebner basis (grevlex, FGLM to
lex) and back substitution, discarding roots of denominators.
Non-polynomial systems use a port of SymPy's substitution method that keeps
`ImageSet` families. Floating coefficients are rationalized, then evaluated.
`solve_poly_system` is the native lex Groebner elimination; it returns
`None` for an inconsistent system and raises `NotImplementedError` for a
positive-dimensional one.

Deviations from SymPy:

- `nonlinsolve([x**2 - y**2], [x, y])` gives `{(-y, y), (y, y)}`; SymPy
  adds the redundant `(-sqrt(y**2), y)` and `(sqrt(y**2), y)`.
- `solve([x**5 - x + 1, y - x], [x, y])` gives the five `CRootOf` tuples;
  SymPy returns `[]`, which is wrong.

## Numerical roots

```python
assert str(nsolve(cos(x) - x, x, 1)) == '0.739085133215161'
assert str(nsolve(x**2 - 2, x, (1, 2), solver='bisect')) == '1.41421356237310'
```

`nsolve` is native (`crates/symox-py/src/solve.rs`). It supports the
`secant` (default for one equation), `newton` (default for systems) and
`bisect` solvers, with `prec`, `tol`, `maxsteps`, `verify` and `dict`.
Solver names are checked like mpmath's `findroot`: an unknown name raises
`ValueError("could not recognize solver")`, a wrong number of starting
points raises `ValueError` (`anderson` with one point: "expected interval
of 2 points, got 1"), a scalar `mdnewton` raises `KeyError('norm')`, and a
system uses the multidimensional Newton method for any valid name.
`method=` raises SymPy's `ValueError` and unknown keywords are ignored.
The valid mpmath solvers `mnewton`, `halley`, `muller`, `illinois`,
`pegasus`, `anderson`, `ridder` and `anewton` on a scalar equation raise
`NotImplementedError`; SymPy computes them.
When the residual does not reach the tolerance, both raise `ValueError`
("Could not find root within given tolerance").
The secant iteration works at `prec + 20` bits like mpmath, so
`nsolve(x**2 + 1, x, 1 + I)` prints `1.0*I`. Below that level the
rounding noise is not mpmath's: `nsolve(x**2 + 3, x, 1 + I)` prints
`1.73205080756888*I`, where SymPy prints
`6.7e-28 + 1.73205080756888*I` (both within the tolerance).
`solver='newton'` on one equation runs the damped multidimensional Newton
method with the symbolic derivative, not mpmath's scalar `Newton`:
`nsolve(x**2 + 1, x, 1 + I, solver='newton')` prints
`8.46365050617523e-23 + 1.0*I`, SymPy `1.0*I`. `verbose=True` and the
`df`, `d1f`, `d2f`, `J`, `norm` and `multidimensional` keywords raise
`NotImplementedError`.

## Remaining gaps

The SymPy 1.14 test suites of the solvers are not yet passed in full; the
counts and the open items are tracked in
[algorithmic parity](algorithmic-parity.md) (solver suites) and in `TODO.md`
(audit of the `solve`, `solveset` and `roots` options). The error messages
of unsolved equations differ from SymPy's even where both raise
`NotImplementedError`.
