# Real inequalities

[README](../README.md) · [Solution checks](solver-checks.md)

`solve_poly_inequality(poly, rel)` returns the real interval or singleton
components for one univariate polynomial with rational coefficients. The
relations are `==`, `!=`, `<`, `<=`, `>` and `>=`. Endpoints and inclusion flags
remain exact; real roots are ordered before returning separate components.
Constant inequalities return `[S.Reals]` or `[S.EmptySet]`.

```python
from symox import Poly, symbols, solve_poly_inequality, Interval

x = symbols('x', real=True)
assert solve_poly_inequality(Poly(x*x-1, x), '<') == [Interval.open(-1, 1)]
```

`solve_rational_inequalities(groups)` takes a disjunction of groups, each group
containing `((numerator_poly, denominator_poly), relation)` conditions joined
by conjunction. The sign test uses the numerator times the denominator, then
excludes zeros of the original denominator. Cancellation therefore does not
erase a specified pole.

```python
from symox import solve_rational_inequalities, oo

conditions = [[((Poly(x*x-1, x), Poly(x-1, x)), '>=')]]
answer = solve_rational_inequalities(conditions)
assert answer == Interval(-1, 1, right_open=True).union(Interval.open(1, oo))
```

The module `symox.solvers.inequalities` also exposes
`reduce_rational_inequalities(groups, gen, relational=True)`. Conditions may be
relations, `(expression, relation)` pairs, booleans or expressions interpreted
as equations to zero. It returns a boolean condition by default, or a native
real set with `relational=False`. Supported floating polynomial coefficients
are converted to exact coefficients for solving, then endpoints are evaluated
numerically. Known-sign real parameters use the native univariate solver.

`reduce_abs_inequality(expr, rel, gen)` reduces a scalar condition against zero.
`reduce_abs_inequalities([(expr, rel), ...], gen)` joins several such conditions.
Nested absolute values, sums, products and integer powers are split into exact
real sign branches. Rational denominator poles are removed in every branch.
`solve_univariate_inequality` and `reduce_inequalities` also use these branches
when the original native solver cannot handle a nested absolute value. Domain
restrictions and the other conditions in a system are intersected exactly.

```python
from symox import Abs, And, Lt, reduce_abs_inequality

assert reduce_abs_inequality(Abs(x-5)-3, '<', x) == And(Lt(2, x), Lt(x, 8))
```

These interfaces use native polynomial arithmetic, real equation solving and
symbolic set operations. They never import SymPy at runtime. Unsupported
algebraic root representations, unordered coefficients and undecidable signs
raise an explicit error. Absolute-value expansion accepts at most 1024 sign
branches, 10,000 distinct graph nodes and 128 nested levels; rational systems
accept at most 4096 conditions. Finite coefficient fields have no ordered real
inequality domain.
