# Finite differences

[README](../README.md) · [Expression calculus](compatibility.md)

`finite_diff_weights`, `apply_finite_diff`, `differentiate_finite` and
`Derivative.as_finite_difference` compute finite stencils with native scalar
arithmetic. The stencil logic is the Python module
`symox.calculus.finite_diff` over native expressions; it requires no SymPy
installation.

```python
from symox import symbols, Function, finite_diff_weights, apply_finite_diff
from symox import Derivative, differentiate_finite

x, h = symbols('x h')
f, g = Function('f'), Function('g')
assert finite_diff_weights(2, [-1, 0, 1], 0)[2][-1] == [1, -2, 1]
assert apply_finite_diff(2, [-1, 0, 1], [1, 0, 1], 0) == 2
stencil = Derivative(f(x), x).as_finite_difference(h)
assert (stencil - (f(x+h/2)-f(x-h/2))/h).expand() == 0
assert differentiate_finite(f(x), x, 2) == f(x-1)-2*f(x)+f(x+1)
```

Weights use the incremental Fornberg recurrence for distinct arbitrary grid
points. The result includes all derivative orders through the requested order,
and all prefixes of the input grid. `weights[order][prefix_length-1]` contains
a full-length coefficient list, padded with zeros for unused grid points.
Order zero computes interpolation weights. A derivative order exceeding the
available polynomial degree has zero weights.

Exact rational and symbolic grids retain exact coefficients. Unknown symbolic
differences between points are assumed nonzero. Known duplicate points and
nonfinite points are rejected. Grid point order controls the shorter prefix
stencils; the full stencil is invariant under matching point/sample permutations.

`as_finite_difference` accepts explicit grid points or a scalar step. A scalar
step produces a centered stencil with `order+1` points. `x0` changes the center,
and `wrt` selects one variable of a mixed derivative. Unselected derivative
variables remain symbolic. Undefined functions can specify variable grid steps,
evaluated at each shifted stencil center.

`differentiate_finite` defaults to forming stencils before expanding analytic
derivatives. This also supports embedded derivatives and their shifted values,
including products and powers. `evaluate=True` requests analytic expansion
first for expressions without embedded derivative nodes. Differences between
these two conventions are intentional and can affect discretized product rules.

Derivatives and grids are limited to order/size 256, and the full weight table
to 100,000 coefficients. Explicit grids need at least `order+1` points when
converting a derivative. As in SymPy, a derivative with a symbolic order
(`Derivative(f(x), (x, n))`) raises `TypeError` instead of producing a stencil.
Substitutions in nested derivative stencils currently
support one substitution variable per `Subs` node; multivariate substitutions
must be expressed with explicit nesting.
