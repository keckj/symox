# Complex branches and polar conversions

[README](../README.md)

`polar_lift`, `periodic_argument`, `principal_branch`, `unbranched_argument`,
`polarify` and `unpolarify` operate on native expression nodes. Their compatible
module path is `symox.functions.elementary.complexes`.

```python
from symox import I, pi, oo, exp_polar, periodic_argument, unpolarify

z = exp_polar(5*I*pi)
assert periodic_argument(z, oo) == 5*pi
assert periodic_argument(z, 2*pi) == pi
assert unpolarify(z) == -1
```

`exp_polar` retains winding information. Its logarithm returns its exponent;
powers multiply that exponent, and `arg(exp_polar(z))` is `im(z)`. In
products, polar exponentials gather their windings by non-numeric term as in
SymPy's `Mul.flatten` (`exp_polar(I*x)*exp_polar(2*I*x)*exp_polar(I)` is
`exp_polar(I)*exp_polar(3*I*x)`). Absolute values and conjugation preserve the
polar exponential. As in SymPy, which gives `exp_polar` no derivative rule,
differentiation is held: `diff(exp_polar(x), x)` is
`Derivative(exp_polar(x), x)` and `diff(exp_polar(2*x), x)` is
`2*Subs(Derivative(exp_polar(_xi_1), _xi_1), _xi_1, 2*x)`. A polar expression
is not automatically identified with an ordinary complex exponential.

`polar_lift` puts positive real magnitudes on the zero sheet and negative real
magnitudes at angle pi. Pure imaginary numbers use angles plus or minus pi/2.
Other unresolved lifts remain explicit. Positive factors can be extracted,
and lifting a polar expression leaves it unchanged. Absolute values ignore
the lift.

`periodic_argument(z,P)` reduces a resolved unwrapped angle into `(-P/2,P/2]`.
Thus the lower endpoint maps to the upper endpoint. `P=oo`, also exposed as
`unbranched_argument`, retains all turns. Products accumulate their factor
arguments; polar powers account for both the real and imaginary parts of the
exponent. Positive real factors do not affect the angle.

`principal_branch(z,P)` combines the magnitude with a polar exponential whose
angle is reduced to that interval. `P=oo` leaves its input unchanged. When
signs, arguments or winding cannot be resolved, constructors keep a symbolic
native node rather than guessing a branch. Numeric evaluation follows SymPy:
`exp_polar` evaluates on the principal sheet only (`exp_polar(3*I*pi)` stays
held), and `exp_polar(I)`, resolved principal branches and periodic arguments
evaluate at any requested precision.

`polarify(expr)` returns `(lifted_expression, reverse_substitutions)`, replacing
free symbols with polar dummy symbols. `subs=False` omits this replacement and
returns just the expression. `lift=True` explicitly lifts symbols and sums,
also returning just the expression. The traversal is a port of SymPy's
`_polarify` (`python/symox/_polar.py`): the integration variable of an
`Integral` is not lifted, and only expression arguments are visited, so the
parameter tuples of `meijerg` and the conditions of `Piecewise` are kept.
Known difference: a numeric `Integral` or `Sum` such as
`Integral(x, (x, 1, 2))` stays `polar_lift(Integral(...))`, where SymPy's
`polar_lift` finds its argument numerically and gives
`exp_polar(0)*Abs(Integral(...))`; symox does not decide the sign of
`re`/`im` of a numeric integral or sum by numerical evaluation.

`unpolarify(expr, subs=reverse_substitutions)` projects through sums, products,
Boolean expressions, inequalities and the heads SymPy marks `unbranched`
(exponential, trigonometric, hyperbolic, gamma, error, Fresnel, Airy, beta and
Mathieu families, `Abs`, `re`, `im`). As in SymPy, an `Eq` or `Ne` is only
projected when one side is zero, and `principal_branch(z, 2*pi)` projects to
`z`. Polar bases of noninteger powers and arguments of branch-sensitive
functions such as log remain lifted. Integer powers can be projected.
`exponents_only=True` restricts conversion, with zero-angle identity factors
still removed; passing `subs` drops `exponents_only`, as in the reference.
These transformations are conservative: they do not implement every SymPy
branch simplification.
