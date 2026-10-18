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
powers multiply that exponent. Differentiation, absolute values and conjugation
preserve the polar exponential. A polar expression is not automatically
identified with an ordinary complex exponential.

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
native node rather than guessing a branch. General numeric approximation of
unresolved branch nodes and arbitrary precision polar evaluation are not
implemented.

`polarify(expr)` returns `(lifted_expression, reverse_substitutions)`, replacing
free symbols with polar dummy symbols. `subs=False` omits this replacement and
returns just the expression. `lift=True` explicitly lifts symbols and sums,
also returning just the expression. Integral, sum and product bound variables
are preserved.

`unpolarify(expr, subs=reverse_substitutions)` projects through sums, products
and supported unbranched functions. Polar bases of noninteger powers and
arguments of branch-sensitive functions such as log remain lifted. Integer
powers can be projected. `exponents_only=True` restricts conversion, with
zero-angle identity factors still removed. These transformations are
conservative: they do not implement every SymPy branch simplification.
