# Differential geometry

[README](../README.md) · [Coordinate vectors](coordinate-vector.md) ·
[Differential equations](differential-equations.md)

`symox.diffgeom` provides the 24 exports from the corresponding SymPy 1.14
package. This is a name-presence measure. The actual operations and their
limits are described here; complete behavioral compatibility is not implied.
The interfaces follow the [SymPy differential-geometry documentation](https://docs.sympy.org/latest/modules/diffgeom.html).
Geometric structure lives in Python objects; scalar derivatives, simplification,
equation solving, matrices and component arrays use the native Symox engine.
There is no runtime dependency on SymPy.

## Manifolds, charts and scalar fields

`Manifold(name, dim)` and `Patch(name, manifold)` represent a named manifold of
positive integer dimension and a local patch. `CoordSystem(name, patch,
symbols=None, relations=None)` supplies local coordinates. Explicit relations
map `(source_name, target_name)` to `(input_symbols, output_expressions)`.

`transform`, `transformation`, `jacobian` and `jacobian_determinant` evaluate
maps, compose declared paths and differentiate them. A missing inverse is
constructed only when the native solver produces a unique complete solution.
Ambiguous or unsupported inverses require an explicit relation with the chosen
branch. The returned transformation is a callable matrix map with `signature`,
`variables` and `expr` attributes. It prints as SymPy's
`Lambda((x, y), Matrix(...))`, but it is not a native scalar `Lambda`
containing a matrix. The legacy `connect_to` method adds explicit relations; gap filling
is replaced by composition on demand.

`Point(coord_sys, coords)` stores coordinate values and converts them through
the same transformation graph. `base_scalars`, `base_vectors` and
`base_oneforms` expose local coordinate functions, derivations and differentials.
`symox.diffgeom.rn` includes rectangular/polar `R2` charts and
rectangular/cylindrical/spherical `R3` charts, together with their declared maps.

```python
from symox.diffgeom.rn import R2_r, R2_p

x, y = R2_r.base_scalars()
e_x, e_y = R2_r.base_vectors()
p = R2_r.point((3, 4))
assert x(p) == 3
assert R2_p.base_scalar(0)(p) == 5
assert e_x(x*y) == y
assert (x*y)(p) == 12
```

`CoordinateSymbol` and `BaseScalarField` have separate identities. The former
represents formal transformation coordinates; the latter is a function on
points. Ordinary scalar symbols are parameters held constant by vector fields.
Coordinate atoms and scalar-field atoms are separate native symbol classes,
one per chart, with Python wrappers carrying `coord_sys`, `index`, `patch`,
`manifold` and `dim`. Both print by their coordinate name, as in SymPy
(`x**2`, LaTeX `\mathbf{x}` for a field), so a coordinate and the field of
the same name print alike while staying distinct. Same-named coordinates of
different charts also stay distinct. `CoordinateSymbol.rewrite(chart)`
expresses a coordinate in another chart (`x.rewrite(R2_p)` is
`rho*cos(theta)`). The `free_symbols` of a scalar expression include these
atoms, unlike SymPy's scalar-field metadata. `expr.atoms(BaseScalarField)`
recovers them; scalar conversion to SymPy produces plain symbols rather than
reconstructing geometric objects. `expr(point)` and `expr.rcall(point)`
evaluate scalar expressions at points after importing this package.

These classes represent local coordinate data. They do not certify atlas
coverage, domains of maps, global inverses, branch choices or smoothness across
chart singularities. Relation graphs use chart names on a manifold; names must
identify the intended charts without conflicting relations.

## Tensors and exterior calculus

Vector fields, covectors and mixed tensors store finite coordinate-component
maps with explicit slot variance. `TensorProduct` builds products; tensor
multiplication with `*` is reserved for scalar coefficients. Calling a tensor
contracts its slots or lets contravariant slots act on scalar functions.
Omitted arguments and `None` leave slots uncontracted. `covariant_order` and
`contravariant_order` report the slot counts. Tensor addition requires the same
slot order and variance. These are structural Python tensor objects, not a
complete subclass hierarchy of native scalar expressions.

`WedgeProduct` performs normalized antisymmetrization, including nested products.
`Differential` implements the exterior derivative of scalar fields and
alternating covariant forms. Higher nonalternating covariant tensors are
rejected. Constant differentials and identically zero tensors normalize to a
callable zero that compares equal to the native scalar zero. `Commutator`
evaluates the Lie bracket of coordinate vector fields, and `LieDerivative`
evaluates scalar, vector, covariant and mixed tensor derivatives. Evaluation is
eager, where SymPy keeps lazy objects (the values are equal); these classes are
not a general algebra of held differential operators.

Tensors print as sums of `TensorProduct` terms of base fields and
differentials, for example `x*dy + y*dx` or
`TensorProduct(dx, dy) - TensorProduct(dy, dx)` for an evaluated wedge
product. SymPy's LaTeX `\otimes`/`\wedge` forms and the unevaluated
`Differential` of a scalar expression are not reproduced.

```python
from symox.diffgeom import Differential, TensorProduct, WedgeProduct, Commutator

dx, dy = R2_r.base_oneforms()
assert Differential(Differential(x*y)) == 0
assert WedgeProduct(dx, dy)(e_x, e_y) == 1
assert WedgeProduct(dx, dy)(e_y, e_x) == -1
assert TensorProduct(dx, dy)(e_x) == dy
assert Commutator(e_x, x*e_y) == e_y
```

`vectors_in_basis` converts contravariant bases using chart Jacobians, retaining
scalar-field coefficients in their original charts. Differential operators
convert coefficients into their working chart as needed. For a tensor written
in several coordinate bases, first choose an unambiguous working basis;
automatic resolution of arbitrary mixed-basis operator expressions is limited.
Coefficient substitution, differentiation, expansion and simplification are
available as methods; public `symox.diff/simplify/expand` dispatch to these
objects after the package is imported.

## Connections, curvature and integral curves

`twoform_to_matrix` converts a rank-two covariant tensor. It follows SymPy's
actual convention: the entry in row `j`, column `i` evaluates `(e_i, e_j)`.
For a symmetric metric this gives its usual coefficient matrix.

`metric_to_Christoffel_1st`, `metric_to_Christoffel_2nd`,
`metric_to_Riemann_components` and `metric_to_Ricci_components` compute the
Levi-Civita connection and its curvature in that coordinate basis. Results are
immutable native-scalar component arrays with multi-index access. Metrics must
be symmetric; raising an index additionally requires invertibility. Coefficients
are differentiated as scalar fields, with native algebraic/trigonometric
simplification. Metric signature, global completeness and coordinate domains
are not inferred.

```python
from symox import sin
from symox.diffgeom import metric_to_Ricci_components

sphere_metric = TensorProduct(dx, dx) + sin(x)**2*TensorProduct(dy, dy)
ricci = metric_to_Ricci_components(sphere_metric)
assert ricci[0, 0] == 1
assert ricci[1, 1] == sin(x)**2
```

`BaseCovarDerivativeOp(coord_sys, index, christoffel)` and
`CovarDerivativeOp(vector, christoffel)` differentiate scalar, vector, covariant
and mixed tensor fields. Supplied Christoffel arrays must have shape
`(dim, dim, dim)` and be expressed in the operator's coordinate basis. A bare
component array carries no coordinate metadata to verify that choice.

`intcurve_series(vector, parameter, start_point, n=6, coord_sys=None,
coeffs=False)` repeatedly applies the vector field to coordinate functions and
evaluates the resulting Taylor terms at the start point. `coeffs=True` returns
the successive matrix terms including parameter powers. `intcurve_diffequ`
returns coordinate ODE residuals and initial-condition residuals for functions
`f_0`, `f_1`, etc.; it does not automatically solve those equations.

Tests cover chart composition/inversion, point evaluation, partial contractions,
graded exterior identities, Jacobi and Cartan identities, metric compatibility,
five independent SymPy curvature comparisons and flat three-dimensional
spherical coordinates. A subprocess verifies calculations while blocking every
SymPy import.
