# Coordinate vector calculus

[README](../README.md) · [Geometry](geometry.md) · [Holonomic functions](holonomic.md)

`symox.vector` provides all 43 exports from the corresponding SymPy 1.14
package. This measures importable names; the supported operations and limits
are documented here. It is separate from the `symox.physics.vector` package
and its mechanics reference frames, described in
[Units and mechanics](physics-mechanics.md).
The conventions follow the [SymPy vector documentation](https://docs.sympy.org/latest/modules/vector/fields.html).
Python objects store geometric structure; scalar arithmetic, matrices,
differentiation, solving and integration use the native engine without
importing SymPy.

## Vectors, dyadics and connected frames

```python
from symox.vector import CoordSys3D, Vector, gradient, divergence, curl

C = CoordSys3D('C')
field = C.x*C.i + C.y*C.j + C.z*C.k
assert divergence(field) == 3
assert curl(gradient(C.x*C.y)) == Vector.zero
assert C.i.cross(C.j) == C.k
assert C.i.dot(C.j) == 0
```

`Vector` expressions are finite maps from base vectors to native scalar
coefficients. Addition, scalar products, division by scalars, substitution,
coefficient differentiation/integration, simplification, magnitude,
normalization, projection and `to_matrix(system)` are available. Use `.dot`
(`&`), `.cross` (`^`) and `.outer` (`|`) for vector products. Multiplication
of two vectors with `*` is rejected. `VectorAdd`, `VectorMul`, `VectorZero`,
`BaseVector` and `VectorKind` expose the corresponding constructors/metadata;
`Add` of two vectors builds a `VectorAdd`. `integrate` and `Integral` act
componentwise. As in SymPy, `equals` compares `(a - b).dot(a - b)` with zero,
so equal vectors written in rotated frames compare equal, and
`projection(other, scalar=True)` is the coefficient `dot(other)/dot(self)`.
Vectors and dyadics print in the reference form.

Dyadics use the same finite-component representation. `BaseDyadic`,
`DyadicAdd`, `DyadicMul` and `DyadicZero` support construction; dot/cross
contractions, transpose, scalar transforms and two-frame matrix conversion
operate on the actual components.

`BaseScalar` creates native symbols and registers their coordinate-system
identity so scalar calculus retains native differentiation and integration.
Distinct connected frames with the same printed name receive distinct scalar
atoms; a later colliding name can display a native dummy-symbol prefix.
These atoms are not Python subclasses with all SymPy BaseScalar attributes.

`CoordSys3D` supplies base scalars/vectors, an origin, a connected parent,
translation and rotations. `locate_new`, `orient_new`,
`orient_new_axis/body/space/quaternion`, `rotation_matrix`, `scalar_map`,
`position_wrt` and transformations to/from a parent implement these relations.
`express(..., variables=True)` changes both bases and coordinate scalars;
without that flag it changes bases only. Dyadics can use a second target frame.
Disconnected frames raise an error.

`AxisOrienter`, `BodyOrienter`, `SpaceOrienter` and `QuaternionOrienter` produce
direction-cosine matrices. The twelve Euler/Tait-Bryan axis sequences and their
body/space multiplication conventions are checked against SymPy. Axis vectors
must be nonzero. As in the reference, quaternions (including symbolic ones)
are used as given, without a unit-norm check, and a supplied 3 by 3 rotation
matrix is stored without an orthogonality check; non-unit or non-orthogonal
input gives a non-rigid frame. `Point` connects displacement
vectors, supports `locate_new`, `position_wrt` and `express_coordinates`.

## Orthogonal coordinates and field operators

Predefined Cartesian, cylindrical and spherical charts supply their metric
(Lamé) coefficients. `create_new` connects a curvilinear chart to a Cartesian
parent. A callable or `(variables, expressions)` transformation can define a
custom orthogonal chart. Nonorthogonal/degenerate charts are rejected; inversion
of a custom chart is not implemented. Predefined chart inverses use principal
`atan2`/`acos` branches. Coordinate-domain restrictions and global charts are
not represented; cylindrical/spherical operators apply away from chart
singularities on the conventional positive-radius angular domains.

```python
P = CoordSys3D('P', transformation='cylindrical')
assert gradient(P.r**2) == 2*P.r*P.i
assert divergence(P.r*P.i) == 2
assert curl(P.r*P.j) == 2*P.k
```

`gradient`, `divergence`, `curl`, scalar/vector `laplacian` and
`directional_derivative` use the metric and coordinate derivatives.
Directional derivatives of vector fields include derivatives of moving basis
vectors, including parented cylindrical/spherical charts. In contrast,
`vector.diff(parameter)` differentiates coefficients with bases held fixed.

`Gradient`, `Divergence`, `Curl`, `Laplacian`, `Dot` and `Cross` are held
operation wrappers with `.doit()`. `gradient/curl/divergence(..., doit=False)`
and the `Del` operator return these wrappers. They are not a complete algebra
of unevaluated vector differential expressions; evaluate them before combining
them with ordinary vector expressions. Vector/dyadic scalar transforms are
available directly as methods; public `symox.diff/integrate/simplify/expand`
also dispatch after the vector package is imported.

`is_conservative` tests zero curl locally; it does not infer domain topology
or global single-valuedness. `is_solenoidal` tests zero divergence.
`scalar_potential` integrates compatible metric components and checks the
resulting gradient; `scalar_potential_difference` evaluates a potential at two
connected points. Unsupported scalar integration can prevent reconstruction.

## Regions and geometric integrals

`ParametricRegion(definition, *bounds)` describes one to three coordinates
using distinct parameters, with bounds `(parameter, lower, upper)` or unbounded
parameters treated as symbolic constants. Its integration dimension is the
number of bounded parameters. Bounds can depend on outer parameters; cyclic
dependencies are rejected.

`ParametricIntegral(field, region)` and `vector_integrate(field, region)`
compute scalar arc length/surface measure, vector circulation/flux, or scalar
volume integrals. Parameter tangents and cross products give the metric and
oriented surface normal. Reversing a surface parametrization reverses vector
flux. Volumes use the signed Jacobian convention. Scalar volume integrands
are required. Curvilinear region definitions are mapped to physical Cartesian
positions before computing tangents. Unsupported antiderivatives remain
native unevaluated integrals.

```python
from symox import Symbol, Rational
from symox.vector import ParametricRegion, vector_integrate

t = Symbol('t', real=True)
curve = ParametricRegion((3*t, 4*t), (t, 0, 1))
assert vector_integrate(1, curve) == 5
assert vector_integrate(C.x*C.i, ParametricRegion((t, t**2), (t, 0, 1))) == Rational(1, 2)
```

`parametric_region_list` adapts geometry points, segments, curves, ellipses
and polygon boundaries. `ImplicitRegion(variables, equation)` accepts an
expression or equality and exposes its polynomial degree. `regular_point`,
`singular_points`, `multiplicity` and `rational_parametrization` are ports of
the reference algorithms; conic regular points use the reference parametric
solution of a ternary quadratic form. For a plane in three variables,
`rational_parametrization` solves the linear equation, where the reference
raises `NotImplementedError`. In three variables `regular_point` skips grid
points where the gradient vanishes and raises `NotImplementedError` when every
candidate is singular (the reference can return a singular point there),
rather than asserting the locus is empty. Positive-dimensional singular loci
are not represented. Automatic bounded integration of implicit regions supports
axis-aligned ellipses/ellipsoids; other cases require an explicit parametric
region with chosen bounds. Arbitrary algebraic-region integration is not
implemented.

Metric simplification handles standard angular ranges such as `0 <= theta <= pi`
for sphere area. General absolute-value integration, parametrization singularities,
multiple coverings and arbitrary branch/domain analysis remain subject to the
scalar engine and the supplied region.

Tests check exact product and differential identities, body/space orientation
conventions, translated/rotated scalar maps, moving bases, geometric boundary
integrals, surface area and the divergence theorem on a sphere. Separate
subprocess tests block SymPy imports during actual vector calculations.
