# Exact Euclidean geometry

[README](../README.md)

`symox.geometry` exports the 30 names in SymPy 1.14's geometry namespace.
This is API name presence. The implementations below cover exact finite
Euclidean constructions; they do not implement every SymPy geometry method.
Geometry objects are immutable Python containers whose coordinates and
calculations use the native Rust scalar engine. Runtime calculations never
import SymPy. Explicit `to_sympy()` is an optional interchange operation.

```python
from symox import Point, Line, Segment, Circle, Triangle, Plane, Point3D

assert Circle((0, 0), 2).intersection(Line((0, 0), (1, 0))) == [Point(-2, 0), Point(2, 0)]
triangle = Triangle((0, 0), (4, 0), (0, 3))
assert triangle.area == 6
assert triangle.incenter == Point(1, 1)
assert triangle.inradius == 1
plane = Plane((0, 0, 0), normal_vector=(0, 0, 1))
assert plane.projection(Point3D(1, 2, 3)) == Point3D(1, 2, 0)
```

## Points and linear entities

`Point`, `Point2D` and `Point3D` accept concrete coordinate sequences or
native symbolic coordinates. Finite coordinates must not be known nonreal.
Default float coordinates are rationalized; `evaluate=False` retains floats.
Points support vector arithmetic, dot products, distances, midpoints, unit
vectors, translations and scaling. Planar rotations and homogeneous
transformation matrices are supported. Spatial points additionally offer a
native cross product and coplanarity checks. As in SymPy, a point's `length`
is zero; `abs(point)` is its Euclidean vector norm.

`Line`, `Ray` and `Segment` select their 2D or 3D variant from the input
dimension. Constructors accept endpoints; lines and rays also accept a slope
or direction ratio. Operations include containment, geometric equality,
parallelism, perpendicularity, angles, projections, point distances and exact
intersections. Collinear overlap returns a point, segment, ray or line as
appropriate. Spatial skew intersections are rejected; segment and ray
distances include interior closest-point candidates. A segment with equal
endpoints becomes that point. A line or ray with equal endpoints is invalid.

## Conics, polygons and planes

`Circle` accepts a center/radius or three noncollinear points. `Ellipse`
accepts a center and horizontal/vertical radii, or a radius and eccentricity.
They expose area, circumference, foci, equations, parametric points, tangent
lines, containment and line/ray/segment intersections. Circle-circle
intersections retain exact radicals. Equal radii produce a circle. Opposite
signed circle radii define the same locus; their intersection returns a
coincident circle rather than reproducing SymPy 1.14's parametrized points.
The point-to-circle distance is also available as a native extension.

`Polygon` provides signed area, perimeter, centroid, sides, angles, convexity,
boundary containment, interior point tests and boundary intersections.
Consecutive repeated vertices and redundant collinear vertices are removed.
Three remaining vertices produce a `Triangle`; one or two produce a point or
segment. `Triangle` adds medians, altitudes, bisectors, orthocenter,
circumcircle, incircle and standard shape/similarity predicates.
`RegularPolygon` stores its center, radius, vertex count and rotation and
generates exact vertices when requested.

`Plane` accepts three noncollinear points or a point and normal vector. It
supports equations, containment, projection, distances, angles, parallelism,
perpendicularity and spatial line/plane intersections. `Parabola` stores a
focus and planar directrix, exposes its vertex, axis, focal length and
equation, and intersects planar linear entities. `Curve` stores two or three
parametric coordinates and a finite parameter interval, with native arc
length integration and point transformations.

## Utilities and limits

The namespace includes `rad`, `deg`, `intersection`, weighted `centroid`,
`convex_hull`, closest/farthest point pairs, `are_similar`, and
`idiff` for repeated implicit differentiation. Geometry objects support
substitutions, exact reconstruction, numerical evaluation, hashing and
pickling.

Intersections and containment that require an undecidable real ordering raise
`NotImplementedError`; conditional geometric sets are not invented. General
ellipse-ellipse intersections, rotated implicit ellipses, spatial rotations
without an explicit transformation matrix, and general conic-conic
intersections remain unimplemented. Parabola scaling requires a nonsingular
similarity. Polygon region boolean operations and holes are not implemented;
polygon intersections concern boundaries. Similarity utilities currently
cover points, linear entities, ellipses and finite polygons (up to 256
vertices for general polygon similarity). General geometry printers and code generation
remain outside this layer.

Points currently support dimensions two and three. Polygons are limited to
4096 vertices, distance-pair searches to 1024 points and implicit derivatives
to order 256. Exact signs and simplification follow the native scalar
engine's supported domains and assumptions.
