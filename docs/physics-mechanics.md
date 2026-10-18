# Units and mechanics

## Units

`physics.units` provides the 361 public bindings of the SymPy 1.14.0 SI catalog:
quantities, physical constants, dimensions, prefixes and aliases. The catalog
contains 127 distinct quantities. Conversions solve dimensional exponent
equations over exact rational coefficients; numerical constants retain the
reference catalog's precision.

```python
from symox.physics.units import kilometer, meter, hour, second, newton, kilogram
from symox.physics.units import convert_to

assert convert_to(kilometer/hour, meter/second) == 5*meter/(18*second)
assert convert_to(newton, [kilogram, meter, second]) == kilogram*meter/second**2
```

`Quantity`, `PhysicalConstant`, `Prefix`, `Dimension`, `DimensionSystem` and
`UnitSystem` support custom units and systems. Use `quantity_simplify` for
unit-aware simplification and `check_dimensions` to reject sums of incompatible
known dimensions. Converting to an incompatible target leaves the input
unchanged, following the reference convention. Unknown ordinary symbols have
no inferred physical dimension. Affine temperature scales are not supported.
The ordinary expression simplifier does not automatically normalize units.

## Moving frames and rigid bodies

`physics.vector` provides `ReferenceFrame`, `Point`, `Vector`, `Dyadic`,
`dynamicsymbols`, differential field operations and vector printing. Frames
support axis, body, space, quaternion and explicit direction cosine matrix
orientations. Frame graphs transform components and angular velocities;
point graphs propagate positions, velocities and accelerations.

```python
from symox import diff, sin, cos
from symox.physics.vector import ReferenceFrame, dynamicsymbols

q = dynamicsymbols('q')
N = ReferenceFrame('N')
A = N.orientnew('A', 'Axis', (q, N.z))
assert A.x.express(N) == cos(q)*N.x + sin(q)*N.y
assert A.x.dt(N) == diff(q, dynamicsymbols._t)*A.y
```

Dot, cross and outer products, vector/dyadic time derivatives, partial
velocities, motion integration, potentials and kinematic equations operate on
native components. Quaternion orientations assume a unit quaternion.
`CoordinateSym` uses native dummy symbols; display and class hierarchies can
differ from SymPy.

`physics.mechanics` provides particles, rigid bodies, inertia dyadics,
momenta, energies, forces and torques. `LagrangesMethod` derives Euler–Lagrange
equations, including applied loads, holonomic/nonholonomic constraints and
Lagrange multipliers. `KanesMethod` derives generalized active and inertia
forces from partial velocities for particles and rigid bodies. Both expose
mass/forcing matrices and first-order right-hand sides.

```python
from symox import symbols, diff
from symox.physics.mechanics import LagrangesMethod, dynamicsymbols

m, k = symbols('m k', positive=True)
q = dynamicsymbols('q')
t = dynamicsymbols._t
model = LagrangesMethod(m*diff(q, t)**2/2-k*q**2/2, [q])
model.form_lagranges_equations()  # m*q'' + k*q = 0
```

Linearization uses native Jacobians. `System` assembles bodies, joints,
actuators and loads; pin, prismatic, cylindrical, planar, spherical and weld
joints generate kinematic differential equations. Pathways include linear and
multi-segment paths, plus spherical and cylindrical wrapping. Springs, dampers,
Duffing springs, friction and force/torque actuators produce actual loads.

Dependent or auxiliary Kane speeds/coordinates, constrained independent-state
linearization and some advanced joint alignment choices are not implemented.
These cases raise explicitly rather than producing incomplete equations.

## Structural mechanics

`physics.continuum_mechanics` implements:

- `Beam`: point forces/couples and polynomial distributed loads, singularity
  functions, pin/roller/fixed supports, reaction solving, shear, bending moment,
  slope, deflection and plots. Distributed loads can terminate at an endpoint.
- `Beam3D`: independent bending planes, axial extension and Saint-Venant torsion.
- `Truss`: planar pin-jointed joint equilibrium with pinned/roller supports,
  member forces and support reactions. Unstable or statically indeterminate
  systems are rejected.
- `Cable`: polygonal point-loaded equilibrium with specified node positions,
  and parabolic cables under uniform load per horizontal span. A uniform-load
  solution needs lowest-point geometry to determine sag.
- `Arch`: three-hinged parabolic arches, point and uniform distributed loads,
  support reactions and internal bending, shear and axial forces.

```python
from symox import Rational
from symox.physics.continuum_mechanics import Beam

beam = Beam(10, 2, 3)
r0 = beam.apply_support(0, 'pin')
r10 = beam.apply_support(10, 'roller')
beam.apply_load(-12, 5, -1)
beam.solve_for_reaction_loads(r0, r10)
assert beam.reaction_loads == {r0: 6, r10: 6}
assert beam.deflection().subs(beam.variable, 5) == Rational(-125, 3)
```

Structural helpers cover the stated models, not every SymPy method. Beam hinges, composite sections, structural
optimization and tied/roller-supported arches are not implemented. `Beam3D`
is an uncoupled linear model, without coupled warping or shear deformation.
Cable loading assumes either point loads or uniform loading, not both.

## Biomechanics

`physics.biomechanics` provides the seven De Groote 2016 characteristic curves,
their explicit derivatives and published default parameters. Passive fiber,
tendon and force–velocity curves have inverse functions.
`ZerothOrderActivation` equates activation and excitation;
`FirstOrderActivationDeGroote2016` implements the smooth first-order activation
ODE with distinct activation/deactivation time constants.

`MusculotendonDeGroote2016` couples those curves and activation to a mechanical
pathway. The rigid-tendon model produces tensile loads. The compliant-tendon
models use normalized fiber length or normalized tendon force as a state.
Explicit formulations invert the force–velocity curve and require zero fiber
damping. Implicit formulations expose normalized fiber velocity as an input
and the force-balance residual through `algebraic_equations`; supply these to
a differential-algebraic solver. They do not silently solve an algebraic
constraint as an explicit ODE. The rigid model uses fixed optimal pennation,
matching the reference model's convention. Elastic formulations use constant
fiber-height geometry.
