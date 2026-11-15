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

`physics.vector` and `physics.mechanics` are the reference SymPy 1.14 sources,
rewritten to import symox; their scalar, matrix, solving and simplification
calls run on the native engine.

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
Frame coordinates are `CoordinateSym` symbols that print as `A_x` and carry
their `.frame`; coordinates of two frames with the same name stay distinct.
`msubs` follows the reference and does not substitute inside derivatives.

`physics.mechanics` provides particles, rigid bodies, `Body`, inertia dyadics,
momenta, energies, forces and torques. `LagrangesMethod` derives Euler-Lagrange
equations, including applied loads, holonomic/nonholonomic constraints and
Lagrange multipliers. `KanesMethod` derives generalized active and inertia
forces from partial velocities, including dependent and auxiliary speeds and
dependent coordinates. Both expose mass/forcing matrices and first-order
right-hand sides.

```python
from symox import symbols, diff
from symox.physics.mechanics import LagrangesMethod, dynamicsymbols

m, k = symbols('m k', positive=True)
q = dynamicsymbols('q')
t = dynamicsymbols._t
model = LagrangesMethod(m*diff(q, t)**2/2-k*q**2/2, [q])
model.form_lagranges_equations()  # m*q'' + k*q = 0
```

`Linearizer` and the `linearize` methods use native Jacobians. `System`
assembles bodies, joints, actuators and loads and forms equations of motion
with `form_eoms`; `JointsMethod` is also available. Pin, prismatic,
cylindrical, planar, spherical and weld joints accept the reference
`parent_point`, `child_point`, `parent_axis`, `child_axis`, `parent_interframe`
and `child_interframe` arguments and generate kinematic differential
equations. Pathways include linear, obstacle-set and wrapping paths on spheres
and cylinders. Springs, dampers, Duffing springs, Coulomb friction and
force/torque actuators produce actual loads. `models` provides the reference
`multi_mass_spring_damper` and `n_link_pendulum_on_cart` examples.

Known differences from SymPy are listed in the
[parity notes](algorithmic-parity.md#mechanics-vector-and-continuum-mechanics-wp-c12-wave-4):
mostly equal closed forms printed differently, set and term ordering, and the
reprs of plot objects.

## Structural mechanics

`physics.continuum_mechanics` implements the following. `Beam`, `Beam3D`,
`Truss` and `Arch` are the reference SymPy 1.14 sources on symox; `Cable` is a
native implementation.

- `Beam`: point forces/couples and polynomial distributed loads, singularity
  functions, pin/roller/fixed supports, rotation and sliding hinges, joined
  beams, reaction solving, shear, bending moment, slope, deflection, their
  `max_*` values, influence lines and plots. Distributed loads can terminate
  at an endpoint.
- `Beam3D`: the reference Timoshenko formulation (shear and bending equations
  solved with `dsolve`), axial extension and torsion.
- `Truss`: planar pin-jointed joint equilibrium with pinned/roller supports,
  member forces, support reactions and drawing.
- `Cable`: polygonal point-loaded equilibrium with specified node positions,
  and parabolic cables under uniform load per horizontal span. A uniform-load
  solution needs lowest-point geometry to determine sag.
- `Arch`: three-hinged parabolic arches with hinge or roller supports, an
  optional tie member, point and distributed loads, support reactions and
  internal bending, shear and axial forces.

Import the arch class with `from symox.physics.continuum_mechanics.arch import Arch`.
The module name is lowercase `arch`; `Arch` is the class name.

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

`Cable` deliberately differs from the reference where the reference
equilibria are wrong: distributed reactions stay balanced, `apply_length`
compares with the true support distance, and `change_support` checks load
positions instead of raising `TypeError`. Cable loading assumes either point
loads or uniform loading, not both. `Beam.max_*` with symbolic `E` and `I` can
be very slow, as in SymPy.

## Biomechanics

`physics.biomechanics` provides the seven De Groote 2016 characteristic curves,
their explicit derivatives and published default parameters. Passive fiber,
tendon and force-velocity curves have inverse functions.
`ZerothOrderActivation` equates activation and excitation;
`FirstOrderActivationDeGroote2016` implements the smooth first-order activation
ODE with distinct activation/deactivation time constants.

`MusculotendonDeGroote2016` couples those curves and activation to a mechanical
pathway. The rigid-tendon model produces tensile loads. The compliant-tendon
models use normalized fiber length or normalized tendon force as a state.
Explicit formulations invert the force-velocity curve and require zero fiber
damping. Implicit formulations expose normalized fiber velocity as an input
and the force-balance residual through `algebraic_equations`; supply these to
a differential-algebraic solver. They do not silently solve an algebraic
constraint as an explicit ODE. The rigid model uses fixed optimal pennation,
matching the reference model's convention. Elastic formulations use constant
fiber-height geometry.
