# Optics and control

## Optics

`physics.optics` provides ABCD ray matrices, geometric rays, free propagation,
flat/curved refraction and mirrors, thin lenses and Gaussian beam parameters.
Ray and beam transformations use actual matrix composition. Gaussian beam
objects expose waist, Rayleigh range, width, radius of curvature, divergence
and Gouy phase. Conjugate helpers solve geometric and Gaussian imaging.

```python
from symox import Rational
from symox.physics.optics import FreeSpace, ThinLens, GeometricRay

ray = GeometricRay(1, 0)
result = FreeSpace(2)*ThinLens(4)*ray
assert result.height == Rational(1, 2)
assert result.angle == Rational(-1, 4)
```

The Gaussian conjugate calculation follows ABCD propagation:
`z_R_out = z_R*f**2 / ((s-f)**2 + z_R**2)`. This deliberately corrects an
inconsistent denominator in the SymPy 1.14 helper; it is tested against beam
propagation rather than reproducing that helper's numerical result.

Jones vectors, Stokes vectors, polarizers, retarders, beam splitters and
Jones-to-Mueller conversion model polarization. `Medium` relates refractive
index, permittivity, permeability, impedance and wave speed using SI quantities.
Snell, deviation, Fresnel amplitude, Brewster/critical-angle, lensmaker,
imaging, hyperfocal and transverse-magnification helpers cover scalar and
explicit vector calculations. `TWave` models symbolic transverse waves and
coherent superposition.

Ray wrappers expose their matrix data but are not native Matrix subclasses.
Vector refraction returns components rather than preserving all Ray3D geometry
metadata. Total internal reflection follows the reference zero-result
convention for refracted rays; complex Fresnel continuation and every
polarizing splitter phase convention are not exhaustively covered.

## Linear control

`physics.control` provides SISO transfer functions, PID controllers, series,
parallel and feedback connections, transfer-function matrices and MIMO
interconnections. Rational arithmetic is symbolic and exact when inputs are
exact. Pole/zero queries, properness, stability, DC/frequency response and
phase/gain margins use native polynomial and scalar calculations.

State-space models validate matrix dimensions and expose controllability,
observability, transfer-function conversion, interconnections and state
solutions. Transfer-function conversion to state space uses a controllable
companion realization. Generalized bilinear, bilinear, forward-difference and
backward-difference transforms return discrete polynomial coefficients.

```python
from symox import symbols
from symox.physics.control import TransferFunction, Feedback

s = symbols('s')
plant = TransferFunction(1, s+1, s)
closed_loop = Feedback(plant).doit()  # 1/(s+2)
```

Pole-zero, impulse, step, ramp, Bode, Nyquist and Nichols numerical data and
Matplotlib plotting functions are available. Time-domain responses use native
inverse Laplace transforms; sampling uses native numerical evaluation, NumPy
arrays and optional Matplotlib. SciPy is not required.

Numerical responses require all parameters other than time/frequency to be
specified. Unsupported inverse transforms (for example some delays or high
order symbolic poles) raise explicitly. Distributional impulse feedthrough is
not rendered as a finite sample. The supported MIMO operations do not reproduce
every fluent operator combination and realization choice in SymPy. Empty
crossover sets give infinite margins; multibranch and transcendental crossover
searches are not a general robust-control analysis engine.
