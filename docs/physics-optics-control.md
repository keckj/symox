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
`conjugate_gauss_beams` also accepts `s_in` (the focal length that places the
output waist for a given input distance), where SymPy 1.14 raises
`NotImplementedError`.

Jones vectors, Stokes vectors, polarizers, retarders, beam splitters and
Jones-to-Mueller conversion model polarization. `Medium` relates refractive
index, permittivity, permeability, impedance and wave speed using SI quantities.
Snell, deviation, Fresnel amplitude, Brewster/critical-angle, lensmaker,
imaging, hyperfocal and transverse-magnification helpers cover scalar and
explicit vector calculations. `TWave` models symbolic transverse waves and
coherent superposition.

`physics.optics.utils` is a port of the reference module: `refraction_angle`
returns a `Ray3D` for a `Ray3D` input, raises `ValueError` on total internal
reflection for an angle input and returns 0 for a vector input, as SymPy does.
Under total internal reflection `fresnel_coefficients` returns the two complex
reflection coefficients. Ray wrappers expose their matrix data but are not
native Matrix subclasses. Every polarizing splitter phase convention is not
exhaustively covered.

## Linear control

`physics.control.lti` is the reference SymPy 1.14 module running on symox. It
provides SISO transfer functions (including `TransferFunction.from_zpk`), PID
controllers, series, parallel and feedback connections (also with
`StateSpace` operands), transfer-function matrices (`from_Matrix`,
`elem_poles`, `elem_zeros`, `eval_frequency`, matrix-style indexing) and the
MIMO interconnections. Rational arithmetic is symbolic and exact when inputs
are exact. Pole/zero queries, properness, stability, DC/frequency response and
phase/gain margins use native polynomial and scalar calculations.

State-space models validate matrix dimensions and expose controllability,
observability, controllable and observable subspaces, transfer-function
conversion, interconnections and state solutions (`dsolve` returns the output
`C*x + D*u`). Generalized bilinear, bilinear, forward-difference and
backward-difference transforms return discrete polynomial coefficients.

```python
from symox import symbols
from symox.physics.control import TransferFunction, Feedback

s = symbols('s')
plant = TransferFunction(1, s+1, s)
closed_loop = Feedback(plant).doit()  # (s + 1)/((s + 1)*(s + 2)), as in SymPy
```

Pole-zero, impulse, step, ramp, Bode, Nyquist and Nichols numerical data and
Matplotlib plotting functions are available. Time-domain responses use native
inverse Laplace transforms; sampling uses native numerical evaluation, NumPy
arrays and optional Matplotlib. SciPy is not required.

Numerical responses require all parameters other than time/frequency to be
specified. Unsupported inverse transforms (for example some delays or high
order symbolic poles) raise explicitly. Distributional impulse feedthrough is
not rendered as a finite sample.

`phase_margin` and `gain_margin` follow the reference rules: the phase margin
is folded below 180 degrees and is `-180` without a gain crossover, and the
gain margin is `oo` without a phase crossover. Symox solves `|G|**2 - 1`
instead of SymPy's logarithmic equation and decides crossover membership and
the `>= 180` test numerically at 50 digits; a tie it cannot decide raises
rather than guessing. Results can be equal closed forms printed differently
from SymPy. Non-SISO systems, time delays and extra symbols raise the
reference errors. This is not a general robust-control analysis engine.
