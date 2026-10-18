# Physics

Physics is implemented locally in `symox.physics`. Commuting scalar expressions,
matrices, derivatives and equations use the native engine. Ordered quantum
products use a Python operator polynomial with native scalar coefficients;
they preserve operator order without treating operators as commuting symbols.
None of these calculations delegates to SymPy at runtime.

The README measures **public name presence** against SymPy 1.14.0. A row at 100%
does not establish complete signature, object protocol or behavioral parity.
The [module inventory](sympy-api.md) includes the deeper import paths as well.

| Guide | Implemented calculations |
|---|---|
| [Mechanics and units](physics-mechanics.md) | SI conversions, moving frames, particles, rigid bodies, Lagrange/Kane equations, structural statics and muscle dynamics |
| [Quantum and atomic physics](physics-quantum.md) | States and ordered operators, qubits, spin, Wigner symbols, atomic wavefunctions, Fock algebra and gamma traces |
| [Optics and control](physics-optics-control.md) | Ray and Gaussian beam propagation, polarization, refraction, transfer functions, state space and response plots |

The implementation targets explicit finite calculations and common symbolic
models. General constrained Kane systems, multispin recoupling, filled Fermi
seas, large circuit simulation and complete structural-mechanics interfaces
remain beyond this implementation. The linked guides describe these limits
next to the supported behavior.

Validation for this implementation: 338 Physics tests pass, including
differential comparisons with SymPy, physical identities, finite quantum
algorithms and checks that block runtime SymPy imports. The complete Python
suite passes 26,566 tests. All eight Python examples in these guides execute.
