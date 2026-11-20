# Physics

Physics is implemented locally in `symox.physics`. Commuting scalar expressions,
matrices, derivatives and equations use the native engine. Ordered quantum
products use a Python operator polynomial with native scalar coefficients;
they preserve operator order without treating operators as commuting symbols.
`physics.vector`, `physics.mechanics`, `physics.control.lti` and the `Beam`,
`Truss` and `Arch` structural modules are the reference SymPy 1.14 sources
rewritten to import symox. None of these calculations delegates to SymPy at
runtime.

The README measures **public name presence** against SymPy 1.14.0. A row at 100%
does not establish complete signature, object protocol or behavioral parity.
The [module inventory](sympy-api.md) includes the deeper import paths as well.

| Guide | Implemented calculations |
|---|---|
| [Mechanics and units](physics-mechanics.md) | SI conversions, moving frames, particles, rigid bodies, Lagrange/Kane equations, structural statics and muscle dynamics |
| [Quantum and atomic physics](physics-quantum.md) | States and ordered operators, qubits, spin, Wigner symbols, atomic wavefunctions, Fock algebra and gamma traces |
| [Optics and control](physics-optics-control.md) | Ray and Gaussian beam propagation, polarization, refraction, transfer functions, state space and response plots |

The implementation targets explicit finite calculations and common symbolic
models. Large circuit simulation, `Sum` over quantum operator expressions and
the remaining printing and ordering differences recorded in the
[parity notes](algorithmic-parity.md) remain beyond this implementation. The
linked guides describe these limits next to the supported behavior.

Validation: the `tests/test_physics_*.py` files check physical identities,
finite quantum algorithms and runtime runs with SymPy imports blocked. The
wave-4 parity files (`tests/test_wrong_*_wave4_parity.py`,
`tests/test_missing_*_wave4_parity.py`) hard-code SymPy 1.14 outputs for
mechanics, control, quantum and second quantization. The Python examples in
these guides execute.
