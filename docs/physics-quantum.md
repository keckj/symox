# Quantum and atomic physics

## Ordered operators and states

`physics.quantum` implements noncommuting operator polynomials, kets, bras,
inner/outer products, adjoints, commutators, anticommutators, tensor products,
wavefunctions and differential operators. `qapply` applies known actions and
retains abstract products when no action is defined. `represent` builds native
matrices; NumPy output is available, and SciPy sparse output needs SciPy.
Hilbert-space objects describe finite complex, function, Fock, direct-sum and
tensor-product spaces.

These objects are Python wrappers, not native scalar `Expr` subclasses.
Coefficients are native expressions. General symbolic noncommutative calculus,
assumptions and serialization protocols do not reproduce all SymPy behavior.

## Qubits, gates and algorithms

Computational qubits use bit 0 as the rightmost bit. The gate module provides
Pauli, Hadamard, phase, T, swap, controlled, CNOT and user-defined matrix gates.
Their finite matrices and actions agree with that bit convention.

```python
from symox.physics.quantum import qapply, represent
from symox.physics.quantum.qubit import Qubit, measure_all
from symox.physics.quantum.gate import H, CNOT

bell = qapply(CNOT(1, 0)*H(1)*Qubit('00'))
probabilities = measure_all(bell)  # |00> and |11>, each with probability 1/2
```

Finite density operators, von Neumann entropy, full/partial Born measurements,
seeded one-shot sampling, QFT/inverse QFT and Grover search operate on explicit
states. Circuit diagrams draw ordered gates, controls and measurement symbols
with Matplotlib; custom named gates can be drawn without inventing a matrix
for an unspecified operation. Measurement diagram symbols describe a drawing,
not a unitary operator.

`quantum.shor` includes an exact reversible modular-exponentiation oracle,
continued-fraction phase reconstruction and modular-order factorization.
`period_find` is a **classical exact simulation**, not a quantum circuit or a
quantum speedup. It requires time proportional to the modular order.
`shor` returns a pair of factors of a composite integer. Explicit oracle
matrices are limited to 12 qubits. Circuit decompositions and hardware execution
are not implemented.

## Continuous states and angular momentum

Cartesian position/momentum overlaps use Dirac deltas and Fourier kernels;
canonical commutators and wavefunction actions are provided. The infinite
square well (`quantum.piab`) provides normalized sine eigenfunctions and
Hamiltonian energies. Oscillator ladder operators and number states are in
`quantum.sho1d`.

`quantum.spin` provides spin eigenstates, raising/lowering actions, Cartesian
angular-momentum operators, rotations and Wigner D elements. Two-spin coupling
and uncoupling use exact Clebsch–Gordan coefficients. Coupling more than two
momenta, general recoupling trees and symbolic unbounded multiplets are not
implemented. Finite representations require explicit integer/half-integer spin.
`operatorset` maps supported observables and eigenstate families.

`physics.wigner` implements finite factorial sums for 3j, 6j and 9j symbols,
Clebsch–Gordan and Racah coefficients, Gaunt/real Gaunt coefficients, rotation
matrices and spherical-harmonic gradient couplings. Held coefficient wrappers
can be substituted and explicitly evaluated with `doit`; arbitrary expressions
containing them do not implement every reference simplification identity.

## Atomic wavefunctions and matrix identities

`physics.hydrogen`, `qho_1d` and `sho` provide normalized hydrogenic, one- and
three-dimensional harmonic-oscillator wavefunctions, associated energies,
coherent-state coefficients and the hydrogenic Dirac energy. Quantum-number
validation requires the stated physical ranges. Native integration may leave
some normalization integrals held; tests additionally check numerical
normalization independently.

`physics.matrices` implements Pauli matrices, Dirac gamma matrices, the
Minkowski metric, normalized discrete Fourier matrices and the parallel-axis
matrix. `physics.paulialgebra` reduces ordered Pauli products using the exact
Clifford relations, including labels for independent Pauli families.

## Second quantization

`physics.secondquant` implements empty-vacuum bosonic/fermionic creation and
annihilation, Fock kets/bras, operator application, finite bosonic bases and
matrix representations. Fermion states list occupied integer orbitals in
ascending order; creation and annihilation include occupation-parity signs.
Boson states list occupation numbers for each mode.

```python
from symox.physics.secondquant import F, Fd, wicks

assert wicks(F(0)*Fd(0) + Fd(0)*F(0)) == 1
assert wicks(F(0)*F(1)*Fd(0)*Fd(1),
             keep_only_fully_contracted=True) == -1
```

Normal ordering and Wick expansion evaluate finite ladder polynomials and
fully contracted vacuum expectations. `NO(...).doit()` discards contractions;
`wicks` includes them. Antisymmetric tensors canonicalize upper/lower index
permutations, delta simplification substitutes repeated indices, and explicit
permutation operators antisymmetrize expressions.

Filled Fermi-sea assumptions, symbolic occupied-orbital ordering, mixed
boson/fermion normal ordering, general diagram expansions and complete dummy
index equivalence classes are not implemented. Dummy substitution renames
native dummy symbols deterministically; it is not a complete canonicalizer for
all abstract many-body expressions.

## High-energy gamma tensors

`physics.hep.gamma_matrices` uses the existing indexed tensor engine for
`LorentzIndex` and `GammaMatrix`. Gamma traces recursively implement the
four-dimensional Clifford relations, including contraction with the Minkowski
metric. Contracted gamma pairs and an adjacent identical slashed-vector pair
can be reduced. Fully scalar contracted gamma products return a scalar
multiple of the explicit 4×4 identity.

The dimension is fixed at four. Gamma-five traces, dimensional regularization,
general slashed-vector factorization and full spinor index protocols are not
implemented.
