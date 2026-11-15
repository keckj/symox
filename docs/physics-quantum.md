# Quantum and atomic physics

## Ordered operators and states

`physics.quantum` implements noncommuting operator polynomials, kets, bras,
inner/outer products, adjoints, commutators, anticommutators, tensor products,
wavefunctions and differential operators. `qapply` applies known actions and
retains abstract products when no action is defined. `represent` builds native
matrices; NumPy output is available, and SciPy sparse output needs SciPy.
Hilbert-space objects describe finite complex, function, Fock, direct-sum and
tensor-product spaces.

`represent` follows SymPy's dispatch (`_represent_<Basis>`, the default
basis, then `rep_innerproduct`/`rep_expectation`). `UnitaryOperator`,
`Operator.inverse`, `Commutator`/`AntiCommutator` `doit` and
`expand(commutator=True)` are available, and `hbar` is a `NumberSymbol`.
`quantum.pauli` provides `SigmaX/Y/Z`, `SigmaMinus/Plus`, `SigmaZKet`/`SigmaZBra`
and `qsimplify_pauli`.

These objects are Python wrappers, not native scalar `Expr` subclasses.
Coefficients are native expressions. General symbolic noncommutative calculus,
assumptions and serialization protocols do not reproduce all SymPy behavior.
For example, `TensorProduct` distributes over sums when it is built, and
operator sums can print their terms in another order.

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
angular-momentum operators, rotations and Wigner D elements. `couple` and
`uncouple` use exact Clebsch-Gordan coefficients for any number of spins, with
an explicit `jcoupling` order; with symbolic spins they return the reference
held `Sum`. Rewriting a spin state of symbolic `j` to another axis raises
`NotImplementedError`. Finite representations require explicit
integer/half-integer spin, and spin operators are represented in the `Jz`
basis only. `operatorset` maps supported observables and eigenstate families.

`physics.wigner` implements finite factorial sums for 3j, 6j and 9j symbols,
Clebsch-Gordan and Racah coefficients, Gaunt/real Gaunt coefficients (real
Gaunt with SymPy's U-matrix convention), rotation
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

`physics.secondquant` follows the reference Fermi-vacuum conventions:
bosonic/fermionic creation and annihilation, Fock kets/bras with Fermi-sea
`fermi_level` and holes, operator application, `VarBosonicBasis` and
`FixedBosonicBasis`. Symbols take the `above_fermi`/`below_fermi` assumptions,
and contractions, `NO`, `wicks` and `evaluate_deltas` respect the Fermi level.

```python
from symox import symbols, KroneckerDelta
from symox.physics.secondquant import F, Fd, NO, wicks

a, b = symbols('a b', above_fermi=True)
i, j = symbols('i j', below_fermi=True)
assert wicks(F(a)*Fd(b)) == KroneckerDelta(a, b) + NO(F(a)*Fd(b))
assert wicks(Fd(i)*F(j), keep_only_fully_contracted=True) == KroneckerDelta(i, j)
```

`NO` sorts its operators with the sign and is zero on a Pauli violation;
`NO(...).doit()` splits general dummy indices above and below the Fermi level.
Antisymmetric tensors canonicalize their index order, and `substitute_dummies`,
`simplify_index_permutations` and `PermutationOperator` are the reference
algorithms. Integer orbital labels have no known Fermi level, so
`wicks(F(0)*Fd(0))` keeps a `KroneckerDelta` with a generated dummy, as in
SymPy.

Bosonic `NO(...).doit()`, `wicks` and contractions return the bosonic results
where SymPy raises. A `KroneckerDelta` built by the user does not vanish across
the Fermi level, `Sum` over operator expressions is not supported, and
`NO.subs` with operator keys works only at the top level.

## High-energy gamma tensors

`physics.hep.gamma_matrices` uses the existing indexed tensor engine for
`LorentzIndex` and `GammaMatrix`. Gamma traces recursively implement the
four-dimensional Clifford relations, including contraction with the Minkowski
metric. Contracted gamma pairs and an adjacent identical slashed-vector pair
can be reduced. Fully scalar contracted gamma products return a scalar
multiple of the explicit 4x4 identity.

The dimension is fixed at four. Gamma-five traces, dimensional regularization,
general slashed-vector factorization and full spinor index protocols are not
implemented.
