# Combinatorics

This module provides native Python combinatorial structures and exact finite
permutation operations, with no SymPy runtime dependency. Every name in the
measured SymPy combinatorics export list is available. This is API presence;
it does not imply complete method parity or large-group algorithms.

## Finite groups and permutations

`SymmetricGroup`, `AlternatingGroup`, `CyclicGroup`, `DihedralGroup`,
`AbelianGroup` and `DirectProduct` construct finite permutation groups.
The lazy `symmetric`, `alternating`, `cyclic` and `dihedral` iterators enumerate
the corresponding permutations in SymPy's order. `SymmetricPermutationGroup`
also exposes exact factorial order and membership without enumeration.

`Permutation` accepts array images or a `Cycle` mapping, with an optional larger
`size`. Multiplication applies the left permutation first. Inversion, integer
powers, cycles, order, parity and signatures are exact. Positional cycle syntax
such as `Permutation(1, 2)` and the wider SymPy Permutation method surface are
not implemented. `Cycle(1, 2, 3)(2, 4)` composes cycles and supports `list(size)`.

Groups expose generators, identity, order, membership, orbit, transitivity,
commutativity, solvability, subgroup construction and element generation.
`Coset(g, H, G=None, dir='+')` supports finite coset lists and membership;
`'+'` denotes the right coset `H*g` and `'-'` the left coset `g*H`.

Generic group operations use bounded exact enumeration, with at most 100,000
elements. Solvability and automatic polycyclic-series construction are limited
to 256 elements. Full Schreier–Sims algorithms are not implemented. Consequently,
large named groups can be constructed while their enumerative operations raise
`NotImplementedError`; use the symbolic symmetric-group interface for its order.

## Ranking and bijections

- `GrayCode(n, start=..., rank=...)` implements binary reflected Gray codes,
  rank/unrank, cyclic advancement, generation and skipping.
- `Subset(subset, superset)` provides binary, Gray and lexicographic ranks and
  neighbours, bitlist conversion and binary/Gray unranking. Elements of the
  superset must be distinct and hashable.
- `Partition` provides disjoint set blocks, restricted growth strings, rank and
  cyclic advancement. `RGS_enum`, `RGS_rank` and `RGS_unrank` use exact completion
  counts; restricted growth strings are limited to length 256.
- `IntegerPartition` accepts a list or multiplicity dictionary, validates an
  optional integer total, and provides lexicographic neighbours, conjugates,
  multiplicities and Ferrers diagrams.
- `Prufer` converts between labelled trees and Prüfer words, with exact rank,
  unrank and cyclic neighbours. Tree inputs are checked for connectivity,
  duplicate edges and valid labels.

Ranking operations do not enumerate preceding subsets, trees or set partitions.
Generator outputs can still be exponentially large when explicitly consumed.

## Free and polycyclic groups

`free_group('x,y')` returns a parent and its free generators. Words reduce
adjacent inverse powers; multiplication, division, inversion, integer powers,
conjugation, commutators, subwords and exponent sums use reduced syllables.
Explicit word powers are limited to 100,000 letters. General finitely presented
group algorithms are not part of this interface.

`group.polycyclic_group()` constructs a certified finite series for supported
solvable permutation groups. `PolycyclicGroup` and `Collector` also accept
explicit generators, series and relative orders. The collector checks normal
forms for uniqueness and coverage, verifies cyclic factors and normality, and
computes power and conjugation relators. It provides exponent vectors, collected
words, depth, leading exponents, induced generating sequences and constructive
membership. Presentations are finite, with at most 100,000 normal forms;
infinite polycyclic presentations are unsupported.

## Polyhedra and Rubik actions

`Polyhedron` tracks named corners, cyclic faces, edges and permutations.
`rotate` mutates the corner arrangement; `reset` restores it and `copy` permits
independent arrangements. `tetrahedron`, `cube`, `octahedron`, `dodecahedron` and
`icosahedron` contain the standard face data. Their rotation groups are generated
by extending oriented face maps across neighbouring faces, yielding orders
12, 24, 24, 60 and 60 respectively.

`RubikGroup(n)` constructs the same `3*(n-1)` slice actions as SymPy on `6*n*n`
facelets, for `2 <= n <= 32`. Construction and composition are supported;
enumerating the Rubik group exceeds the generic group bound.
