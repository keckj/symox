# Combinatorics

This module provides Python combinatorial structures and exact finite
permutation operations, with no SymPy runtime dependency. Every name in the
measured SymPy combinatorics export list is available. This is API presence;
it does not imply complete method parity. The group modules (`perm_groups`,
`free_groups`, `fp_groups`, `coset_table`, `homomorphisms`, `pc_groups`,
`rewritingsystem`, `tensor_can`, `util`, `testutil`) are ports of the SymPy
1.14 algorithms. Everything in this package runs in Python; there is no Rust
kernel for permutation groups. `Permutation` and `PermutationGroup` are not
subclasses of the scalar `Basic` hierarchy, so generic expression methods
such as `subs`, `xreplace` or `has` are absent.

## Finite groups and permutations

`SymmetricGroup`, `AlternatingGroup`, `CyclicGroup`, `DihedralGroup`,
`AbelianGroup` and `DirectProduct` construct finite permutation groups.
The lazy `symmetric`, `alternating`, `cyclic` and `dihedral` iterators enumerate
the corresponding permutations in SymPy's order. `SymmetricPermutationGroup`
also exposes exact factorial order and membership without enumeration.

`Permutation` accepts array images, cycle lists, positional cycle syntax
(`Permutation(1, 2)`, `Permutation(3)(0, 1)`) or a `Cycle` mapping, with an
optional larger `size`. Multiplication applies the left permutation first.
Inversion, integer powers, cycles, order, parity, signatures, ranks and the
other methods of the SymPy `Permutation` class are available, and the repr is
the reference cyclic form. `Cycle` is a port of the reference class:
`Cycle(1, 2, 3)(2, 4)` composes cycles, `list(Cycle(1, 2))` is the array form
and `list(size)` pads it.

Order and membership use an exact incremental stabilizer chain (Schreier-Sims),
without enumerating elements, so `RubikGroup(3).order()` and membership in
`SymmetricGroup(20)` are immediate. Bases, strong generating sets,
transversals and the structural algorithms are the reference ones: centers,
centralizers, normal closures, derived, lower central and composition series,
solvability, primitivity and blocks, Sylow subgroups, subgroup search, coset
factorization and ranking, presentations and abelian invariants.
`Coset(g, H, G=None, dir='+')` supports finite coset lists and membership;
`'+'` denotes the right coset `H*g` and `'-'` the left coset `g*H`.

As in SymPy, `normal_closure`, `sylow_subgroup`, `random_stab`,
`schreier_sims_random`, `polycyclic_group` and the presentations are
randomized, so their generators can differ between runs, and
`schreier_sims_random` is a Monte Carlo algorithm. `generate` and
`elements` enumerate the whole group and are only practical for groups of
moderate order.

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

## Free, finitely presented and polycyclic groups

`free_group('x,y')` returns a parent and its free generators. Words reduce
adjacent inverse powers; multiplication, division, inversion, integer powers,
conjugation, commutators, subwords and exponent sums follow the reference
`free_groups` module. Free group words are pure Python.

`FpGroup(F, relators)` provides the reference finitely presented group
algorithms: Todd-Coxeter coset enumeration (`coset_table`), order and index
computations, presentation simplification, subgroup presentations, rewriting
systems and homomorphisms. Abelian invariants use an exact integer Smith
normal form. Algorithms that need a finite group (permutation presentations,
kernels, abelian invariants, composition series) raise `NotImplementedError`
on infinite groups. Coset enumeration of hard presentations is as slow as in
the reference.

`group.polycyclic_group()` builds a polycyclic sequence of a solvable
permutation group from its derived series, as the reference does; the chosen
sequence is randomized. `PolycyclicGroup` and `Collector` also accept explicit
generators, series and relative orders, and provide power and conjugation
relators, exponent vectors, collected words, depth, leading exponents, induced
generating sequences and constructive membership. Three reference results are
deliberately corrected: `collected_word` keeps the whole power relator,
`induced_pcgs` inverts leading exponents modulo the relative order, and
`constructive_membership_test` returns the exponent vector of a member
instead of `False`. A leading exponent that is not invertible modulo its
relative order raises `NotImplementedError`.

## Polyhedra and Rubik actions

`Polyhedron` tracks named corners, cyclic faces, edges and permutations.
`rotate` mutates the corner arrangement; `reset` restores it and `copy` permits
independent arrangements. `tetrahedron`, `cube`, `octahedron`, `dodecahedron` and
`icosahedron` contain the standard face data. Their rotation groups are generated
by extending oriented face maps across neighbouring faces, yielding orders
12, 24, 24, 60 and 60 respectively.

`RubikGroup(n)` constructs the same `3*(n-1)` slice actions as SymPy on `6*n*n`
facelets, for `2 <= n <= 32`. Its order and membership come from the
stabilizer chain; enumerating its elements is impractical. The polyhedron
rotation groups use other generators than the reference, so element
generation order can differ.

`tensor_can.canonicalize` is the reference Butler-Portugal algorithm. The
`galois` module enumerates the transitive subgroups of `S1` through `S6`, and
`group_numbers` provides the order criteria (`is_nilpotent_number`,
`is_abelian_number`, `is_cyclic_number`) and `groups_count`.
