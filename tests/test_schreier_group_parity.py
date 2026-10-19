"""Exact group cardinalities/membership, checked against complete small closures."""
from collections import deque
import itertools
import math
import random

import pytest
from sympy.combinatorics import Permutation as RP,PermutationGroup as RG
from symox.combinatorics import (Permutation as P,PermutationGroup as G,Cycle,
    SymmetricGroup,DihedralGroup,AbelianGroup)


def closure(generators,degree):
    identity=tuple(range(degree));result={identity};pending=deque([identity])
    while pending:
        element=pending.popleft()
        for generator in generators:
            candidate=tuple(generator[image] for image in element)
            if candidate not in result:result.add(candidate);pending.append(candidate)
    return result


def arrays(group):return {tuple(p.array_form) for p in group.generate_schreier_sims()}


@pytest.mark.parametrize('degree',range(5))
def test_every_small_generator_pair_against_complete_independent_closure(degree):
    permutations=list(itertools.permutations(range(degree)))
    for first,second in itertools.product(permutations,repeat=2):
        expected=closure([first,second],degree)
        group=G(P(first),P(second))
        assert group.order()==len(expected)
        assert arrays(group)==expected
        for candidate in permutations:
            assert group.contains(P(candidate))==(candidate in expected)
        chain=group._chain()
        assert math.prod(len(level) for level in chain.transversals)==len(expected)
        prefix=[]
        for level,basepoint in enumerate(group.base):
            stabilizer={p for p in expected if all(p[point]==point for point in prefix)}
            assert arrays(group.basic_stabilizers[level])==stabilizer
            for image,representative in group.basic_transversals[level].items():
                assert representative(basepoint)==image
                assert all(representative(point)==point for point in prefix)
            prefix.append(basepoint)
        assert {p for p in expected if all(p[point]==point for point in prefix)}=={tuple(range(degree))}


def test_random_larger_groups_against_reference_order_and_membership():
    rng=random.Random(1729)
    for degree in range(5,10):
        for _ in range(12):
            generators=[]
            for _ in range(3):
                array=list(range(degree));rng.shuffle(array);generators.append(array)
            native,reference=G([P(g) for g in generators]),RG([RP(g) for g in generators])
            assert native.order()==reference.order()
            for _ in range(10):
                candidate=list(range(degree));rng.shuffle(candidate)
                assert native.contains(P(candidate))==reference.contains(RP(candidate))
            for action in ['tuples','sets','union']:
                for points in [(0,degree-1),(0,0)]:
                    # The reference union path appends to its input list.
                    assert native.orbit(list(points),action=action)==reference.orbit(list(points),action=action)
            assert native.orbits()==reference.orbits()
            assert native.is_transitive()==reference.is_transitive()
            assert native.is_transitive(strict=False)==reference.is_transitive(strict=False)
            assert native.stabilizer(0).order()==reference.stabilizer(0).order()
            base,strong=native.schreier_sims_incremental(base=[degree-1])
            assert base[0]==degree-1
            assert RG([RP(g.array_form) for g in strong]).order()==reference.order()


def test_large_orders_and_membership_do_not_enumerate_group_elements():
    degree=20
    group=G(SymmetricGroup(degree).generators)
    assert group.order()==math.factorial(degree)
    candidate=list(reversed(range(degree)))
    assert group.contains(P(candidate))
    assert group.stabilizer(0).order()==math.factorial(degree-1)
    assert group.orbit(0)==set(range(degree))
    first=list(itertools.islice(group.generate_schreier_sims(),10))
    assert len(set(first))==10 and all(group.contains(p) for p in first)
    assert group._elements is None
    F=P(2,19,21,8)(3,17,20,10)(4,6,7,5)
    R=P(1,5,21,14)(3,7,23,12)(8,10,11,9)
    D=P(6,18,14,10)(7,19,15,11)(20,22,23,21)
    rubik=G(F,R,D)
    assert rubik.order()==3674160
    assert rubik.contains(P(F*R*D,size=24))
    assert rubik._elements is None


def test_cyclic_membership_uses_compatible_shifts_without_building_transversals():
    generator=P(0,1,2,3)(4,5,6,7,8,9)
    group=G(generator)
    assert group.order()==12
    assert group.contains(P(0,1,2,3)(4,7)(5,8)(6,9))
    assert not group.contains(P(0,1,2,3)(4,6,8)(5,7,9))
    assert group._chain_cache is None
    degree=10000;cycle=P(list(range(1,degree))+[0]);large=G(cycle)
    assert large.order()==degree and large.contains(cycle**321)
    assert large.orbit(list(range(degree)),action='union')==set(range(degree))
    assert large._chain_cache is None and large._elements is None


def test_constructor_cycles_duplicates_empty_and_generator_words():
    assert G().degree==0 and G().order()==1
    native=G(Cycle(2,3)(4,5),Cycle(1,2,3,4,5))
    reference=RG(RP(2,3)(4,5),RP(1,2,3,4,5))
    assert native.order()==reference.order()
    for strict in [True,False]:
        for p in [P(1,2),P(0,1,2),P([0]),P(1,2,3,4,5)]:
            assert native.contains(p,strict=strict)==reference.contains(RP(p.array_form),strict=strict)
    native=G(P(2,1),P(1,0));reference=RG(RP(2,1),RP(1,0))
    for indices in [[],[0,1],[1,0],[0,0,1,1,0]]:
        assert native.make_perm(indices).array_form==reference.make_perm(indices).array_form
    for seed in [0,42,1729]:
        assert native.make_perm(20,seed=seed).array_form==reference.make_perm(20,seed=seed).array_form
    p=P([1,0]);assert len(G(p,p))==1 and len(G(p,p,dups=False))==2
    assert G(P([0]),P([1,0]),degree=3).degree==3


def test_derived_and_lower_central_series_without_element_enumeration():
    groups=[SymmetricGroup(n) for n in range(2,9)]+[DihedralGroup(4),DihedralGroup(6),AbelianGroup(2,3)]
    for original in groups:
        native=G(original.generators);reference=RG([RP(g.array_form) for g in native.generators])
        assert [g.order() for g in native.derived_series()]==[g.order() for g in reference.derived_series()]
        assert [g.order() for g in native.lower_central_series()]==[g.order() for g in reference.lower_central_series()]
        assert native.is_solvable==reference.is_solvable
        assert native.is_nilpotent==reference.is_nilpotent
        assert native._elements is None
