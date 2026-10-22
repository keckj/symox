"""Exhaustive finite comparisons, exact rank arithmetic and large input bounds."""
import itertools
import math
import random

import pytest
from sympy.combinatorics import Permutation as Reference
from symox.combinatorics import Permutation


def images(value):return None if value is None else value.array_form


@pytest.mark.parametrize('degree',range(7))
def test_all_small_permutations_ranks_successors_and_discrete_invariants(degree):
    for array in itertools.permutations(range(degree)):
        native,reference=Permutation(array),Reference(array)
        for name in ['ascents','descents','index','inversions','inversion_vector',
                     'runs','support','length','parity','transpositions']:
            assert getattr(native,name)()==getattr(reference,name)(),(array,name)
        for name in ['cardinality','cycles','cycle_structure','full_cyclic_form',
                     'is_Empty','is_Singleton','is_Permutation','is_identity','is_odd']:
            assert getattr(native,name)==getattr(reference,name),(array,name)
        assert list(native)==list(reference)
        assert [i^native for i in range(degree)]==[i^reference for i in range(degree)]
        for shift in [-1,0,1]:
            assert (native+shift).array_form==(reference+shift).array_form
            assert (native-shift).array_form==(reference-shift).array_form
        assert native.copy()==native and native.copy() is not native
        assert native.atoms()==reference.atoms()==set(array)
        for name in ['min','max']:
            try:expected=getattr(reference,name)()
            except ValueError:
                with pytest.raises(ValueError):getattr(native,name)()
            else:assert getattr(native,name)()==expected
        for rank,unrank,successor in [('rank','unrank_lex','next_lex'),
            ('rank_nonlex','unrank_nonlex','next_nonlex'),
            ('rank_trotterjohnson','unrank_trotterjohnson','next_trotterjohnson')]:
            value=getattr(native,rank)()
            assert value==getattr(reference,rank)(),(array,rank)
            assert getattr(Permutation,unrank)(degree,value).array_form==list(array)
            if degree>=2:
                assert images(getattr(native,successor)())==images(getattr(reference,successor)()),(array,successor)
        if degree:
            assert Permutation.from_inversion_vector(native.inversion_vector())==native


@pytest.mark.parametrize('degree',range(7))
def test_small_group_operations_and_distances(degree):
    other=list(reversed(range(degree)))
    q,Q=Permutation(other),Reference(other)
    for array in itertools.permutations(range(degree)):
        p,P=Permutation(array),Reference(array)
        for name in ['commutator','mul_inv']:
            assert getattr(p,name)(q).array_form==getattr(P,name)(Q).array_form
        assert p.commutes_with(q)==P.commutes_with(Q)
        assert (p^q).array_form==(P^Q).array_form
        assert p.get_positional_distance(q)==P.get_positional_distance(Q)
        assert p.get_precedence_distance(q)==P.get_precedence_distance(Q)
        assert p.get_adjacency_distance(q)==P.get_adjacency_distance(Q)
        assert Permutation.rmul(p,q,p).array_form==Reference.rmul(P,Q,P).array_form
        assert Permutation.rmul_with_af(p,q,p).array_form==Reference.rmul_with_af(P,Q,P).array_form


@pytest.mark.parametrize('degree',range(9))
def test_matrix_representations_and_nonlex_inverse_side_effect(degree):
    array=list(range(degree));random.Random(1729).shuffle(array)
    p,P=Permutation(array),Reference(array)
    assert p.get_adjacency_matrix().to_sympy()==P.get_adjacency_matrix()
    assert p.get_precedence_matrix().to_sympy()==P.get_precedence_matrix()
    inv_perm,reference_inverse=(~p).array_form,(~P).array_form
    assert p.rank_nonlex(inv_perm=inv_perm)==P.rank_nonlex(inv_perm=reference_inverse)
    assert inv_perm==reference_inverse


def test_josephus_sequence_construction_and_degree_changes():
    for degree in range(12):
        for step in [-1,0,1,2,8,32]:
            for survivors in [0,1,2,6,100]:
                assert Permutation.josephus(step,degree,survivors).array_form==Reference.josephus(step,degree,survivors).array_form
    for sequence in ['caba',[3,1,4,1,5],[]]:
        for key in [None,lambda value:-ord(value) if isinstance(value,str) else -value]:
            assert Permutation.from_sequence(i=sequence,key=key).array_form==Reference.from_sequence(i=sequence,key=key).array_form
    for array in [[],[0],[1,0,2,3],[2,1,0,4,3]]:
        p,P=Permutation(array),Reference(array)
        assert p.list()==P.list()
        for degree in range(9):
            assert p.list(size=degree)==P.list(size=degree)
            try:expected=P.resize(degree)
            except ValueError:
                with pytest.raises(ValueError):p.resize(n=degree)
            else:assert p.resize(n=degree).array_form==expected.array_form


def test_large_ranks_inversions_and_powers_stay_exact():
    degree=2000;reverse=Permutation(range(degree-1,-1,-1))
    last=math.factorial(degree)-1
    assert reverse.rank()==last
    assert Permutation.unrank_lex(size=degree,rank=last)==reverse
    assert Permutation.unrank_nonlex(n=degree,r=reverse.rank_nonlex())==reverse
    assert Permutation.unrank_trotterjohnson(size=degree,rank=reverse.rank_trotterjohnson())==reverse
    assert Permutation(range(100000-1,-1,-1)).inversions()==100000*99999//2
    exponent=10**2000+1
    assert reverse**exponent==reverse
    assert reverse**(-exponent)==reverse


def test_singleton_orderings_terminate_and_out_of_range_ranks_remain_permutations():
    # There is exactly one permutation of degree zero/one. SymPy 1.14's empty
    # next_lex raises IndexError and empty next_trotterjohnson returns itself.
    for degree in [0,1]:
        p=Permutation(range(degree))
        assert p.next_lex() is p.next_nonlex() is p.next_trotterjohnson() is None
    # Upstream Trotter-Johnson unranking can return duplicate images for ranks
    # outside [0,n!-1]. Native unranking extends the ordering periodically.
    for rank in [-10**30,-1,6,8,10**30]:
        assert sorted(Permutation.unrank_trotterjohnson(3,rank).array_form)==[0,1,2]


def test_degree_errors_and_invalid_inversion_vectors():
    p,q=Permutation([1,0]),Permutation([1,2,0])
    for method in ['commutator','get_adjacency_distance','get_precedence_distance','get_positional_distance']:
        with pytest.raises(ValueError):getattr(p,method)(q)
    with pytest.raises(ValueError):p^q
    for entries in [[2],[0,2],[-3]]:
        with pytest.raises(ValueError):Permutation.from_inversion_vector(entries)
