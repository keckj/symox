"""SLP certificates must reconstruct strong generators and certify stabilizers."""
import itertools
import math
import random
import subprocess
import sys

import pytest
from sympy.combinatorics import Permutation as RP,PermutationGroup as RG
from sympy.combinatorics.testutil import _verify_bsgs
from symox.combinatorics import Permutation as P,PermutationGroup as G,SymmetricGroup


def verify(generators,degree,base=()):
    native=G([P(g) for g in generators] or [P(range(degree))])
    reference=RG([RP(g) for g in generators] or [RP(range(degree))])
    supplied_base=list(base);supplied_gens=[P(g) for g in generators]
    saved_base=supplied_base[:];saved_gens=supplied_gens[:]
    result,strong,witnesses=native.schreier_sims_incremental(
        base=supplied_base,gens=supplied_gens,slp_dict=True)
    assert supplied_base==saved_base and supplied_gens==saved_gens
    assert result[:len(base)]==list(base)
    assert set(witnesses)==set(strong)
    earlier=set()
    original=set(supplied_gens)
    for generator in strong:
        word=witnesses[generator]
        assert P.rmul(*word)==generator
        if generator not in original:
            assert all(factor in earlier or ~factor in earlier for factor in word)
        earlier.add(generator)
    if strong and result:
        assert _verify_bsgs(reference,result,[RP(g.array_form) for g in strong])
    else:assert reference.order()==1
    plain_base,plain_strong=native.schreier_sims_incremental(base=list(base),gens=supplied_gens)
    assert (plain_base,plain_strong)==(result,strong)
    return result,strong,witnesses


@pytest.mark.parametrize('degree',range(5))
def test_all_small_generator_pairs_have_valid_witnesses_and_bsgs(degree):
    permutations=list(itertools.permutations(range(degree)))
    for generators in itertools.product(permutations,repeat=2):
        verify(generators,degree)


def test_larger_random_groups_and_redundant_base_points():
    rng=random.Random(1729)
    for degree in range(5,10):
        for _ in range(10):
            generators=[]
            for _ in range(3):
                array=list(range(degree));rng.shuffle(array);generators.append(array)
            verify(generators,degree)
            verify(generators,degree,[0,0,1])


@pytest.mark.parametrize('degree',[0,1,4])
def test_identity_empty_and_duplicate_generators(degree):
    identity=list(range(degree))
    for generators in [[],[identity],[identity,identity]]:
        base=[0,0] if degree else []
        result=verify(generators,degree,base)
        if len(generators)==1:
            reference=RG(RP(identity)).schreier_sims_incremental(
                base=base,gens=[RP(g) for g in generators],slp_dict=True)
            assert result[0]==reference[0]
            assert [g.array_form for g in result[1]]==[g.array_form for g in reference[1]]
        else:
            # Reference 1.14 indexes an empty filtered generator list. Native
            # can still return the exact, empty certificate of this group.
            with pytest.raises(IndexError):
                RG(RP(identity)).schreier_sims_incremental(
                    base=base,gens=[RP(g) for g in generators],slp_dict=True)
            assert result==(base,[],{})
    if degree>1:
        transposition=identity[:];transposition[:2]=[1,0]
        _,strong,_=verify([transposition,transposition,identity],degree)
        assert strong[:2]==[P(transposition)]*2


@pytest.mark.parametrize('degree',[10,20,30])
def test_large_symmetric_group_witnesses_do_not_enumerate_group(degree):
    group=G(SymmetricGroup(degree).generators)
    base,strong,witnesses=group.schreier_sims_incremental(slp_dict=True)
    assert group._elements is None and group._chain_cache is None
    assert all(P.rmul(*word)==g for g,word in witnesses.items())
    reference=RG([RP(g.array_form) for g in group.generators])
    assert reference.order()==math.factorial(degree)
    assert _verify_bsgs(reference,base,[RP(g.array_form) for g in strong])


def test_slp_implementation_never_imports_reference():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('reference import')
sys.meta_path.insert(0,Block())
from symox.combinatorics import Permutation,PermutationGroup,SymmetricGroup
group=PermutationGroup(SymmetricGroup(20).generators)
base,strong,words=group.schreier_sims_incremental(slp_dict=True)
assert all(Permutation.rmul(*word)==g for g,word in words.items())
assert group._elements is None
'''
    subprocess.run([sys.executable,'-c',source],check=True,timeout=20)
