"""Exact group criteria, representations and signed tensor orbit checks."""
import importlib
import itertools
import random
import pytest
import symox as S
import sympy as R


@pytest.mark.parametrize('n',range(1,101))
def test_group_order_criteria_and_counts(n):
    native=importlib.import_module('symox.combinatorics.group_numbers');reference=importlib.import_module('sympy.combinatorics.group_numbers')
    for name in ['is_nilpotent_number','is_cyclic_number','is_abelian_number']:assert getattr(native,name)(n)==getattr(reference,name)(n)
    try:expected=reference.groups_count(n)
    except ValueError:
        with pytest.raises(ValueError):native.groups_count(n)
    else:assert native.groups_count(n)==expected


def test_gray_codes_partitions_and_subset_helpers():
    from symox.combinatorics.graycode import bin_to_gray,gray_to_bin,graycode_subsets
    from symox.combinatorics.partitions import RGS_generalized,random_integer_partition
    from symox.combinatorics.subsets import ksubsets
    from sympy.combinatorics.partitions import RGS_generalized as reference
    for width in range(1,7):
        for i in range(1<<width):
            value=format(i,'0%db'%width);assert gray_to_bin(bin_to_gray(value))==value
    assert len(list(graycode_subsets('abcd')))==16
    assert list(ksubsets('abc',2))==[('a','b'),('a','c'),('b','c')]
    for size in range(8):assert RGS_generalized(size).tolist()==[[S.sympify(v) for v in row] for row in reference(size).tolist()]
    assert random_integer_partition(100,seed=[1,1,12,1,2,1,85,1])==[85,12,2,1]
    assert sum(random_integer_partition(123,seed=42))==123


def test_named_transitive_groups_and_rubik_representations():
    native=importlib.import_module('symox.combinatorics.galois');reference=importlib.import_module('sympy.combinatorics.galois')
    from symox.combinatorics.generators import rubik,rubik_cube_generators
    for enum in ['S5TransitiveSubgroups','S6TransitiveSubgroups']:
        for name in getattr(native,enum):
            group=name.get_perm_group();expected=getattr(reference,enum)[name.name].get_perm_group()
            assert group.order()==expected.order() and group.is_transitive()
            assert {p.signature() for p in group.generators}=={p.signature() for p in expected.generators}
    groups=native.find_transitive_subgroups_of_S6(native.S6TransitiveSubgroups.S4p)
    assert groups[native.S6TransitiveSubgroups.S4p].order()==24
    assert all(p.size==48 and p.order()==4 for p in rubik_cube_generators())
    assert len(rubik(3))==6 and all(p.size==54 and p.order()==4 for p in rubik(3))


def test_symbolic_permutation_application_and_alias_identity():
    from symox.combinatorics.permutations import Permutation,Perm,AppliedPermutation
    from symox.combinatorics.perm_groups import PermutationGroup,PermGroup
    x=S.Symbol('x');p=Permutation([1,2,0])
    assert Perm is Permutation and PermGroup is PermutationGroup
    application=AppliedPermutation(p,x)
    assert application.subs(x,1)==2
    assert (2*application).subs(x,2)==0
    assert p.apply(x).subs(x,0)==1


def test_tensor_canonicalization_against_reference_for_random_signed_products():
    native=importlib.import_module('symox.combinatorics.tensor_can');reference=importlib.import_module('sympy.combinatorics.tensor_can')
    from symox.combinatorics import Permutation
    from sympy.combinatorics import Permutation as RefPermutation
    rng=random.Random(42)
    for antisymmetric in [False,True]:
        base,gens=native.get_symmetric_group_sgs(2,antisymmetric);rb,rg=reference.get_symmetric_group_sgs(2,antisymmetric)
        for metric in [0,1,None]:
            for exchange in [0,1,None]:
                for _ in range(5):
                    values=list(range(4));rng.shuffle(values);values.extend([4,5])
                    g=Permutation(values);expected=reference.canonicalize(RefPermutation(values),[2,3],metric,(rb,rg,2,exchange))
                    assert native.canonicalize(g,[2,3],metric,(base,gens,2,exchange))==expected
    base,gens=native.get_symmetric_group_sgs(2,True)
    assert native.canonicalize(Permutation(range(4)),[0,1],0,(base,gens,1,0))==0


def test_tensor_helpers_and_graph_certificate():
    native=importlib.import_module('symox.combinatorics.tensor_can');reference=importlib.import_module('sympy.combinatorics.tensor_can')
    from symox.combinatorics import Permutation,PermutationGroup
    from symox.combinatorics.testutil import graph_certificate
    for signed in [False,True]:
        a=[[1,0,2,3]] if signed else [[1,0]];b=[[0,1,3,2]] if signed else [[1,0]]
        assert native.perm_af_direct_product(a,b,signed)==reference.perm_af_direct_product(a,b,signed)
    assert native.dummy_sgs(list(range(2,8)),0,8)==reference.dummy_sgs(list(range(2,8)),0,8)
    base,gens=native.get_symmetric_group_sgs(3,True);transversals=native.get_transversals(base,gens)
    assert set(transversals[0])=={0,1,2} and set(transversals[1])=={1,2}
    newbase,newgens=native.get_minimal_bsgs(base,gens)
    assert newbase==base and PermutationGroup(newgens).order()==6
    a={0:[1],1:[0,2],2:[1]};b={3:[7,5],5:[3],7:[3]}
    assert graph_certificate(a)==graph_certificate(b)
