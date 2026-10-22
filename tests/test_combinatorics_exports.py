"""Finite combinatorial bijections and independently checked group actions."""
import itertools
import pickle
import pytest
import symox as s
import symox.combinatorics as c
import sympy.combinatorics as sc


def arrays(group):return {tuple(p.array_form) for p in group.elements}


@pytest.mark.parametrize('name',['SymmetricGroup','AlternatingGroup','CyclicGroup','DihedralGroup'])
@pytest.mark.parametrize('n',range(1,7))
def test_named_groups(name,n):
    a=getattr(c,name)(n);b=getattr(sc,name)(n)
    assert a.degree==b.degree
    assert a.order()==b.order()
    assert arrays(a)==arrays(b)


@pytest.mark.parametrize('name',['symmetric','alternating','cyclic','dihedral'])
@pytest.mark.parametrize('n',range(1,6))
def test_group_enumerators(name,n):
    assert [p.array_form for p in getattr(c,name)(n)]==[p.array_form for p in getattr(sc,name)(n)]


@pytest.mark.parametrize('orders',[(2,3),(2,2,3),(1,),(4,2)])
def test_direct_products(orders):
    a=c.AbelianGroup(*orders);b=sc.AbelianGroup(*orders)
    assert a.order()==b.order() and a.is_abelian
    assert arrays(a)==arrays(b)
    assert arrays(c.DirectProduct(*(c.CyclicGroup(n) for n in orders)))==arrays(a)


@pytest.mark.parametrize('n',range(1,9))
def test_gray_codes_and_subsets(n):
    a=c.GrayCode(n);b=sc.GrayCode(n)
    codes=list(a.generate_gray());expected=list(b.generate_gray())
    assert codes==expected and len(set(codes))==1<<n
    assert all(sum(x!=y for x,y in zip(u,v))==1 for u,v in zip(codes,codes[1:]))
    for rank,bits in enumerate(codes):
        assert c.GrayCode(n,start=bits).rank==rank
        assert c.GrayCode.unrank(n,rank)==bits
        subset=c.Subset.unrank_binary(rank,list(range(n)))
        oracle=sc.Subset(subset.subset,subset.superset)
        assert subset.rank_binary==rank
        assert subset.rank_gray==oracle.rank_gray
        assert subset.rank_lexicographic==oracle.rank_lexicographic
        assert subset.next_lexicographic().subset==oracle.next_lexicographic().subset
        assert subset.prev_lexicographic().subset==oracle.prev_lexicographic().subset
        assert subset.next_gray().subset==oracle.next_gray().subset
    assert c.GrayCode(n,rank=(1<<n)-1).next().rank==0


@pytest.mark.parametrize('n',range(2,7))
def test_prufer_bijection(n):
    limit=min(n**(n-2),200)
    for rank in range(limit):
        a=c.Prufer.unrank(rank,n);b=sc.Prufer.unrank(rank,n)
        assert a.prufer_repr==b.prufer_repr
        assert a.tree_repr==b.tree_repr
        assert a.rank==rank
        assert c.Prufer(a.tree_repr,n).rank==rank
        assert a.next().prev().rank==rank
    assert c.Prufer.edges([1,2,3],[2,4,5])==([[0,1],[1,2],[1,3],[3,4]],5)
    for edges in [[(0,1),(1,2),(0,2)],[(0,1),(0,1)]]:
        with pytest.raises(ValueError):c.Prufer(edges,4)


@pytest.mark.parametrize('n',range(1,8))
def test_restricted_growth_partitions(n):
    assert c.RGS_enum(n)==sc.RGS_enum(n)
    for rank in range(c.RGS_enum(n)):
        values=c.RGS_unrank(rank,n)
        assert values==sc.RGS_unrank(rank,n)
        assert c.RGS_rank(values)==rank
        p=c.Partition.from_rgs(values,list(range(n)))
        assert p.rank==rank and p.RGS==tuple(values)
        assert (p+1-1)==p
    with pytest.raises(ValueError):c.RGS_rank([0,2])


@pytest.mark.parametrize('n',range(1,16))
def test_integer_partition_neighbors(n):
    a=c.IntegerPartition([1]*n);b=sc.IntegerPartition([1]*n)
    for _ in range(int(__import__('sympy').partition(n))):
        assert a.partition==b.partition
        assert a.prev_lex().partition==b.prev_lex().partition
        assert a.as_dict()==b.as_dict() and a.conjugate==b.conjugate
        assert a.as_ferrers()==b.as_ferrers()
        a,b=a.next_lex(),b.next_lex()
    assert a.partition==[1]*n


@pytest.mark.parametrize('name',['tetrahedron','cube','octahedron','dodecahedron','icosahedron'])
def test_regular_polyhedra(name):
    a=getattr(c,name).copy();b=getattr(sc,name)
    assert a.faces==set(map(tuple,b.faces)) and a.edges==set(map(tuple,b.edges))
    assert arrays(a.pgroup)==arrays(b.pgroup)
    for g in a.pgroup.generators:
        original=a.corners;a.rotate(g)
        assert a.corners==tuple(original[g(i)] for i in range(a.size))
        a.rotate(~g);assert a.corners==original
    a.reset();assert a.array_form==list(range(a.size))


@pytest.mark.parametrize('n',[2,3,4,5])
def test_rubik_slice_actions(n):
    a=c.RubikGroup(n);b=sc.RubikGroup(n)
    assert a.degree==6*n*n
    assert [g.array_form for g in a.generators]==[g.array_form for g in b.generators]
    assert all(g.order()==4 for g in a.generators)


def test_cycles_free_words_and_cosets():
    a=c.Cycle(1,2,3)(2,4);b=sc.Cycle(1,2,3)(2,4)
    assert a.list()==b.list() and a.size==b.size
    assert c.Permutation(a).array_form==sc.Permutation(b).array_form
    F,x,y=c.free_group('x,y');G,u,v=sc.free_group('x,y')
    word=x*y*x**-1*y**2
    assert word.to_sympy()==u*v*u**-1*v**2
    assert (x*y*y**-1*x**-1)==F.identity
    assert word.inverse()*word==F.identity
    assert word.exponent_sum(x)==0 and word.exponent_sum(y)==3
    assert word.subword(1,3).to_sympy()==(u*v*u**-1*v**2).subword(1,3)
    H=c.CyclicGroup(3);g=c.Permutation([1,0,2])
    for direction in ['+','-']:
        a=c.Coset(g,H,dir=direction);b=sc.Coset(g.to_sympy(),sc.CyclicGroup(3),dir=direction)
        assert a.is_left_coset==b.is_left_coset
        assert {tuple(p.array_form) for p in a.as_list()}=={tuple(p.array_form) for p in b.as_list()}
    a=c.SymmetricPermutationGroup(12)
    assert a.order()==479001600 and a.contains(c.Permutation(range(12)))
    assert c.Permutation(range(12)) in a


@pytest.mark.parametrize('group',[c.SymmetricGroup(3),c.SymmetricGroup(4),c.DihedralGroup(4),c.AbelianGroup(4,2)])
def test_polycyclic_collection_reconstructs_every_element(group):
    pc=group.polycyclic_group();collector=pc.collector
    direct=c.Collector(pc.pcgs,pc.pc_series,pc.relative_order)
    assert c.PolycyclicGroup(pc.pcgs,pc.pc_series,pc.relative_order).is_solvable
    assert pc.length()==len(pc.relative_order)
    # As in the reference, the series refines the derived series and its
    # relative orders need not be prime (AbelianGroup(4, 2) gives [2, 4]).
    assert pc.is_prime_order()==all(s.isprime(n) for n in pc.relative_order)
    for element in group.elements:
        exponents=collector.exponent_vector(element)
        result=group.identity;word=collector.free_group.identity
        for g,x,n in zip(pc.pcgs,collector.free_group.generators,exponents):
            result*=g**n;word*=x**n
        assert result==element
        assert direct.exponent_vector(element)==exponents
        assert collector.collected_word(word)==word
        assert collector.collected_word(word**-1*word)==collector.free_group.identity
        assert collector.constructive_membership_test(pc.pcgs,element) is not False
    for relation,value in collector.pc_presentation.items():
        assert collector._permutation(relation)==collector._permutation(value)
    subgroup=c.CyclicGroup(3)
    with pytest.raises(ValueError):c.Collector(subgroup.generators,[],[2])


def test_trivial_polycyclic_group_and_prufer_optional_count():
    pc=c.SymmetricGroup(1).polycyclic_group()
    assert pc.relative_order==[] and pc.length()==0
    assert pc.collector.exponent_vector(c.Permutation([0]))==[]
    assert c.Prufer([0,2],4).tree_repr==sc.Prufer([0,2],4).tree_repr
