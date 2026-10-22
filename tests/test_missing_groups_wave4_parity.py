"""Reference regressions for the ported permutation, finitely presented and
polycyclic group algorithms (wave 4, missing groups).

Expected values were computed with SymPy 1.14 and are compared through
``sstr``. Randomized algorithms are checked by invariants only.
"""
import pytest
import symox as S
import symox.combinatorics as C
import symox.combinatorics.free_groups as FG
import symox.combinatorics.fp_groups as FP
import symox.combinatorics.coset_table as CT
import symox.combinatorics.tensor_can as TC
import symox.combinatorics.testutil as TU
from symox.combinatorics.testutil import _verify_bsgs, _verify_centralizer, _verify_normal_closure

Permutation = C.Permutation
PG = C.PermutationGroup
Sym, Alt, Cyc, Dih, Ab = (C.SymmetricGroup, C.AlternatingGroup, C.CyclicGroup,
                          C.DihedralGroup, C.AbelianGroup)
F, a, b = FG.free_group('a b')


def _value(f):
    try:
        v = f()
    except Exception as e:
        return '%s: %s' % (type(e).__name__, e)
    return v if isinstance(v, str) else S.sstr(v)


CASES = [
    ('center_D6_gens', lambda: Dih(6).center().generators,
     '[(0 3)(1 4)(2 5)]'),
    ('center_D6_order', lambda: Dih(6).center().order(),
     '2'),
    ('center_S4_order', lambda: Sym(4).center().order(),
     '1'),
    ('centralizer_S6_C6', lambda: Sym(6).centralizer(Cyc(6)).generators,
     '[(0 1 2 3 4 5)]'),
    ('centralizer_S6_C6_order', lambda: Sym(6).centralizer(Cyc(6)).order(),
     '6'),
    ('centralizer_perm', lambda: Sym(5).centralizer(Permutation(0,1,2)).order(),
     'IndexError: list index out of range'),
    ('normal_closure_S5_C5', lambda: Sym(5).normal_closure(Cyc(5)).order(),
     '60'),
    ('normal_closure_S4_perm', lambda: Sym(4).normal_closure(Permutation(0,1)(2,3)).order(),
     '4'),
    ('normal_closure_is_normal', lambda: Sym(4).normal_closure(Permutation(0,1)(2,3)).is_normal(Sym(4)),
     'True'),
    ('derived_series_S4', lambda: [G.order() for G in Sym(4).derived_series()],
     '[24, 12, 4, 1]'),
    ('derived_series_S5', lambda: [G.order() for G in Sym(5).derived_series()],
     '[120, 60]'),
    ('lower_central_D8', lambda: [G.order() for G in Dih(8).lower_central_series()],
     '[16, 4, 2, 1]'),
    ('lower_central_A4', lambda: [G.order() for G in Alt(4).lower_central_series()],
     '[12, 4]'),
    ('composition_series_S4', lambda: [G.order() for G in Sym(4).composition_series()],
     '[24, 12, 4, 2, 1]'),
    ('composition_series_C12', lambda: [G.order() for G in Cyc(12).composition_series()],
     '[12, 6, 3, 1]'),
    ('composition_series_A5', lambda: [G.order() for G in Alt(5).composition_series()],
     'NotImplementedError: Group should be solvable'),
    ('is_primitive', lambda: [Dih(10).is_primitive(),Sym(5).is_primitive(),Cyc(7).is_primitive(),Cyc(6).is_primitive(),Alt(4).is_primitive(randomized=False)],
     '[False, True, True, False, True]'),
    ('minimal_block_D10', lambda: Dih(10).minimal_block([0,5]),
     '[0, 1, 2, 3, 4, 0, 1, 2, 3, 4]'),
    ('minimal_block_D10b', lambda: Dih(10).minimal_block([0,1]),
     '[0, 0, 0, 0, 0, 0, 0, 0, 0, 0]'),
    ('minimal_blocks_D6', lambda: Dih(6).minimal_blocks(),
     '[[0, 1, 0, 1, 0, 1], [0, 1, 2, 0, 1, 2]]'),
    ('minimal_blocks_S5', lambda: Sym(5).minimal_blocks(),
     '[[0, 0, 0, 0, 0]]'),
    ('minimal_blocks_intrans', lambda: PG(Permutation(1,2,5,6)).minimal_blocks(),
     'False'),
    ('sylow_S6', lambda: [Sym(6).sylow_subgroup(p).order() for p in (2,3,5)],
     '[16, 9, 5]'),
    ('sylow_D6', lambda: Dih(6).sylow_subgroup(2).order(),
     '4'),
    ('sylow_S10_2', lambda: Sym(10).sylow_subgroup(2).order(),
     '256'),
    ('sylow_A5_7', lambda: Alt(5).sylow_subgroup(7).order(),
     '1'),
    ('subgroup_search_even', lambda: Sym(6).subgroup_search(lambda x:x.is_even).order(),
     '360'),
    ('subgroup_search_fix', lambda: Sym(6).subgroup_search(lambda x:x(0)==0 and x(1)==1, tests={0:lambda c:c(0)==0,1:lambda c:c(1)==1}).order(),
     'KeyError: 4'),
    ('ssi_S5', lambda: Sym(5).schreier_sims_incremental(),
     '([0, 1, 3, 2], [(0 1 2 3 4), (4)(0 1), (1 4 3 2), (3 4), (2 4)])'),
    ('ssi_D6_base', lambda: Dih(6).schreier_sims_incremental(base=[2,1]),
     '([2, 1], [(0 1 2 3 4 5), (0 5)(1 4)(2 3), (5)(0 4)(1 3)])'),
    ('base_strong_S5', lambda: (Sym(5).base,Sym(5).strong_gens),
     '([0, 1, 3, 2], [(0 1 2 3 4), (4)(0 1), (1 4 3 2), (3 4), (2 4)])'),
    ('base_strong_A5', lambda: (Alt(5).base,Alt(5).strong_gens),
     '([0, 2, 1], [(4)(0 1 2), (0 1 2 3 4), (2 3 4), (1 3 4)])'),
    ('base_strong_D7', lambda: (Dih(7).base,Dih(7).strong_gens),
     '([0, 1], [(0 1 2 3 4 5 6), (0 6)(1 5)(2 4), (1 6)(2 5)(3 4)])'),
    ('basic_orbits_A5', lambda: Alt(5).basic_orbits,
     '[[0, 1, 2, 3, 4], [1, 2, 3, 4], [1, 3, 4]]'),
    ('coset_factor', lambda: Sym(5).coset_factor(Permutation(0,3,1,4)),
     '[(0 3 1 4 2), (4), (4), (2 4)]'),
    ('coset_factor_af', lambda: Sym(5).coset_factor(Permutation(0,3,1,4),True),
     '[3, 1, 3, 4]'),
    ('coset_factor_out', lambda: Alt(5).coset_factor(Permutation(0,1)),
     'ValueError: g should be the same size as permutations of G'),
    ('coset_rank', lambda: [Sym(5).coset_rank(Permutation(0,3,1,4)),Alt(5).coset_rank(Permutation(0,1))],
     'ValueError: g should be the same size as permutations of G'),
    ('coset_unrank', lambda: [Sym(4).coset_unrank(k) for k in (0,5,23,24)],
     '[(3), (0 1 3 2), (0 3 1), None]'),
    ('coset_unrank_af', lambda: Sym(4).coset_unrank(7,af=True),
     '[3, 2, 0, 1]'),
    ('coset_ranks_bij', lambda: sorted(Dih(5).coset_rank(g) for g in Dih(5).elements),
     '[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]'),
    ('abelian_invariants', lambda: [Sym(3).abelian_invariants(),Ab(2,4,6).abelian_invariants(),Cyc(12).abelian_invariants(),Dih(8).abelian_invariants()],
     '[[2], [2, 2, 3, 4], [3, 4], [2, 2]]'),
    ('baseswap', lambda: Sym(4).baseswap(Sym(4).base,Sym(4).strong_gens,0,randomized=False),
     '([1, 0, 2], [(0 1 2 3), (3)(0 1), (1 3 2), (2 3), (0 2 3)])'),
    ('commutator', lambda: Sym(5).commutator(Sym(5),Alt(5)).order(),
     '60'),
    ('equals', lambda: [PG(Permutation(0,1,2,3,4),Permutation(1,4)(2,3)).equals(Dih(5)),Sym(4).equals(Alt(4))],
     '[True, False]'),
    ('index', lambda: [Sym(5).index(Alt(5)),Sym(4).index(Dih(4))],
     '[2, 3]'),
    ('is_alt_sym', lambda: [Dih(10).is_alt_sym(),Sym(7).is_alt_sym(),Alt(8).is_alt_sym()],
     '[False, True, True]'),
    ('is_elementary', lambda: [Ab(2,2).is_elementary(2),Ab(2,4).is_elementary(2),Cyc(3).is_elementary(3)],
     '[True, False, True]'),
    ('is_normal', lambda: [Alt(5).is_normal(Sym(5)),Cyc(3).is_normal(Sym(3)),PG(Permutation(0,1)).is_normal(Sym(3))],
     '[True, True, False]'),
    ('is_subgroup', lambda: [Alt(4).is_subgroup(Sym(4)),Sym(3).is_subgroup(Sym(4)),Sym(3).is_subgroup(Sym(4),strict=False)],
     '[True, False, True]'),
    ('orbit_rep', lambda: Dih(6).orbit_rep(0,4),
     '(5)(0 4)(1 3)'),
    ('orbit_rep_none', lambda: PG(Permutation(0,1),Permutation(2,3)).orbit_rep(0,3),
     'False'),
    ('orbit_transversal', lambda: Dih(5).orbit_transversal(0,pairs=True),
     '[(0, (4)), (1, (0 1 2 3 4)), (4, (0 4)(1 3)), (2, (0 2 4 1 3)), (3, (4)(0 3)(1 2))]'),
    ('pointwise_stabilizer', lambda: Sym(7).pointwise_stabilizer([2,3,5]).order(),
     '24'),
    ('stabilizer', lambda: Dih(6).stabilizer(2).generators,
     '[(5)(0 4)(1 3)]'),
    ('schreier_vector', lambda: PG(Permutation(0,1,3,4),Permutation(0,4,1,2,3)).schreier_vector(0),
     '[-1, 0, 1, 0, 1]'),
    ('schreier_vector2', lambda: PG(Permutation(0,1)(5,6),Permutation(5,6,7)).schreier_vector(4),
     '[None, None, None, None, -1, None, None, None]'),
    ('predicates', lambda: [(G.is_abelian,G.is_nilpotent,G.is_solvable,G.is_perfect,G.is_cyclic,G.is_dihedral,G.is_transitive(),G.transitivity_degree,G.is_trivial) for G in (Sym(4),Alt(5),Dih(4),Cyc(6),Ab(2,2))],
     '[(False, False, True, False, False, False, True, 4, False), (False, False, False, True, False, False, True, 3, False), (False, True, True, False, False, True, True, 1, False), (True, True, True, False, True, False, True, 1, False), (True, True, True, False, False, True, False, 0, False)]'),
    ('conj_classes', lambda: [len(G.conjugacy_classes()) for G in (Sym(4),Alt(5),Dih(5))],
     '[5, 5, 4]'),
    ('max_div', lambda: [Sym(5).max_div,Dih(6).max_div],
     '[1, 3]'),
    ('elements_count', lambda: len(Sym(5).elements),
     '120'),
    ('presentation_S4', lambda: (lambda P:(P.order(),len(P.generators)))(Sym(4).presentation()),
     '(24, 2)'),
    ('strong_presentation_D4', lambda: (lambda P:(P.order(),len(P.generators),len(P.relators)))(Dih(4).strong_presentation()),
     '(8, 3, 4)'),
    ('strong_presentation_A4', lambda: (lambda P:(P.order(),len(P.generators),len(P.relators)))(Alt(4).strong_presentation()),
     '(12, 2, 4)'),
    ('polycyclic_S4', lambda: Sym(4).polycyclic_group().relative_order,
     '[2, 3, 2, 2]'),
    ('polycyclic_D6', lambda: Dih(6).polycyclic_group().relative_order,
     '[2, 2, 3]'),
    ('fp_order', lambda: [FP.FpGroup(F,[a**2,b**3,(a*b)**4]).order(),FP.FpGroup(F,[a**2,b**3,(a*b)**5]).order(),FP.FpGroup(F,[a**3,b**2,(a*b)**3]).order()],
     '[24, 60, 12]'),
    ('fp_index', lambda: FP.FpGroup(F,[a**2,b**3,(a*b)**4]).index([a]),
     '12'),
    ('fp_subgroup_order', lambda: FP.FpGroup(F,[a**2,b**3,(a*b)**4]).subgroup([a*b]).order(),
     '4'),
    ('fp_infinite', lambda: FP.FpGroup(F,[a*b*a**-1*b**-1]).order(),
     'oo'),
    ('fp_abelian_invariants', lambda: FP.FpGroup(F,[a**2,b**4,a*b*a**-1*b**-1]).abelian_invariants(),
     '[2, 4]'),
    ('low_index', lambda: len(FP.low_index_subgroups(FP.FpGroup(F,[a**2,b**3,(a*b)**4]),4)),
     '4'),
    ('reidemeister', lambda: (lambda R:(len(R[0]),len(R[1])))(FP.reidemeister_presentation(FP.FpGroup(F,[a**3,b**5,(a*b)**2]),[a])),
     '(1, 1)'),
    ('simplify_presentation', lambda: FP.simplify_presentation(FP.FpGroup(F,[a**2,b**2,a*b*a**-1*b**-1,(a*b)**2])).relators,
     '[a**2, b**2, a*b*a**-1*b**-1, a*b*a*b]'),
    ('fp_is_abelian', lambda: FP.FpGroup(F,[a**2,b**3,a*b*a**-1*b**-1]).is_abelian,
     'True'),
    ('fp_center', lambda: len(FP.FpGroup(F,[a**4,b**2,(a*b)**2]).center()),
     '1'),
    ('fp_derived_series', lambda: [G.order() for G in FP.FpGroup(F,[a**2,b**3,(a*b)**4]).derived_series()],
     "AttributeError: 'list' object has no attribute 'order'"),
    ('fp_permutation_order', lambda: FP.FpGroup(F,[a**2,b**3,(a*b)**4])._to_perm_group()[0].order(),
     '24'),
    ('coset_enum_r', lambda: CT.coset_enumeration_r(FP.FpGroup(F,[a**3,b**2,(a*b)**3]),[a]).table,
     '[[0, 0, 1, 1], [2, 3, 0, 0], [3, 1, 3, 3], [1, 2, 2, 2]]'),
    ('coset_enum_c2', lambda: CT.coset_enumeration_c(FP.FpGroup(F,[a**2,b**3,(a*b)**4]),[a*b]).table,
     '[[1, 1, 2, 1], [0, 0, 0, 2], [3, 3, 1, 0], [2, 2, 4, 5], [5, 5, 5, 3], [4, 4, 3, 4], [4, None, None, None]]'),
    ('canon', lambda: TC.canonicalize(Permutation([2,1,0,3,4,5]),[0,1,2,3],0,(TC.riemann_bsgs[0],TC.riemann_bsgs[1],1,0)),
     '[0, 2, 1, 3, 5, 4]'),
    ('canon2', lambda: TC.canonicalize(Permutation([0,2,1,3,4,5]),[0,1,2,3],0,(TC.riemann_bsgs[0],TC.riemann_bsgs[1],1,0)),
     '[0, 2, 1, 3, 4, 5]'),
    ('canon3', lambda: TC.canonicalize(Permutation(list(range(12))+[12,13])*Permutation(1,6)(3,9),list(range(12)),0,(TC.riemann_bsgs[0],TC.riemann_bsgs[1],3,0)),
     '0'),
    ('graph_cert', lambda: TU.graph_certificate({0:[1,2,3,5],1:[0,2,4],2:[0,1,3,4],3:[0,2,4],4:[1,2,3,5],5:[0,4]}),
     '[0, 2, 4, 6, 1, 8, 10, 12, 3, 14, 16, 18, 5, 9, 15, 7, 11, 17, 13, 19, 20, 21]'),
    ('min_bsgs', lambda: TC.get_minimal_bsgs(*TC.riemann_bsgs),
     '([0, 2], [(0 1)(4 5), (2 3)(4 5), (5)(0 2)(1 3)])'),
    ('rmul', lambda: Permutation.rmul(Permutation([1,2,0]),[0,2,1]),
     '(2)(0 1)'),
    ('partition_key_sort', lambda: [C.Partition(*p).members for p in sorted([[[1,2]],[[1],[2]],[[1,2],[3]],[[3,4]]],key=lambda p:C.Partition(*p).sort_key())],
     '[(1, 2), (1, 2), (3, 4), (1, 2, 3)]'),
    ('direct_product', lambda: C.DirectProduct(Cyc(2),Cyc(3),Dih(3)).order(),
     '36'),
    ('rubik', lambda: C.RubikGroup(2).order(),
     '3674160'),
    ('mul_group', lambda: (Cyc(2)*Cyc(3)).order(),
     '6'),
    ('generator_product', lambda: Sym(4).generator_product(Permutation(0,2,1,3),original=True),
     '[(0 1 2 3), (0 1 2 3), (0 3 2 1), (0 3 2 1), (3)(0 1), (0 1 2 3), (0 1 2 3)]'),
    ('make_perm_seq', lambda: Sym(4).make_perm(3,[0,1,0]),
     '(0 1 3 2)'),
    ('elementary_abelian', lambda: Ab(2,2,2).is_elementary(2),
     'True'),
    ('degree_order', lambda: [(G.degree,G.order()) for G in (Sym(0) if False else Sym(1),Cyc(1),Dih(1),Dih(2),Ab(1))],
     '[(1, 1), (1, 1), (2, 2), (4, 4), (1, 1)]'),
    ('is_symmetric_alt', lambda: [Sym(5).is_symmetric,Alt(5).is_alternating,PG(Permutation(0,1,2,3,4),Permutation(0,1)).is_symmetric],
     '[True, True, True]'),
]


@pytest.mark.parametrize('name,f,expected', CASES, ids=[c[0] for c in CASES])
def test_reference_value(name, f, expected):
    assert _value(f) == expected


@pytest.mark.parametrize('group', [Sym(5), Alt(6), Dih(7), Ab(2, 4), Cyc(7), PG(Permutation(0, 1, 2)(3, 4), Permutation(5, 6))])
def test_randomized_bsgs_invariants(group):
    order = group.order()
    # Monte Carlo: about 1% of runs are incomplete with consec_succ=8.
    base, gens = group.schreier_sims_random(consec_succ=60)
    assert _verify_bsgs(group, base, gens)
    assert PG(gens).order() == order
    base, gens = group.schreier_sims_incremental()
    assert _verify_bsgs(group, base, gens)
    g = group.random_stab(0)
    assert g(0) == 0 and group.contains(g)
    assert group.order() // group.stabilizer(0).order() == len(group.orbit(0))


def test_schreier_sims_random_single_point_base():
    # Divergence: the reference raises IndexError when the base has one point.
    base, gens = Cyc(7).schreier_sims_random()
    assert base == [0] and PG(gens).order() == 7


def test_randomized_subgroups_by_invariants():
    S6 = Sym(6)
    for p, n in ((2, 16), (3, 9), (5, 5)):
        P = S6.sylow_subgroup(p)
        assert P.order() == n and P.is_subgroup(S6)
    closure = Sym(5).normal_closure(Cyc(5))
    assert _verify_normal_closure(Sym(5), Cyc(5), closure)
    assert closure.is_normal(Sym(5)) and closure.equals(Alt(5))
    assert _verify_centralizer(Sym(6), Cyc(6))
    assert _verify_centralizer(Dih(6), Dih(6))
    pc = Sym(4).polycyclic_group()
    assert pc.relative_order == [2, 3, 2, 2] and pc.length() == 4


def test_coset_rank_unrank_round_trip():
    G = Dih(5)
    for g in G.elements:
        assert G.coset_unrank(G.coset_rank(g)) == g
        assert G.coset_factor(g) and G.coset_factor(g, True)


def test_collector_membership_reports_exponents():
    # Divergence: the reference compares the sifted permutation with 1 and
    # always reports False, symox returns the exponent vector of members.
    pc = Sym(4).polycyclic_group()
    collector = pc.collector
    for g in Sym(4).elements:
        assert collector.constructive_membership_test(pc.pcgs, g) is not False


def test_strong_presentation_and_fp_isomorphism():
    for G in (Dih(4), Alt(4), Sym(4)):
        P = G.strong_presentation()
        assert P.order() == G.order()
        T, _ = P._to_perm_group()
        assert T.order() == G.order()


def test_induced_pcgs_is_complete():
    # Divergence: the reference sifts with a float inverse, never sifts the
    # last depth and omits powers, so induced_pcgs([(0 1 2 3)]) is
    # [(0 1 2 3)] and the membership test rejects every element.
    pcgs = [Permutation(0, 1, 2, 3), Permutation(3)(0, 2, 1),
            Permutation(0, 1)(2, 3), Permutation(0, 3)(1, 2)]
    c = C.Collector(pcgs, [PG(pcgs[i:] or [Permutation(3)]) for i in range(5)], [2, 3, 2, 2])
    gens = [Permutation(0, 1, 2, 3)]
    ipcgs = c.induced_pcgs(gens)
    assert gens == [Permutation(0, 1, 2, 3)]
    assert ipcgs == [Permutation(0, 1, 2, 3), Permutation(0, 2)(1, 3)]
    assert c.constructive_membership_test(ipcgs, Permutation(0, 2)(1, 3)) == [0, 1]
    assert c.constructive_membership_test(ipcgs, Permutation(0, 3, 2, 1)) == [1, 1]
    assert c.constructive_membership_test(ipcgs, Permutation(3)(0, 1)) is False
    S4 = Sym(4)
    for gens in ([Permutation(0, 1)(2, 3), Permutation(0, 3)(1, 2)], [Permutation(3)(0, 2, 1)],
                 [Permutation(0, 1, 2, 3), Permutation(3)(0, 2, 1)]):
        H = PG(gens)
        ipcgs = c.induced_pcgs(gens)
        assert PG(ipcgs).order() == H.order()
        for g in S4.elements:
            assert (c.constructive_membership_test(ipcgs, g) is not False) == H.contains(g)


def test_induced_pcgs_non_invertible_exponent_is_unsupported():
    # Relative orders such as [2, 2, 4] need an inverse of 2 modulo 4: the
    # sift raises NotImplementedError instead of returning a wrong sequence.
    pc = Dih(8).polycyclic_group()
    c = pc.collector
    for g in Dih(8).elements:
        for h in Dih(8).elements:
            try:
                ipcgs = c.induced_pcgs([g, h])
            except NotImplementedError:
                continue
            H = PG([g, h])
            assert PG(ipcgs or [Dih(8).identity]).order() == H.order()
