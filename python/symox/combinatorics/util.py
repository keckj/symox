"""Finite stabilizer and transversal utilities for native permutation groups."""
from .permutations import Permutation
from .perm_groups import PermutationGroup


def _base_ordering(base,degree):
    order=list(base)+[i for i in range(degree) if i not in base]
    return [order.index(i) for i in range(degree)]

def _distribute_gens_by_base(base,gens):
    return [[p for p in gens if all(p(i)==i for i in base[:level])] for level in range(len(base))]

def _orbits_transversals_from_bsgs(base,strong_gens_distr,transversals_only=False):
    transversals=[];orbits=[]
    for point,gens in zip(base,strong_gens_distr):
        group=PermutationGroup(gens or [Permutation(range(max(base,default=-1)+1))]);representatives={}
        for p in group._enumerate():representatives.setdefault(p(point),p)
        transversals.append(representatives);orbits.append(list(representatives))
    return transversals if transversals_only else (orbits,transversals)
