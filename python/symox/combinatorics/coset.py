"""Finite permutation cosets."""
from .permutations import Permutation
from .perm_groups import PermutationGroup


class Coset:
    def __init__(self,g,H,G=None,dir='+'):
        if not isinstance(g,Permutation) or not isinstance(H,PermutationGroup):raise TypeError('expected permutation and subgroup')
        if dir not in ('+','-'):raise ValueError('coset direction must be + or -')
        if G is not None and (not G.contains(g,strict=False) or any(not G.contains(h,strict=False) for h in H.generators)):
            raise ValueError('coset or subgroup is outside the ambient group')
        self.g,self.H,self.G,self.dir=g,H,G,dir
    @property
    def is_left_coset(self):return self.dir=='-'
    @property
    def is_right_coset(self):return self.dir=='+'
    def as_list(self):
        return [self.g*h if self.is_left_coset else h*self.g for h in self.H._enumerate()]
    def __contains__(self,g):return g in self.as_list()
    def __eq__(self,other):return isinstance(other,Coset) and set(self.as_list())==set(other.as_list())
