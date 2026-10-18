"""Exact finite groups with compact permutation generators."""
from .permutations import Permutation
from .perm_groups import PermutationGroup
from .generators import _degree


def _cycle(n,indices):
    a=list(range(n));indices=list(indices)
    for i,j in zip(indices,indices[1:]+indices[:1]):
        a[i]=j
    return Permutation(a)


def SymmetricGroup(n):
    n=_degree(n)
    gens=[_cycle(n,range(n))]
    if n>1:
        gens.append(_cycle(n,[0,1]))
    return PermutationGroup(gens)


def AlternatingGroup(n):
    n=_degree(n)
    if n<3:
        return PermutationGroup(Permutation([0]))
    return PermutationGroup([_cycle(n,[0,1,i]) for i in range(2,n)] or [Permutation(range(n))])


def CyclicGroup(n):
    n=_degree(n,1)
    return PermutationGroup(_cycle(n,range(n)))


def DihedralGroup(n):
    n=_degree(n,1)
    if n==1:
        return SymmetricGroup(2)
    if n==2:
        return PermutationGroup(Permutation([1,0,3,2]),Permutation([2,3,0,1]))
    return PermutationGroup(_cycle(n,range(n)),Permutation(list(range(n))[::-1]))


def DirectProduct(*groups):
    n=sum(g.degree for g in groups)
    generators=[];offset=0
    for group in groups:
        for p in group.generators:
            images=list(range(n))
            for i in range(group.degree):
                images[offset+i]=offset+p(i)
            generators.append(Permutation(images))
        offset+=group.degree
    return PermutationGroup(generators or [Permutation([])])


def AbelianGroup(*cyclic_orders):
    return DirectProduct(*(CyclicGroup(n) for n in cyclic_orders))


class SymmetricPermutationGroup(PermutationGroup):
    def __init__(self,deg):
        super().__init__(SymmetricGroup(deg).generators)

    def order(self):
        from math import factorial
        return factorial(self.degree)

    def contains(self,permutation,strict=True):
        return isinstance(permutation,Permutation) and (permutation.size==self.degree or not strict and permutation.size<self.degree)

    __contains__ = contains


def RubikGroup(n):
    """Slice rotations on 6*n*n facelets, with fixed cube orientation."""
    n=_degree(n,2)
    if n>32:raise NotImplementedError('Rubik cubes are limited to side length 32')
    m=n-1;stickers=[]
    for f in range(6):
        for r in range(n):
            for c in range(n):
                a,b=-m+2*c,m-2*r
                position=((a,m,-b),(a,b,m),(m,b,-a),(-a,b,-m),(-m,b,a),(a,-m,b))[f]
                normal=((0,1,0),(0,0,1),(1,0,0),(0,0,-1),(-1,0,0),(0,-1,0))[f]
                stickers.append((position,normal))
    lookup={s:i for i,s in enumerate(stickers)};gens=[]
    def turn(v,axis):
        x,y,z=v
        if axis==2:return y,-x,z
        if axis==0:return x,z,-y
        return z,y,-x
    for axis,sign in ((2,1),(0,1),(1,-1)):
        for layer in range(n-1):
            plane=sign*(m-2*layer);images=list(range(len(stickers)))
            for i,(position,normal) in enumerate(stickers):
                if position[axis]==plane:
                    dest=lookup[(turn(position,axis),turn(normal,axis))]
                    images[dest]=i
            gens.append(Permutation(images))
    return PermutationGroup(gens)
