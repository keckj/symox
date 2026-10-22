"""Named finite permutation groups (generators as in the reference module)."""
from .permutations import Permutation
from .perm_groups import PermutationGroup
from .generators import _degree
from ..utilities.iterables import uniq
_af_new = Permutation._af_new


def DirectProduct(*groups):
    degrees = []
    gens_count = []
    total_degree = 0
    total_gens = 0
    for group in groups:
        current_deg = group.degree
        current_num_gens = len(group.generators)
        degrees.append(current_deg)
        total_degree += current_deg
        gens_count.append(current_num_gens)
        total_gens += current_num_gens
    array_gens = []
    for i in range(total_gens):
        array_gens.append(list(range(total_degree)))
    current_gen = 0
    current_deg = 0
    for i in range(len(gens_count)):
        for j in range(current_gen, current_gen + gens_count[i]):
            gen = groups[i].generators[j - current_gen].array_form
            array_gens[j][current_deg:current_deg + degrees[i]] = [x + current_deg for x in gen]
        current_gen += gens_count[i]
        current_deg += degrees[i]
    perm_gens = list(uniq([_af_new(list(a)) for a in array_gens]))
    return PermutationGroup(perm_gens, dups=False)


def AbelianGroup(*cyclic_orders):
    groups = []
    degree = 0
    order = 1
    for size in cyclic_orders:
        degree += size
        order *= size
        groups.append(CyclicGroup(size))
    G = DirectProduct(*groups)
    G._is_abelian = True
    G._degree = degree
    G._order = order
    return G

def AlternatingGroup(n):
    n = _degree(n)
    if n == 0:
        return PermutationGroup(Permutation([0]))
    if n in (1, 2):
        return PermutationGroup([Permutation([0])])
    a = list(range(n))
    (a[0], a[1], a[2]) = (a[1], a[2], a[0])
    gen1 = a
    if n % 2:
        a = list(range(1, n))
        a.append(0)
        gen2 = a
    else:
        a = list(range(2, n))
        a.append(1)
        a.insert(0, 0)
        gen2 = a
    gens = [gen1, gen2]
    if gen1 == gen2:
        gens = gens[:1]
    G = PermutationGroup([_af_new(a) for a in gens], dups=False)
    set_alternating_group_properties(G, n, n)
    G._is_alt = True
    return G

def set_alternating_group_properties(G, n, degree):
    if n < 4:
        G._is_abelian = True
        G._is_nilpotent = True
    else:
        G._is_abelian = False
        G._is_nilpotent = False
    if n < 5:
        G._is_solvable = True
    else:
        G._is_solvable = False
    G._degree = degree
    G._is_transitive = True
    G._is_dihedral = False

def CyclicGroup(n):
    n = _degree(n, 1)
    a = list(range(1, n))
    a.append(0)
    gen = _af_new(a)
    G = PermutationGroup([gen])
    G._is_abelian = True
    G._is_nilpotent = True
    G._is_solvable = True
    G._degree = n
    G._is_transitive = True
    G._order = n
    G._is_dihedral = n == 2
    return G

def DihedralGroup(n):
    n = _degree(n, 1)
    if n == 1:
        return PermutationGroup([Permutation([1, 0])])
    if n == 2:
        return PermutationGroup([Permutation([1, 0, 3, 2]), Permutation([2, 3, 0, 1]), Permutation([3, 2, 1, 0])])
    a = list(range(1, n))
    a.append(0)
    gen1 = _af_new(a)
    a = list(range(n))
    a.reverse()
    gen2 = _af_new(a)
    G = PermutationGroup([gen1, gen2])
    if n & n - 1 == 0:
        G._is_nilpotent = True
    else:
        G._is_nilpotent = False
    G._is_dihedral = True
    G._is_abelian = False
    G._is_solvable = True
    G._degree = n
    G._is_transitive = True
    G._order = 2 * n
    return G

def SymmetricGroup(n):
    n = _degree(n)
    if n == 0:
        G = PermutationGroup(Permutation([]))
    elif n == 1:
        G = PermutationGroup([Permutation([0])])
    elif n == 2:
        G = PermutationGroup([Permutation([1, 0])])
    else:
        a = list(range(1, n))
        a.append(0)
        gen1 = _af_new(a)
        a = list(range(n))
        (a[0], a[1]) = (a[1], a[0])
        gen2 = _af_new(a)
        G = PermutationGroup([gen1, gen2])
    set_symmetric_group_properties(G, n, n)
    G._is_sym = True
    return G

def set_symmetric_group_properties(G, n, degree):
    if n < 3:
        G._is_abelian = True
        G._is_nilpotent = True
    else:
        G._is_abelian = False
        G._is_nilpotent = False
    if n < 5:
        G._is_solvable = True
    else:
        G._is_solvable = False
    G._degree = degree
    G._is_transitive = True
    G._is_dihedral = n in [2, 3]


class SymmetricPermutationGroup(PermutationGroup):
    def __init__(self,deg):
        super().__init__(SymmetricGroup(deg).generators)

    def order(self):
        from math import factorial
        return factorial(self.degree)

    def contains(self,g,strict=True):
        permutation=g
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


