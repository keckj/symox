"""Finite polyhedra and exact orientation-preserving vertex symmetries."""
import copy
from .permutations import Permutation
from .perm_groups import PermutationGroup
from ..utilities.iterables import minlex


def _face(face):
    face=tuple(face)
    return min(seq[i:]+seq[:i] for seq in (face,face[::-1]) for i in range(len(face)))


class Polyhedron:
    def __init__(self,corners,faces=(),pgroup=()):
        from .. import sympify
        self._original=tuple(map(sympify,corners));self._corners=self._original
        if len(set(self._original))!=len(self._original):raise ValueError('corners must be distinct')
        self.faces=frozenset(_face(f) for f in faces)
        if any(len(f)<3 or len(set(f))!=len(f) or any(not 0<=i<self.size for i in f) for f in self.faces):
            raise ValueError('invalid polyhedron face')
        self.pgroup=pgroup if isinstance(pgroup,PermutationGroup) else PermutationGroup(list(pgroup) or [Permutation(range(self.size))])
        if self.pgroup.degree!=self.size:raise ValueError('permutation degree differs from corner count')
    @property
    def corners(self):return self._corners
    vertices=corners
    @property
    def size(self):return len(self._corners)
    @property
    def array_form(self):return [self._original.index(v) for v in self._corners]
    @property
    def cyclic_form(self):return Permutation(self.array_form).cyclic_form
    @property
    def edges(self):return frozenset(tuple(sorted((f[i-1],f[i]))) for f in self.faces for i in range(len(f)))
    def rotate(self,perm):
        if isinstance(perm,int):perm=self.pgroup.generators[perm]
        if not isinstance(perm,Permutation) or perm.size!=self.size:raise ValueError('invalid vertex permutation')
        self._corners=tuple(self._corners[perm(i)] for i in range(self.size))
    def reset(self):self._corners=self._original
    def copy(self):
        # SymPy's Basic.copy rebuilds from args, which hold the original corners.
        new=copy.copy(self);new._corners=new._original;return new


def _rotations(n,faces):
    # Orient neighbouring faces with opposite traversal on the common edge.
    faces=[tuple(f) for f in faces];oriented=[faces.pop(0)]
    while faces:
        for i,f in enumerate(faces):
            directions={(a,b) for g in oriented for a,b in zip(g,g[1:]+g[:1])}
            edges=list(zip(f,f[1:]+f[:1]))
            if any((a,b) in directions or (b,a) in directions for a,b in edges):
                if any((a,b) in directions for a,b in edges):f=f[::-1]
                oriented.append(f);faces.pop(i);break
        else:raise ValueError('disconnected face complex')
    # A map of one oriented face extends uniquely across adjacent faces.
    edge_face={(a,b):f for f in oriented for a,b in zip(f,f[1:]+f[:1])}
    source=oriented[0];perms=[]
    for target in oriented:
        for shift in range(len(target)):
            image=target[shift:]+target[:shift];mapping=dict(zip(source,image));queue=[source];visited=set()
            valid=True
            for f in queue:
                if f in visited:continue
                visited.add(f)
                for a,b in zip(f,f[1:]+f[:1]):
                    neighbor=edge_face[(b,a)]
                    mapped=edge_face[(mapping[b],mapping[a])]
                    k=neighbor.index(b);l=mapped.index(mapping[b])
                    for j in range(len(neighbor)):
                        old,new=neighbor[(k+j)%len(neighbor)],mapped[(l+j)%len(mapped)]
                        if old in mapping and mapping[old]!=new:valid=False;break
                        mapping[old]=new
                    if not valid:break
                    queue.append(neighbor)
                if not valid:break
            if valid and len(mapping)==n and len(set(mapping.values()))==n:
                perms.append(Permutation([mapping[i] for i in range(n)]))
    generators=[];covered={Permutation(range(n))}
    for p in perms:
        if p not in covered:
            generators.append(p);covered=set(PermutationGroup(generators)._enumerate())
    return generators


# Generators of the reference pgroups, in the reference order, so that
# rotate(i) uses the same rotation as SymPy.
_REFERENCE_GENERATORS={'tetrahedron':[[0,2,3,1],[1,2,0,3],[3,1,0,2],[3,0,2,1],[1,0,3,2],[2,3,0,1],[3,2,1,0]],'cube':[[1,2,3,0,5,6,7,4],[4,0,3,7,5,1,2,6],[4,5,1,0,7,6,2,3],[1,0,4,5,2,3,7,6],[6,2,1,5,7,3,0,4],[6,7,3,2,5,4,0,1],[3,7,4,0,2,6,5,1],[4,7,6,5,0,3,2,1],[6,5,4,7,2,1,0,3],[0,3,7,4,1,2,6,5],[5,1,0,4,6,2,3,7],[5,6,2,1,4,7,3,0],[7,4,0,3,6,5,1,2]],'octahedron':[[0,2,3,4,1,5],[4,1,0,3,5,2],[1,5,2,0,4,3],[1,0,4,5,2,3],[2,3,0,1,5,4],[3,5,4,0,2,1],[4,3,5,1,0,2],[5,4,3,2,1,0],[5,2,1,4,3,0],[4,0,3,5,1,2],[1,2,0,4,5,3],[2,5,3,0,1,4],[4,5,1,0,3,2]],'dodecahedron':[[1,2,3,4,0,6,7,8,9,5,11,12,13,14,10,16,17,18,19,15],[5,0,4,9,14,10,1,3,13,15,6,2,8,19,16,17,11,7,12,18],[10,6,1,0,5,16,11,2,4,14,17,7,3,9,15,19,18,12,8,13],[6,11,7,2,1,10,17,12,3,0,16,18,8,4,5,14,15,19,13,9],[2,7,12,8,3,1,11,18,13,4,6,17,19,9,0,5,10,16,15,14],[4,3,8,13,9,0,2,12,19,14,1,7,18,15,5,10,6,11,17,16],[1,0,5,10,6,2,4,14,16,11,3,9,15,17,7,12,8,13,19,18],[7,2,1,6,11,12,3,0,10,17,8,4,5,16,18,19,13,9,14,15],[12,8,3,2,7,18,13,4,1,11,19,9,0,6,17,16,15,14,5,10],[8,13,9,4,3,12,19,14,0,2,18,15,5,1,7,11,17,16,10,6],[4,9,14,5,0,3,13,15,10,1,8,19,16,6,2,7,12,18,17,11],[5,14,15,16,10,0,9,19,17,6,4,13,18,11,1,2,3,8,12,7],[11,6,10,16,17,7,1,5,15,18,2,0,14,19,12,8,3,4,9,13],[18,12,7,11,17,19,8,2,6,16,13,3,1,10,15,14,9,4,0,5],[18,19,13,8,12,17,15,9,3,7,16,14,4,2,11,6,10,5,0,1],[13,19,15,14,9,8,18,16,5,4,12,17,10,0,3,2,7,11,6,1],[16,15,19,18,17,10,14,13,12,11,5,9,8,7,6,1,0,4,3,2],[18,17,16,15,19,12,11,10,14,13,7,6,5,9,8,3,2,1,0,4],[15,19,18,17,16,14,13,12,11,10,9,8,7,6,5,0,4,3,2,1],[17,16,15,19,18,11,10,14,13,12,6,5,9,8,7,2,1,0,4,3],[19,18,17,16,15,13,12,11,10,14,8,7,6,5,9,4,3,2,1,0],[0,4,9,14,5,1,3,13,15,10,2,8,19,16,6,11,7,12,18,17],[6,1,0,5,10,11,2,4,14,16,7,3,9,15,17,18,12,8,13,19],[11,7,2,1,6,17,12,3,0,10,18,8,4,5,16,15,19,13,9,14],[7,12,8,3,2,11,18,13,4,1,17,19,9,0,6,10,16,15,14,5],[3,8,13,9,4,2,12,19,14,0,7,18,15,5,1,6,11,17,16,10],[14,9,13,19,15,5,4,8,18,16,0,3,12,17,10,6,1,2,7,11],[16,10,5,14,15,17,6,0,9,19,11,1,4,13,18,12,7,2,3,8],[16,17,11,6,10,15,18,7,1,5,19,12,2,0,14,9,13,8,3,4],[11,17,18,12,7,6,16,19,8,2,10,15,13,3,1,0,5,14,9,4],[8,12,18,19,13,3,7,17,15,9,2,11,16,14,4,0,1,6,10,5]],'icosahedron':[[0,2,3,4,5,1,7,8,9,10,6,11],[5,1,0,4,10,6,7,2,3,9,11,8],[1,7,2,0,5,6,11,8,3,4,10,9],[2,7,8,3,0,1,6,11,9,4,5,10],[3,2,8,9,4,0,1,7,11,10,5,6],[4,0,3,9,10,5,1,2,8,11,6,7],[1,0,5,6,7,2,3,4,10,11,8,9],[2,3,0,1,7,8,9,4,5,6,11,10],[3,9,4,0,2,8,11,10,5,1,7,6],[4,9,10,5,0,3,8,11,6,1,2,7],[5,4,10,6,1,0,3,9,11,7,2,8],[6,5,10,11,7,1,0,4,9,8,2,3],[7,2,1,6,11,8,3,0,5,10,9,4],[8,9,3,2,7,11,10,4,0,1,6,5],[9,11,10,4,3,8,7,6,5,0,2,1],[10,9,11,6,5,4,3,8,7,1,0,2],[11,6,10,9,8,7,1,5,4,3,2,0],[11,8,7,6,10,9,3,2,1,5,4,0],[11,10,9,8,7,6,5,4,3,2,1,0],[11,7,6,10,9,8,2,1,5,4,3,0],[11,9,8,7,6,10,4,3,2,1,5,0],[5,0,4,10,6,1,2,3,9,11,7,8],[1,2,0,5,6,7,8,3,4,10,11,9],[2,8,3,0,1,7,11,9,4,5,6,10],[3,8,9,4,0,2,7,11,10,5,1,6],[4,3,9,10,5,0,2,8,11,6,1,7],[10,5,4,9,11,6,1,0,3,8,7,2],[6,7,1,5,10,11,8,2,0,4,9,3],[7,11,8,2,1,6,10,9,3,0,5,4],[8,7,11,9,3,2,1,6,10,4,0,5],[9,3,8,11,10,4,0,2,7,6,5,1]]}


def _regular(n,faces,name=None):
    generators=[Permutation(a) for a in _REFERENCE_GENERATORS[name]] if name else _rotations(n,faces)
    return Polyhedron(range(n),faces,generators)


tetrahedron=_regular(4,[(0,1,2),(0,1,3),(0,2,3),(1,2,3)],'tetrahedron')
cube=_regular(8,[(0,1,2,3),(0,1,5,4),(0,3,7,4),(1,2,6,5),(2,3,7,6),(4,5,6,7)],'cube')
octahedron=_regular(6,[(0,1,2),(0,1,4),(0,2,3),(0,3,4),(1,2,5),(1,4,5),(2,3,5),(3,4,5)],'octahedron')
dodecahedron=_regular(20,[(0,1,2,3,4),(0,1,6,10,5),(0,4,9,14,5),(1,2,7,11,6),(2,3,8,12,7),(3,4,9,13,8),(5,10,16,15,14),(6,10,16,17,11),(7,11,17,18,12),(8,12,18,19,13),(9,13,19,15,14),(15,16,17,18,19)],'dodecahedron')
icosahedron=_regular(12,[(0,1,2),(0,1,5),(0,2,3),(0,3,4),(0,4,5),(1,2,7),(1,5,6),(1,6,7),(2,3,8),(2,7,8),(3,4,9),(3,8,9),(4,5,10),(4,9,10),(5,6,10),(6,7,11),(6,10,11),(7,8,11),(8,9,11),(9,10,11)],'icosahedron')
