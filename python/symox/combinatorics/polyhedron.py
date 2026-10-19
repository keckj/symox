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
        self._original=tuple(corners);self._corners=self._original
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
    def copy(self):return copy.copy(self)


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


def _regular(n,faces):return Polyhedron(range(n),faces,_rotations(n,faces))


tetrahedron=_regular(4,[(0,1,2),(0,1,3),(0,2,3),(1,2,3)])
cube=_regular(8,[(0,1,2,3),(0,1,5,4),(0,3,7,4),(1,2,6,5),(2,3,7,6),(4,5,6,7)])
octahedron=_regular(6,[(0,1,2),(0,1,4),(0,2,3),(0,3,4),(1,2,5),(1,4,5),(2,3,5),(3,4,5)])
dodecahedron=_regular(20,[(0,1,2,3,4),(0,1,6,10,5),(0,4,9,14,5),(1,2,7,11,6),(2,3,8,12,7),(3,4,9,13,8),(5,10,16,15,14),(6,10,16,17,11),(7,11,17,18,12),(8,12,18,19,13),(9,13,19,15,14),(15,16,17,18,19)])
icosahedron=_regular(12,[(0,1,2),(0,1,5),(0,2,3),(0,3,4),(0,4,5),(1,2,7),(1,5,6),(1,6,7),(2,3,8),(2,7,8),(3,4,9),(3,8,9),(4,5,10),(4,9,10),(5,6,10),(6,7,11),(6,10,11),(7,8,11),(8,9,11),(9,10,11)])
