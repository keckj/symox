"""Native finite combinatorics and exact word arithmetic."""
from .permutations import Permutation, Cycle
from .perm_groups import PermutationGroup
from .generators import cyclic, alternating, symmetric, dihedral
from .named_groups import (SymmetricGroup, AlternatingGroup, CyclicGroup,
    DihedralGroup, AbelianGroup, DirectProduct, RubikGroup, SymmetricPermutationGroup)
from .graycode import GrayCode
from .subsets import Subset
from .partitions import Partition, IntegerPartition, RGS_enum, RGS_rank, RGS_unrank
from .prufer import Prufer
from .coset import Coset
from .polyhedron import Polyhedron, tetrahedron, cube, octahedron, dodecahedron, icosahedron

__all__ = ['Permutation','PermutationGroup','Cycle','cyclic','alternating','symmetric',
    'dihedral','SymmetricGroup','AlternatingGroup','CyclicGroup','DihedralGroup',
    'AbelianGroup','DirectProduct','RubikGroup','SymmetricPermutationGroup',
    'GrayCode','Subset','Partition','IntegerPartition','RGS_enum','RGS_rank',
    'RGS_unrank','Prufer','free_group','Coset','PolycyclicGroup','Collector',
    'Polyhedron','tetrahedron','cube','octahedron','dodecahedron','icosahedron']

# Free, finitely presented and polycyclic groups import the core, which is
# still initializing when this package loads; resolve them on first access.
from symox._api_exports import export as _export
_export(globals(), 'symox.combinatorics.free_groups', 'free_group')
_export(globals(), 'symox.combinatorics.pc_groups', 'PolycyclicGroup Collector')
