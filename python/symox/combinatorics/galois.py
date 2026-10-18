"""Names and permutation representatives for transitive groups of degree <= 4."""
from enum import Enum

from .permutations import Permutation
from .perm_groups import PermutationGroup


class _TransitiveSubgroups(Enum):
    def get_perm_group(self):
        arrays = _GENERATORS[type(self).__name__, self.name]
        return PermutationGroup([Permutation(array) for array in arrays])


class S1TransitiveSubgroups(_TransitiveSubgroups):
    S1 = 'S1'


class S2TransitiveSubgroups(_TransitiveSubgroups):
    S2 = 'S2'


class S3TransitiveSubgroups(_TransitiveSubgroups):
    A3 = 'A3'
    S3 = 'S3'


class S4TransitiveSubgroups(_TransitiveSubgroups):
    C4 = 'C4'
    V = 'V'
    D4 = 'D4'
    A4 = 'A4'
    S4 = 'S4'


_GENERATORS = {
    ('S1TransitiveSubgroups', 'S1'): [[0]],
    ('S2TransitiveSubgroups', 'S2'): [[1, 0]],
    ('S3TransitiveSubgroups', 'A3'): [[1, 2, 0]],
    ('S3TransitiveSubgroups', 'S3'): [[1, 2, 0], [1, 0, 2]],
    ('S4TransitiveSubgroups', 'C4'): [[1, 2, 3, 0]],
    ('S4TransitiveSubgroups', 'V'): [[1, 0, 3, 2], [2, 3, 0, 1]],
    ('S4TransitiveSubgroups', 'D4'): [[1, 2, 3, 0], [3, 2, 1, 0]],
    ('S4TransitiveSubgroups', 'A4'): [[1, 2, 0, 3], [0, 2, 3, 1]],
    ('S4TransitiveSubgroups', 'S4'): [[1, 2, 3, 0], [1, 0, 2, 3]],
}
