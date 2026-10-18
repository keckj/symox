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

from ._group_data import GENERATORS as _GROUP_DATA


def _group_constructor(name):
    def construct():return PermutationGroup([Permutation(array) for array in _GROUP_DATA[name]])
    construct.__name__=name;construct.__doc__='Construct the standard permutation representation of '+name+'.'
    return construct
for _name in ['four_group','M20','S3_in_S6','A4_in_S6','S4m','S4p','A4xC2','S4xC2','G18','G36m','G36p','G72','PSL2F5','PGL2F5']:
    globals()[_name]=_group_constructor(_name)

class S5TransitiveSubgroups(_TransitiveSubgroups):
    C5='C5';D5='D5';M20='M20';A5='A5';S5='S5'
class S6TransitiveSubgroups(_TransitiveSubgroups):
    C6='C6';S3='S3';D6='D6';A4='A4';G18='G18';A4xC2='A4 x C2';S4m='S4-';S4p='S4+';G36m='G36-';G36p='G36+';S4xC2='S4 x C2';PSL2F5='PSL2(F5)';G72='G72';PGL2F5='PGL2(F5)';A6='A6';S6='S6'
_GENERATORS.update({key:value for key,value in _GROUP_DATA.items() if isinstance(key,tuple)})


def find_transitive_subgroups_of_S6(*targets,print_report=False):
    """Find deterministic catalog representatives and verify transitivity."""
    found={}
    for target in targets:
        if not isinstance(target,S6TransitiveSubgroups):raise ValueError('targets must be names of transitive S6 subgroups')
        group=target.get_perm_group()
        if not group.is_transitive():raise ArithmeticError('invalid transitive subgroup representation')
        found[target]=group
        if print_report:print(str(target)+': '+str(group.generators))
    return found
