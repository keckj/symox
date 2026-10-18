"""Lazy enumeration of the elementary finite permutation groups."""
import operator
from itertools import permutations
from .permutations import Permutation


def _degree(n,minimum=0):
    n=operator.index(n)
    if n<minimum:
        raise ValueError('invalid group degree')
    return n


def symmetric(n):
    for images in permutations(range(_degree(n))):
        yield Permutation(images)


def alternating(n):
    yield from (g for g in symmetric(n) if g.is_even)


def cyclic(n):
    n=_degree(n,1)
    for k in range(n):
        yield Permutation([(i+k)%n for i in range(n)])


def dihedral(n):
    n=_degree(n,1)
    if n==1:
        yield from symmetric(2)
    elif n==2:
        for a in ((0,1,2,3),(1,0,3,2),(2,3,0,1),(3,2,1,0)):
            yield Permutation(a)
    else:
        for p in cyclic(n):
            yield p
            yield Permutation(p.array_form[::-1])
