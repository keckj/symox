"""Independent finite graph and tensor checks for combinatorics clients."""
from itertools import permutations
from math import factorial


def canonicalize_naive(g,dummies,sym,*v):
    # This public utility shares the explicit orbit algorithm; both paths are
    # bounded exact enumeration, with no claim of a second independent solver.
    from .tensor_can import canonicalize
    return canonicalize(g,dummies,sym,*v)


def graph_certificate(gr):
    vertices=sorted(gr);n=len(vertices)
    if factorial(n)>1000000:raise ValueError('graph labeling exceeds one million permutations')
    edges={vertex:set(gr[vertex]) for vertex in vertices}
    if any(neighbor not in edges for values in edges.values() for neighbor in values):raise ValueError('all adjacency vertices must be declared')
    return min(tuple(int(b in edges[a]) for a in order for b in order) for order in permutations(vertices))

__all__=['canonicalize_naive','graph_certificate']
