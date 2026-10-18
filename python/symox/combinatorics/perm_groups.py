"""Small finite permutation groups for Galois results.

Group enumeration is exact and bounded at 100,000 elements. This is not an
implementation of the full Schreier-Sims machinery for large groups.
"""
from .permutations import Permutation


class PermutationGroup:
    def __init__(self, *generators):
        if len(generators) == 1 and isinstance(generators[0], (list, tuple)):
            generators = tuple(generators[0])
        if not generators or any(not isinstance(g, Permutation) for g in generators):
            raise ValueError('Expected one or more permutation generators')
        self.degree = max(g.size for g in generators)
        self._generators = tuple(Permutation(g, size=self.degree) for g in generators)
        self._elements = None

    @property
    def generators(self):
        return list(self._generators)

    @property
    def identity(self):
        return Permutation(range(self.degree))

    def _enumerate(self):
        if self._elements is None:
            elements = {self.identity}
            queue = [self.identity]
            for value in queue:
                for generator in self._generators:
                    result = value*generator
                    if result not in elements:
                        if len(elements) >= 100_000:
                            raise NotImplementedError('Permutation group enumeration limit exceeded')
                        elements.add(result)
                        queue.append(result)
            self._elements = tuple(queue)
        return self._elements

    @property
    def elements(self):
        return list(self._enumerate())

    def order(self):
        return len(self._enumerate())

    def contains(self, permutation, strict=True):
        if not isinstance(permutation, Permutation):
            return False
        if permutation.size != self.degree:
            if strict or permutation.size > self.degree:
                return False
            permutation = Permutation(permutation, size=self.degree)
        return permutation in set(self._enumerate())

    __contains__ = contains

    def is_transitive(self, strict=True):
        orbit = {g(0) for g in self._enumerate()} if self.degree else set()
        return len(orbit) == self.degree

    @property
    def is_abelian(self):
        return all(a*b == b*a for a in self._generators for b in self._generators)

    @property
    def is_solvable(self):
        group = self
        while group.order() != 1:
            if group.order() > 256:
                raise NotImplementedError('Derived subgroup computation limit exceeded')
            elements = group._enumerate()
            # All commutators generate the derived subgroup.
            derived = PermutationGroup([~a*~b*a*b for a in elements for b in elements])
            if derived.order() == group.order():
                return False
            group = derived
        return True

    def __eq__(self, other):
        if not isinstance(other, PermutationGroup):
            return NotImplemented
        return self._generators == other._generators

    def __repr__(self):
        return 'PermutationGroup(%r)' % self.generators
