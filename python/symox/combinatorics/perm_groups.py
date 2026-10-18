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

    is_group = True

    def __getitem__(self,index):
        return self.generators[index]

    def generate_schreier_sims(self,af=False):
        for element in self._enumerate():
            yield element.array_form if af else element

    def orbit(self,alpha):
        return {g(alpha) for g in self._enumerate()}

    def subgroup(self,gens):
        group=PermutationGroup(gens or [self.identity])
        if any(not self.contains(g,strict=False) for g in group.generators):raise ValueError('generators are outside the group')
        return group

    def polycyclic_group(self):
        from .pc_groups import PolycyclicGroup
        if self.order()>256:raise NotImplementedError('polycyclic series construction is limited to 256 elements')
        if not self.is_solvable:raise ValueError('group is not solvable')
        # Refine a derived series through maximal normal subgroups with
        # cyclic quotient, checking each finite subgroup exactly.
        sequence=[];series=[self];orders=[];group=self
        while group.order()>1:
            elems=group._enumerate();derived=PermutationGroup([~a*~b*a*b for a in elems for b in elems])
            tail=derived
            while True:
                candidate=None
                for g in elems:
                    if tail.contains(g):continue
                    enlarged=PermutationGroup(tail.generators+[g])
                    if enlarged.order()<group.order():candidate=enlarged;break
                if candidate is None:break
                tail=candidate
            generator=next(g for g in elems if not tail.contains(g))
            sequence.append(generator);orders.append(group.order()//tail.order());series.append(tail);group=tail
        return PolycyclicGroup(sequence,series,orders)

# Compatible imports reuse the existing native objects.
from symox._api_exports import export as _export
_export(globals(), 'symox.combinatorics.coset', 'Coset')
_export(globals(), 'symox.combinatorics.named_groups', 'SymmetricPermutationGroup')
