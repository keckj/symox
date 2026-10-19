"""Exact finite permutation groups with incremental stabilizer chains."""
from collections import deque
import math
import operator
import random

from .permutations import Permutation,Cycle
from ._schreier import StabilizerChain,multiply,inverse


class PermutationGroup:
    def __init__(self, *generators,dups=True,**kwargs):
        if generators and isinstance(generators[0], (list, tuple)):
            generators = tuple(generators[0])
        if not generators:generators=(Permutation([]),)
        if any(isinstance(g,Cycle) for g in generators):generators=tuple(Permutation(g) for g in generators)
        if any(not isinstance(g, Permutation) for g in generators):
            raise ValueError('Expected one or more permutation generators')
        self.degree = max(g.size for g in generators)
        if len({g.size for g in generators})>1:
            degree=kwargs.pop('degree',None)
            if degree is not None:self.degree=operator.index(degree)
        if kwargs:raise TypeError('unexpected group constructor keyword: '+next(iter(kwargs)))
        generators=[Permutation(g,size=self.degree) for g in generators]
        if dups:generators=list(dict.fromkeys(generators))
        if len(generators)>1:generators=[g for g in generators if not g.is_Identity]
        self._generators=tuple(generators or [Permutation(range(self.degree))])
        self._elements = None
        self._chain_cache=None
        self._order=None

    @property
    def generators(self):
        return list(self._generators)

    @property
    def identity(self):
        return Permutation(range(self.degree))

    def _enumerate(self):
        if self._elements is None:
            self._elements=tuple(self.generate_schreier_sims())
        return self._elements

    def _chain(self):
        if self._chain_cache is None:
            self._chain_cache=StabilizerChain([g._images for g in self._generators],self.degree)
        return self._chain_cache

    def schreier_sims(self):self._chain()

    @property
    def base(self):return list(self._chain().base)
    @property
    def strong_gens(self):return [Permutation(g) for g in self._chain().strong] or [self.identity]
    @property
    def basic_orbits(self):return [sorted(level) for level in self._chain().transversals]
    @property
    def basic_transversals(self):return [{point:Permutation(g) for point,g in level.items()}
        for level in self._chain().transversals]
    @property
    def basic_stabilizers(self):return [PermutationGroup([Permutation(g) for g in generators]
        or [self.identity]) for generators in self._chain().level_generators]

    def schreier_sims_incremental(self,base=None,gens=None,slp_dict=False):
        generators=self.generators if gens is None else list(gens)
        generators=[Permutation(g,size=self.degree) for g in generators]
        if len(generators)==1 and generators[0].is_Identity:
            result=([] if base is None else list(base),generators)
            return result+({generators[0]:[generators[0]]},) if slp_dict else result
        originals=[g for g in generators if not g.is_Identity]
        if base is None and gens is None and not slp_dict:chain=self._chain()
        else:
            chain=StabilizerChain([g._images for g in originals],self.degree,
                () if base is None else base,track_words=bool(slp_dict))
        initial=set(g._images for g in originals)
        strong=[Permutation(g) for g in chain.strong]
        result=(list(chain.base),originals+[g for g in strong if g._images not in initial])
        if not slp_dict:return result
        # Public SLPs use rmul order (opposite to our right-action words).
        inverses={}
        def factor(term):
            if term>0:return strong[term-1]
            if term not in inverses:inverses[term]=~strong[-term-1]
            return inverses[term]
        witnesses={g:[factor(term) for term in reversed(word)]
            for g,word in zip(strong,chain._words.words)}
        return result+(witnesses,)

    @property
    def elements(self):
        return list(self._enumerate())

    def order(self):
        if self._order is None:
            self._order=(self._generators[0].order() if len(self._generators)==1 else self._chain().order)
        return self._order

    def contains(self, g, strict=True):
        permutation=g
        if not isinstance(permutation, Permutation):
            return False
        if permutation.size != self.degree:
            if strict:
                return False
            permutation = Permutation(permutation, size=self.degree)
        if permutation.is_Identity or permutation in self._generators:return True
        if getattr(self,'_is_symmetric',False):return True
        if getattr(self,'_is_alternating',False):return permutation.is_even
        if len(self._generators)==1:
            # Powers of one generator act by a uniform rotation on each
            # disjoint cycle. Generalized CRT checks that the shifts describe
            # one common exponent, including noncoprime cycle lengths.
            residue,modulus=0,1
            for cycle in self._generators[0].full_cyclic_form:
                positions={point:index for index,point in enumerate(cycle)}
                shift=positions.get(permutation(cycle[0]))
                if shift is None:return False
                size=len(cycle)
                if any(permutation(point)!=cycle[(index+shift)%size] for index,point in enumerate(cycle)):
                    return False
                divisor=math.gcd(modulus,size)
                if (shift-residue)%divisor:return False
                reduced=size//divisor
                multiplier=((shift-residue)//divisor*pow(modulus//divisor,-1,reduced))%reduced
                residue+=modulus*multiplier;modulus*=reduced
            return True
        return self._chain().contains(permutation._images)

    __contains__ = contains

    def is_transitive(self, strict=True):
        if strict and hasattr(self,'_is_transitive'):return self._is_transitive
        if not strict:return sum(len(orbit)>1 for orbit in self.orbits())==1
        self._is_transitive=not self.degree or len(self.orbit(0))==self.degree
        return self._is_transitive

    @property
    def is_abelian(self):
        if hasattr(self,'_is_abelian'):return self._is_abelian
        if len(self._generators)==1:return True
        return all(a*b == b*a for a in self._generators for b in self._generators)

    @property
    def is_solvable(self):
        if hasattr(self,'_is_solvable'):return self._is_solvable
        if len(self._generators)==1:return True
        self._is_solvable=self.derived_series()[-1].order()==1
        return self._is_solvable

    def _normal_closure(self,generators):
        generators=list(dict.fromkeys(g for g in generators if not g.is_Identity))
        subgroup=PermutationGroup(generators or [self.identity]);pending=deque(generators)
        conjugators=self.generators+[~g for g in self._generators]
        while pending:
            generator=pending.popleft()
            for conjugator in conjugators:
                candidate=~conjugator*generator*conjugator
                if not subgroup.contains(candidate):
                    generators.append(candidate);pending.append(candidate)
                    subgroup=PermutationGroup(generators)
        return subgroup

    def derived_subgroup(self):
        return self._normal_closure([~a*~b*a*b for a in self._generators for b in self._generators])

    def derived_series(self):
        series=[self];group=self
        while group.order()!=1:
            derived=group.derived_subgroup()
            if derived.order()==group.order():break
            series.append(derived);group=derived
        return series

    def lower_central_series(self):
        series=[self];group=self
        while group.order()!=1:
            commutators=[~a*~b*a*b for a in self._generators for b in group._generators]
            following=self._normal_closure(commutators)
            if following.order()==group.order():break
            series.append(following);group=following
        return series

    @property
    def is_nilpotent(self):
        if hasattr(self,'_is_nilpotent'):return self._is_nilpotent
        if len(self._generators)==1:return True
        self._is_nilpotent=self.lower_central_series()[-1].order()==1
        return self._is_nilpotent

    def __eq__(self, other):
        if not isinstance(other, PermutationGroup):
            return NotImplemented
        return self._generators == other._generators

    def __repr__(self):
        return 'PermutationGroup(%r)' % self.generators

    is_group = True

    def __getitem__(self,index):
        return self.generators[index]

    def __len__(self):return len(self._generators)

    def generate_schreier_sims(self,af=False):
        for element in self._chain().generate():
            yield list(element) if af else Permutation(element)

    def generate(self, method='coset', af=False):
        """Iterate group elements by stabilizer cosets or Dimino's algorithm."""
        if method == 'coset':
            return self.generate_schreier_sims(af)
        if method == 'dimino':
            return self.generate_dimino(af)
        raise NotImplementedError('No generation defined for %s' % method)

    def generate_dimino(self, af=False):
        """Enumerate incremental subgroup cosets with Dimino's algorithm."""
        identity = tuple(range(self.degree))
        elements, seen = [identity], {identity}
        yield list(identity) if af else Permutation(identity)
        generators = [g._images for g in self._generators]
        for index in range(len(generators)):
            subgroup = elements.copy()
            pending = [identity]
            while pending:
                frontier, pending = pending, []
                for representative in frontier:
                    for generator in generators[:index+1]:
                        coset = multiply(generator, representative)
                        if coset in seen:
                            continue
                        for element in subgroup:
                            image = multiply(coset, element)
                            seen.add(image)
                            elements.append(image)
                            pending.append(image)
                            yield list(image) if af else Permutation(image)
        self._order = len(elements)

    def random(self, af=False):
        """Sample uniformly from the unique stabilizer-transversal products."""
        element = tuple(range(self.degree))
        for transversal in self._chain().transversals:
            element = multiply(random.choice(tuple(transversal.values())), element)
        return list(element) if af else Permutation(element)

    def conjugacy_class(self, x):
        """Find the conjugation orbit without enumerating the entire group."""
        result, pending = {x}, deque([x])
        conjugators = [(g, ~g) for g in self._generators]
        while pending:
            element = pending.popleft()
            for generator, inverse_generator in conjugators:
                image = generator*element*inverse_generator
                if image not in result:
                    result.add(image)
                    pending.append(image)
        return result

    def conjugacy_classes(self):
        """Partition the group into conjugation orbits, identity first."""
        seen = {self.identity}
        classes = [seen.copy()]
        for element in self.generate():
            if element not in seen:
                orbit = self.conjugacy_class(element)
                classes.append(orbit)
                seen.update(orbit)
        return classes

    def _orbit_transversal(self,alpha):
        alpha=operator.index(alpha)
        if not 0<=alpha<self.degree:raise IndexError('point outside the permutation domain')
        representatives={alpha:tuple(range(self.degree))};queue=deque([alpha])
        for_point_generators=[g._images for g in self._generators]
        while queue:
            point=queue.popleft()
            for g in for_point_generators:
                image=g[point]
                if image not in representatives:
                    representatives[image]=multiply(representatives[point],g);queue.append(image)
        return representatives

    def orbit_transversal(self,alpha,pairs=False):
        representatives=self._orbit_transversal(alpha)
        if pairs:return [(point,Permutation(g)) for point,g in representatives.items()]
        return [Permutation(g) for g in representatives.values()]

    def _point_orbit(self,points):
        orbit={operator.index(point) for point in points}
        if any(not 0<=point<self.degree for point in orbit):raise IndexError('point outside the permutation domain')
        queue=deque(orbit)
        while queue:
            value=queue.popleft()
            for g in self._generators:
                image=g(value)
                if image not in orbit:orbit.add(image);queue.append(image)
        return orbit

    def orbit(self,alpha,action='tuples'):
        try:point=operator.index(alpha)
        except TypeError:pass
        else:return self._point_orbit((point,))
        if action=='union':return self._point_orbit(alpha)
        if action not in ('tuples','sets'):raise ValueError('unknown orbit action')
        normalize=frozenset if action=='sets' else tuple
        state=normalize(alpha);orbit={state};queue=deque([state])
        while queue:
            state=queue.popleft()
            for g in self._generators:
                image=normalize(g(point) for point in state)
                if image not in orbit:orbit.add(image);queue.append(image)
        return {tuple(state) for state in orbit} if action=='sets' else orbit

    def orbits(self,rep=False):
        unseen=set(range(self.degree));result=[]
        while unseen:
            orbit=self.orbit(min(unseen));unseen.difference_update(orbit)
            result.append(min(orbit) if rep else orbit)
        return result

    def stabilizer(self,alpha):
        representatives=self._orbit_transversal(alpha)
        inverses={point:inverse(g) for point,g in representatives.items()}
        identity=tuple(range(self.degree));generators=[]
        for point,representative in representatives.items():
            for g in self._generators:
                schreier=multiply(multiply(representative,g._images),inverses[g(point)])
                if schreier!=identity:generators.append(Permutation(schreier))
        return PermutationGroup(generators or [self.identity])

    def make_perm(self,n,seed=None):
        if isinstance(n,(tuple,list,range)):
            if seed is not None:raise ValueError('a sequence of generator indices requires seed=None')
            n,seed=len(n),n
        else:
            try:n=int(n)
            except TypeError:raise ValueError('n must be an integer or sequence') from None
        if isinstance(seed,(list,tuple)):
            stream=iter(seed)
            def choose():
                for index in stream:
                    if 0<=index<len(self):return index
                raise ValueError('generator-index sequence is too short')
        elif seed is None:
            from ..core.random import randrange
            choose=lambda:randrange(len(self))
        else:
            rng=random.Random(seed);choose=lambda:rng.randrange(len(self))
        result=self.identity
        for _ in range(n):result=self._generators[choose()]*result
        return result

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

PermGroup=PermutationGroup
