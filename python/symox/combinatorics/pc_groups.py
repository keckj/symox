"""Collection in certified finite polycyclic presentations.

Normal forms are enumerated up to 100000 elements. Infinite polycyclic
presentations and large-group collection require a different backend.
"""
from itertools import product
from math import prod
import operator
from .free_groups import free_group, FreeGroupElement
from .perm_groups import PermutationGroup
from .permutations import Permutation


class Collector:
    def __init__(self,pcgs,pc_series,relative_order,free_group_=None,pc_presentation=None):
        self.pcgs=list(pcgs);self.pc_series=list(pc_series)
        self.relative_order=[operator.index(v) for v in relative_order]
        if len(self.pcgs)!=len(self.relative_order) or any(v<1 for v in self.relative_order):
            raise ValueError('finite relative orders must be positive and match the generators')
        if prod(self.relative_order)>100000:raise NotImplementedError('polycyclic normal form enumeration exceeds 100000 elements')
        if any(not isinstance(g,Permutation) for g in self.pcgs):raise TypeError('finite permutation generators required')
        self.group=PermutationGroup(self.pcgs or [self.pc_series[0].identity if self.pc_series else Permutation([])])
        self.free_group=free_group('x:%d'%len(self.pcgs))[0] if free_group_ is None else free_group_
        if self.free_group.rank!=len(self.pcgs):raise ValueError('free group has the wrong rank')
        self.index={s:i for i,s in enumerate(self.free_group.symbols)}
        self._vectors={}
        for vector in product(*(range(n) for n in self.relative_order)):
            value=self.group.identity
            for g,e in zip(self.pcgs,vector):value=value*g**e
            if value in self._vectors:raise ValueError('polycyclic normal forms are not unique')
            self._vectors[value]=vector
        if len(self._vectors)!=self.group.order():raise ValueError('normal forms do not cover the group')
        # Certify cyclic quotient factors and normality of each tail subgroup.
        tails=[PermutationGroup(self.pcgs[i:] or [self.group.identity]) for i in range(len(self.pcgs)+1)]
        for i,g in enumerate(self.pcgs):
            if tails[i].order()!=self.relative_order[i]*tails[i+1].order() or any(not tails[i+1].contains(~g*h*g,strict=False) for h in tails[i+1].generators):
                raise ValueError('the generators are not a polycyclic sequence')
        if self.pc_series:
            if len(self.pc_series)!=len(tails) or any(set(a._enumerate())!=set(b._enumerate()) for a,b in zip(self.pc_series,tails)):
                raise ValueError('series does not match the polycyclic sequence')
        self.pc_series=tails
        self.pc_presentation=self.pc_relators()
        if pc_presentation is not None and pc_presentation!=self.pc_presentation:
            raise ValueError('supplied relators disagree with the permutation presentation')

    def exponent_vector(self,element):
        if not isinstance(element,Permutation):raise TypeError('expected a permutation')
        element=Permutation(element,size=self.group.degree)
        try:return list(self._vectors[element])
        except KeyError:raise ValueError('element is outside the polycyclic group') from None

    def _word(self,element):
        word=self.free_group.identity
        for g,e in zip(self.free_group.generators,self.exponent_vector(element)):word*=g**e
        return word

    def _permutation(self,word):
        if isinstance(word,tuple) and not word:return self.group.identity
        if not isinstance(word,FreeGroupElement) or word.group!=self.free_group:raise ValueError('word belongs to another free group')
        value=self.group.identity
        for s,e in word.array_form:value*=self.pcgs[self.index[s]]**e
        return value

    def collected_word(self,word):return self._word(self._permutation(word))

    def pc_relators(self):
        # Reference layout: generators are processed from the last one,
        # identity values are stored as the empty tuple ().
        result={}
        gens=self.free_group.generators
        def value(element):
            word=self._word(element)
            return word if word else ()
        for i in reversed(range(len(self.pcgs))):
            g,n=self.pcgs[i],self.relative_order[i]
            result[gens[i]**n]=value(g**n)
            for j in reversed(range(i+1,len(gens))):
                result[gens[i]**-1*gens[j]*gens[i]]=value(~g*self.pcgs[j]*g)
        return result

    def relations(self):
        return ({k:v for k,v in self.pc_presentation.items() if len(k.array_form)==1},
                {k:v for k,v in self.pc_presentation.items() if len(k.array_form)!=1})

    def minimal_uncollected_subword(self,word):
        self._permutation(word)
        for s,e in word.array_form:
            if not 0<e<self.relative_order[self.index[s]]:return ((s,e),)
        for (s,e),(t,f) in zip(word.array_form,word.array_form[1:]):
            if self.index[s]>self.index[t]:return ((s,e),(t,1 if f>0 else -1))
        return None

    def subword_index(self,word,w):
        # Reference: letter positions (low, high) of the first occurrence of
        # ``w`` in ``word``, or (-1, -1).
        low=-1
        high=-1
        for i in range(len(word)-len(w)+1):
            if word.subword(i,i+len(w))==w:
                low=i
                high=i+len(w)
                break
        return low,high

    def map_relation(self,w):
        # Reference: for w = x_i*x_j return the conjugate relator value of
        # x_j**-1*x_i*x_j from the polycyclic presentation.
        array=w.array_form
        s1=array[0][0]
        s2=array[1][0]
        key=self.free_group.dtype(((s2,-1),(s1,1),(s2,1)))
        return self.pc_presentation[key]

    def depth(self,element):
        return next((i+1 for i,v in enumerate(self.exponent_vector(element)) if v),len(self.pcgs)+1)

    def leading_exponent(self,element):
        return next((v for v in self.exponent_vector(element) if v),None)

    def induced_pcgs(self,gens):
        subgroup=PermutationGroup(gens)
        chosen=[]
        for d in range(1,len(self.pcgs)+1):
            candidates=[g for g in subgroup._enumerate() if self.depth(g)==d]
            if candidates:
                chosen.append(min(candidates,key=lambda g:self.leading_exponent(g)))
        return chosen

    def constructive_membership_test(self,ipcgs,g):
        group=PermutationGroup(ipcgs or [self.group.identity])
        return group.contains(g,strict=False)


class PolycyclicGroup:
    is_group=is_solvable=True
    def __init__(self,pc_sequence,pc_series,relative_order,collector=None):
        self.collector=Collector(pc_sequence,pc_series,relative_order) if collector is None else collector
        self.pcgs=self.collector.pcgs;self.pc_series=self.collector.pc_series;self.relative_order=self.collector.relative_order
    def length(self):return len(self.pcgs)
    def is_prime_order(self):
        from .. import isprime
        return all(isprime(n) for n in self.relative_order)
