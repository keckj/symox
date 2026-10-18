"""Finite subset ranking in binary, lexicographic and Gray orders."""
import operator
from .graycode import _rank


class Subset:
    def __init__(self,subset,superset):
        self.superset=list(superset)
        if len(set(self.superset))!=len(self.superset):
            raise ValueError('superset must have distinct elements')
        indices=self.subset_indices(subset,self.superset)
        if indices is None or len(set(indices))!=len(indices):
            raise ValueError('subset must contain distinct superset elements')
        self.subset=[self.superset[i] for i in sorted(indices)]

    def __eq__(self,other):
        return isinstance(other,Subset) and (self.subset,self.superset)==(other.subset,other.superset)

    @property
    def size(self):
        return len(self.subset)

    @property
    def superset_size(self):
        return len(self.superset)

    @property
    def cardinality(self):
        return 1<<self.superset_size

    @staticmethod
    def subset_indices(subset,superset):
        try:
            return [list(superset).index(v) for v in subset]
        except ValueError:
            return None

    @staticmethod
    def bitlist_from_subset(subset,superset):
        if any(v not in superset for v in subset):
            raise ValueError('not a subset')
        return ''.join('1' if v in subset else '0' for v in superset)

    @classmethod
    def subset_from_bitlist(cls,super_set,bitlist):
        if len(super_set)!=len(bitlist) or set(bitlist)-{'0','1'}:
            raise ValueError('invalid bitlist')
        return cls([v for v,b in zip(super_set,bitlist) if b=='1'],super_set)

    @property
    def rank_binary(self):
        return int(self.bitlist_from_subset(self.subset,self.superset) or '0',2)

    @property
    def rank_gray(self):
        return _rank(self.bitlist_from_subset(self.subset,self.superset) or '0')

    @property
    def rank_lexicographic(self):
        rank=0;start=0;n=self.superset_size
        for i in self.subset_indices(self.subset,self.superset):
            rank+=1+sum(1<<(n-j-1) for j in range(start,i))
            start=i+1
        return rank

    @classmethod
    def unrank_binary(cls,rank,superset):
        n=len(superset);rank=operator.index(rank)%(1<<n)
        return cls([v for i,v in enumerate(superset) if rank>>(n-i-1)&1],superset)

    @classmethod
    def unrank_gray(cls,rank,superset):
        rank=operator.index(rank)%(1<<len(superset))
        return cls.unrank_binary(rank^(rank>>1),superset)

    def iterate_binary(self,k):
        return self.unrank_binary(self.rank_binary+operator.index(k),self.superset)

    def next_binary(self):return self.iterate_binary(1)
    def prev_binary(self):return self.iterate_binary(-1)

    def iterate_graycode(self,k):
        return self.unrank_gray(self.rank_gray+operator.index(k),self.superset)

    def next_gray(self):return self.iterate_graycode(1)
    def prev_gray(self):return self.iterate_graycode(-1)

    def _unrank_lex(self,rank):
        n=self.superset_size;rank%=1<<n;indices=[];start=0
        while rank:
            rank-=1
            for i in range(start,n):
                count=1<<(n-i-1)
                if rank<count:
                    indices.append(i);start=i+1;break
                rank-=count
        return type(self)([self.superset[i] for i in indices],self.superset)

    def next_lexicographic(self):return self._unrank_lex(self.rank_lexicographic+1)
    def prev_lexicographic(self):return self._unrank_lex(self.rank_lexicographic-1)


def ksubsets(superset,k):
    from itertools import combinations
    return combinations(superset,operator.index(k))
