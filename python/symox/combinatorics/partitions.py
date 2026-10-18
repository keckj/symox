"""Set partitions and integer partitions with exact ranking."""
from functools import lru_cache
import operator
from collections import Counter


@lru_cache(maxsize=8192)
def _completions(remaining,blocks):
    if not remaining:return 1
    return blocks*_completions(remaining-1,blocks)+_completions(remaining-1,blocks+1)


def _rgs(values):
    values=[operator.index(v) for v in values]
    maximum=-1
    for v in values:
        if v<0 or v>maximum+1:raise ValueError('invalid restricted growth string')
        maximum=max(maximum,v)
    return values


def RGS_enum(m):
    m=operator.index(m)
    if m<0:raise ValueError('negative string length')
    if m>256:raise NotImplementedError('restricted growth strings are limited to length 256')
    return _completions(m-1,1) if m else 0


def RGS_rank(rgs):
    values=_rgs(rgs);rank=0;blocks=1
    if len(values)>256:raise NotImplementedError('restricted growth strings are limited to length 256')
    for i,v in enumerate(values[1:],1):
        rank+=v*_completions(len(values)-i-1,blocks)
        blocks=max(blocks,v+1)
    return rank


def RGS_unrank(rank,m):
    m=operator.index(m);rank=operator.index(rank)
    if m<1 or not 0<=rank<RGS_enum(m):raise ValueError('restricted growth rank out of range')
    values=[0];blocks=1
    for i in range(1,m):
        count=_completions(m-i-1,blocks)
        v=min(rank//count,blocks)
        rank-=v*count;values.append(v);blocks=max(blocks,v+1)
    return values


class Partition:
    def __init__(self,*partition):
        from .. import default_sort_key
        blocks=[frozenset(block) for block in partition]
        if any(not b for b in blocks) or sum(map(len,blocks))!=len(set().union(*blocks)):
            raise ValueError('partition blocks must be nonempty and disjoint')
        self.members=tuple(sorted(set().union(*blocks),key=default_sort_key))
        self.size=len(self.members)
        self._blocks=tuple(sorted(blocks,key=lambda b:min(self.members.index(v) for v in b)))

    @property
    def partition(self):
        return [[v for v in self.members if v in block] for block in self._blocks]

    @property
    def RGS(self):
        return tuple(next(i for i,b in enumerate(self._blocks) if v in b) for v in self.members)

    @property
    def rank(self):return RGS_rank(self.RGS)

    @classmethod
    def from_rgs(cls,rgs,elements):
        values=_rgs(rgs);elements=list(elements)
        if len(values)!=len(elements) or len(set(elements))!=len(elements):raise ValueError('invalid partition elements')
        return cls(*([v for v,r in zip(elements,values) if r==i] for i in range(max(values,default=-1)+1)))

    def __add__(self,k):
        if not self.size:return self
        return self.from_rgs(RGS_unrank((self.rank+operator.index(k))%RGS_enum(self.size),self.size),self.members)

    def __sub__(self,k):return self+(-operator.index(k))
    def __eq__(self,other):return isinstance(other,Partition) and self._blocks==other._blocks
    def __hash__(self):return hash(frozenset(self._blocks))
    def __lt__(self,other):return self.RGS<other.RGS
    def __le__(self,other):return self==other or self<other
    def __repr__(self):return 'Partition(%s)'%', '.join(map(str,self.partition))


class IntegerPartition:
    def __init__(self,partition,integer=None):
        if integer is not None:partition,integer=integer,operator.index(partition)
        if isinstance(partition,dict):
            parts=[]
            for k,v in partition.items():
                k,v=operator.index(k),operator.index(v)
                if v<0:raise ValueError('negative multiplicity')
                parts.extend([k]*v)
        else:parts=[operator.index(v) for v in partition]
        self.partition=sorted(parts,reverse=True)
        self.integer=sum(parts)
        if any(v<1 for v in parts) or integer is not None and self.integer!=integer:
            raise ValueError('invalid integer partition')

    def as_dict(self):return dict(Counter(self.partition))
    @property
    def conjugate(self):return [sum(v>=i for v in self.partition) for i in range(1,max(self.partition,default=0)+1)]
    def as_ferrers(self,char='#'):return '\n'.join(char*v for v in self.partition)
    def __eq__(self,other):return isinstance(other,IntegerPartition) and self.partition==other.partition
    def __lt__(self,other):return self.partition[::-1]<other.partition[::-1]
    def __le__(self,other):return self==other or self<other
    def __repr__(self):return 'IntegerPartition(%r)'%self.partition
    def __str__(self):return str(self.partition)

    def next_lex(self):
        parts=list(self.partition)
        if not parts:return self
        if len(parts)==1:return type(self)([1]*self.integer)
        # Ascending lexicographic order of the decreasing parts.
        tail=parts.pop();j=len(parts)-1
        while j>=0 and (j>0 and parts[j]==parts[j-1]):
            tail+=parts.pop();j-=1
        parts[j]+=1;tail-=1
        parts.extend([1]*tail)
        return type(self)(parts)

    def prev_lex(self):
        parts=list(self.partition)
        if not parts:return self
        if parts[0]==1:return type(self)([self.integer])
        tail=0
        while parts[-1]==1:tail+=parts.pop()
        value=parts.pop()-1;tail+=value+1
        parts.extend([value]*(tail//value))
        if tail%value:parts.append(tail%value)
        return type(self)(parts)


def RGS_generalized(m):
    from symox import zeros
    m=operator.index(m)
    if m<0:raise ValueError('growth table order must be nonnegative')
    if m>256:raise ValueError('growth table is limited to order 256')
    table=zeros(m+1)
    for j in range(m+1):table[0,j]=1
    for i in range(1,m+1):
        for j in range(m-i+1):table[i,j]=j*table[i-1,j]+table[i-1,j+1]
    return table


def random_integer_partition(n,seed=None):
    from random import Random
    from ..core.random import rng
    n=operator.index(n)
    if n<1:raise ValueError('partition size must be positive')
    if isinstance(seed,(list,tuple)):
        sequence=iter(seed)
        def randint(a,b):
            while True:
                value=operator.index(next(sequence))
                if a<=value<=b:return value
    else:randint=rng.randint if seed is None else Random(seed).randint
    result=[]
    while n:
        part=randint(1,n);count=randint(1,n//part);result.extend([part]*count);n-=part*count
    return sorted(result,reverse=True)
