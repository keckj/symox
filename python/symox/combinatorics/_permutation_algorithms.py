"""Exact finite-permutation algorithms over arbitrary Python integers.

Order-statistic trees avoid repeated array searches in inversion vectors,
lexicographic and Trotter--Johnson ranking. No symbolic reference runtime is
used here. Rank arithmetic retains the full factorial-sized integer.
"""
import math
import operator
import random as _random
from collections import deque


class _Fenwick:
    def __init__(self,size,ones=False):
        self.tree=[0]+([i & -i for i in range(1,size+1)] if ones else [0]*size)

    def add(self,index,delta):
        index+=1
        while index<len(self.tree):
            self.tree[index]+=delta
            index+=index & -index

    def prefix(self,end):
        result=0
        while end:
            result+=self.tree[end]
            end-=end & -end
        return result

    def select(self,rank):
        """Index of a zero-based order statistic among nonnegative counts."""
        index=0
        bit=1 << ((len(self.tree)-1).bit_length()-1)
        while bit:
            candidate=index+bit
            if candidate<len(self.tree) and self.tree[candidate]<=rank:
                index=candidate
                rank-=self.tree[candidate]
            bit>>=1
        return index


class _Algorithms:
    is_Permutation=True
    default_assumptions={}
    print_cyclic=None

    def __iter__(self):return iter(self._images)
    def __rxor__(self,i):return self(i)
    def __add__(self,other):return type(self).unrank_lex(self.size,self.rank()+operator.index(other))
    def __sub__(self,other):return self+-operator.index(other)
    def atoms(self):return {self}
    def copy(self):return type(self)(self._images)
    def support(self):return [i for i,value in enumerate(self._images) if i!=value]
    def length(self):return sum(i!=value for i,value in enumerate(self._images))
    def parity(self):
        if getattr(self,'_parity_cache',None) is None:
            self._parity_cache=sum(len(cycle)-1 for cycle in self.cyclic_form)%2
        return self._parity_cache
    def signature(self):return -1 if self.parity() else 1
    is_even=property(lambda self:not self.parity())
    is_odd=property(lambda self:bool(self.parity()))
    is_identity=property(lambda self:self.is_Identity)
    is_Empty=property(lambda self:self.size==0)
    is_Singleton=property(lambda self:self.size==1)
    cardinality=property(lambda self:math.factorial(self.size))
    cycles=property(lambda self:len(self.full_cyclic_form))

    @property
    def full_cyclic_form(self):
        visited=[False]*self.size
        result=[]
        for start in range(self.size):
            if visited[start]:continue
            cycle=[];value=start
            while not visited[value]:
                visited[value]=True;cycle.append(value);value=self._images[value]
            result.append(cycle)
        return result

    @property
    def cycle_structure(self):
        if getattr(self,'_cycle_structure_cache',None) is not None:return dict(self._cycle_structure_cache)
        counts={};fixed=self.size
        for cycle in self.cyclic_form:
            counts[len(cycle)]=counts.get(len(cycle),0)+1
            fixed-=len(cycle)
        if fixed:counts[1]=fixed
        self._cycle_structure_cache=counts
        return dict(counts)

    def order(self):return math.lcm(*self.cycle_structure)

    def max(self):
        return max((value for i,value in enumerate(self._images) if i!=value)) if self.size else 0

    def min(self):
        return min((value for i,value in enumerate(self._images) if i!=value)) if self.size else 0

    def ascents(self):return [i for i in range(self.size-1) if self._images[i]<self._images[i+1]]
    def descents(self):return [i for i in range(self.size-1) if self._images[i]>self._images[i+1]]
    def index(self):return sum(self.descents())

    def runs(self):
        result=[]
        for value in self._images:
            if not result or result[-1][-1]>value:result.append([])
            result[-1].append(value)
        return result

    def transpositions(self):
        return [(cycle[0],value) for cycle in self.cyclic_form for value in reversed(cycle[1:])]

    def list(self,size=None):
        result=self.array_form
        if size is not None:
            size=operator.index(size)
            if size>self.size:result.extend(range(self.size,size))
            else:
                while result and result[-1]==len(result)-1:result.pop()
        return result

    def resize(self,n):
        n=operator.index(n)
        if n>=self.size:
            return self if n==self.size else type(self)(self.array_form+list(range(self.size,n)))
        if n<0:raise ValueError('permutation size must be nonnegative')
        if any((i<n)!=(value<n) for i,value in enumerate(self._images)):
            raise ValueError('resizing would break a permutation cycle')
        return type(self)(self._images[:n])

    def inversion_vector(self):
        tree=_Fenwick(self.size)
        result=[0]*self.size
        for i in range(self.size-1,-1,-1):
            value=self._images[i]
            result[i]=tree.prefix(value)
            tree.add(value,1)
        return result[:-1]

    def inversions(self):return sum(self.inversion_vector())

    @classmethod
    def from_inversion_vector(cls,inversion):
        entries=list(inversion);degree=len(entries)+1
        tree=_Fenwick(degree,ones=True);images=[]
        for i,entry in enumerate(entries):
            entry=operator.index(entry);remaining=degree-i
            if entry<0:entry+=remaining
            if not 0<=entry<remaining:raise ValueError('invalid inversion vector')
            value=tree.select(entry);images.append(value);tree.add(value,-1)
        images.append(tree.select(0))
        return cls(images)

    def rank(self):
        if getattr(self,'_rank',None) is not None:return self._rank
        result=0
        for i,value in enumerate(self.inversion_vector()):result=result*(self.size-i)+value
        self._rank=result
        return result

    @classmethod
    def unrank_lex(cls,size,rank):
        size=operator.index(size);rank=operator.index(rank)
        if size<0:raise ValueError('permutation size must be nonnegative')
        rank%=math.factorial(size)
        digits=[0]*size
        for radix in range(1,size+1):rank,digits[size-radix]=divmod(rank,radix)
        tree=_Fenwick(size,ones=True);images=[]
        for digit in digits:
            value=tree.select(digit);images.append(value);tree.add(value,-1)
        return cls(images)

    def next_lex(self):
        images=self.array_form
        pivot=self.size-2
        while pivot>=0 and images[pivot]>images[pivot+1]:pivot-=1
        if pivot<0:return None
        successor=self.size-1
        while images[successor]<images[pivot]:successor-=1
        images[pivot],images[successor]=images[successor],images[pivot]
        images[pivot+1:]=reversed(images[pivot+1:])
        return type(self)(images)

    def rank_nonlex(self,inv_perm=None):
        inverse=(~self).array_form if inv_perm is None else inv_perm
        if not inverse:return 0
        images=self.array_form;digits=[]
        for degree in range(self.size,1,-1):
            value=images[degree-1];position=inverse[degree-1]
            images[degree-1],images[position]=images[position],value
            inverse[degree-1],inverse[value]=inverse[value],position
            digits.append((degree,value))
        rank=0
        for degree,value in reversed(digits):rank=value+degree*rank
        return rank

    @classmethod
    def unrank_nonlex(cls,n,r):
        n=operator.index(n);r=operator.index(r)
        if n<0:raise ValueError('permutation size must be nonnegative')
        r%=math.factorial(n)
        images=list(range(n))
        # Digits depend only on r modulo n!, including for negative r.
        for degree in range(n,0,-1):
            r,index=divmod(r,degree)
            images[degree-1],images[index]=images[index],images[degree-1]
        return cls(images)

    def next_nonlex(self):
        rank=self.rank_nonlex()+1
        return None if rank==self.cardinality else type(self).unrank_nonlex(self.size,rank)

    def rank_trotterjohnson(self):
        positions=(~self)._images;tree=_Fenwick(self.size);rank=0
        for value,position in enumerate(positions):
            before=tree.prefix(position)
            rank=(value+1)*rank+(value-before if rank%2==0 else before)
            tree.add(position,1)
        return rank

    @classmethod
    def unrank_trotterjohnson(cls,size,rank):
        size=operator.index(size);rank=operator.index(rank)
        if size<0:raise ValueError('permutation size must be nonnegative')
        # Out-of-range reference inputs can produce duplicate images. Extend
        # the ordering periodically, always returning an actual permutation.
        rank%=math.factorial(size)
        positions=[0]*size
        for degree in range(size,1,-1):
            rank,digit=divmod(rank,degree)
            positions[degree-1]=degree-1-digit if rank%2==0 else digit
        tree=_Fenwick(size,ones=True);images=[0]*size
        for value in range(size-1,-1,-1):
            position=tree.select(positions[value]);images[position]=value;tree.add(position,-1)
        return cls(images)

    def next_trotterjohnson(self):
        # A value moves left when its lower-valued subpermutation is even,
        # right when it is odd. Move the largest mobile value one position.
        positions=(~self)._images;tree=_Fenwick(self.size);directions=[];parity=0
        for value,position in enumerate(positions):
            directions.append(1 if parity else -1)
            parity^=(value-tree.prefix(position))%2
            tree.add(position,1)
        for value in range(self.size-1,-1,-1):
            position=positions[value];neighbor=position+directions[value]
            if 0<=neighbor<self.size and self._images[neighbor]<value:
                images=self.array_form
                images[position],images[neighbor]=images[neighbor],images[position]
                return type(self)(images)
        return None

    def _same_degree(self,other):
        if self.size!=other.size:raise ValueError('permutations must have equal sizes')

    def commutes_with(self,other):return self*other==other*self

    def commutator(self,x):
        self._same_degree(x)
        return ~x * ~self * x * self

    def __xor__(self,other):
        if not isinstance(other,type(self)):return NotImplemented
        self._same_degree(other)
        return ~other * self * other

    def mul_inv(self,other):return other * ~self

    @classmethod
    def rmul(cls,*args):
        result=args[0]
        for value in args[1:]:result=value*result
        return result

    @classmethod
    def rmul_with_af(cls,*args):return cls.rmul(*args)

    @classmethod
    def from_sequence(cls,i,key=None):
        sequence=sorted(enumerate(i),key=(lambda pair:key(pair[1])) if key else lambda pair:(pair[1],pair[0]))
        return ~cls([index for index,_ in sequence])

    @classmethod
    def random(cls,n):
        images=list(range(n));_random.shuffle(images);return cls(images)

    @classmethod
    def josephus(cls,m,n,s=1):
        m=max(operator.index(m),1);n=max(operator.index(n),0)
        threshold=max(s,1)
        if n<=threshold:return cls(range(n))
        if m<=max(n.bit_length(),1):
            queue=deque(range(n));images=[]
            while len(queue)>threshold:
                queue.rotate(-(m-1));images.append(queue.popleft())
            images.extend(queue)
            return cls(images)
        tree=_Fenwick(n,ones=True);remaining=n;index=0;images=[]
        alive=bytearray([1])*n
        while remaining>threshold:
            index=(index+m-1)%remaining
            value=tree.select(index);images.append(value);tree.add(value,-1);remaining-=1;alive[value]=0
        if remaining:
            survivors=[i for i in range(n) if alive[i]]
            index%=remaining
            images.extend(survivors[index:]);images.extend(survivors[:index])
        return cls(images)

    def get_adjacency_matrix(self):
        from .. import zeros
        matrix=zeros(self.size)
        for a,b in zip(self._images,self._images[1:]):matrix[a,b]=1
        return matrix

    def get_precedence_matrix(self):
        from .. import zeros
        matrix=zeros(self.size)
        for i,a in enumerate(self._images):
            for b in self._images[i+1:]:matrix[a,b]=1
        return matrix

    def get_adjacency_distance(self,other):
        self._same_degree(other)
        pairs=set(zip(self._images,self._images[1:]))
        shared=sum(pair in pairs for pair in zip(other._images,other._images[1:]))
        return self.size-1-shared

    def get_positional_distance(self,other):
        self._same_degree(other)
        return sum(abs(a-b) for a,b in zip(self._images,other._images))

    def get_precedence_distance(self,other):
        self._same_degree(other)
        inverse=(~other)._images
        return type(self)([inverse[value] for value in self._images]).inversions()


def install(permutation):
    for name,value in vars(_Algorithms).items():
        if name not in {'__dict__','__weakref__','__doc__','__module__'}:
            setattr(permutation,name,value)
