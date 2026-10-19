"""Exact incremental stabilizer chains for finite permutation groups.

Permutations act on the right: (a*b)[i] = b[a[i]]. Persistent orbit
transversals and queued Schreier pairs certify the stabilizer at each base
point. No group-element enumeration, floating-point arithmetic or reference
implementation is used to compute order or membership.
"""
from collections import deque
from heapq import heappush,heappop
from math import prod


def multiply(a,b):return tuple(b[i] for i in a)


def inverse(a):
    result=[0]*len(a)
    for i,image in enumerate(a):result[image]=i
    return tuple(result)


class _PairBatch:
    """A lazy rectangle of Schreier pairs over append-only point/generator lists."""
    __slots__=('row','end','column','start_column','end_column')
    def __init__(self,start,end,column,end_column):
        self.row=start;self.end=end
        self.column=self.start_column=column;self.end_column=end_column

    def take(self):
        result=self.row,self.column
        self.column+=1
        if self.column==self.end_column:
            self.row+=1;self.column=self.start_column
        return result


class _WordTracker:
    """Persistent transversal paths and short words in earlier strong generators.

    Signed one-based indices encode generators and their inverses. Paths share
    their prefixes; words reference strong generators rather than recursively
    expanding them into the original generators.
    """
    def __init__(self,strong):
        self.ids={g:i+1 for i,g in enumerate(strong)}
        self.words=[[i+1] for i in range(len(strong))]
        self.paths=[]

    @staticmethod
    def invert(word):return [-term for term in reversed(word)]

    def path(self,level,point):
        path=self.paths[level][point];word=[]
        while path is not None:
            path,term=path;word.append(term)
        word.reverse();return word

    def append(self,g,word):
        reduced=[]
        for term in word:
            if reduced and reduced[-1]==-term:reduced.pop()
            else:reduced.append(term)
        self.ids[g]=len(self.words)+1;self.words.append(reduced)


class StabilizerChain:
    def __init__(self,generators,degree,base=(),track_words=False):
        self.degree=degree;self.identity=tuple(range(degree));self.base=list(base)
        if any(not 0<=i<degree for i in self.base):
            raise ValueError('base points must belong to the permutation domain')
        self.strong=list(dict.fromkeys(g for g in generators if g!=self.identity))
        if any(len(g)!=degree for g in self.strong):raise ValueError('inconsistent generator degrees')
        self._words=_WordTracker(self.strong) if track_words else None
        for g in self.strong:
            if all(g[b]==b for b in self.base):
                self.base.append(next(i for i,image in enumerate(g) if i!=image))
        self.transversals=[];self.inverse_transversals=[];self.orbit_points=[]
        self.level_generators=[];self._pending=[];self._active=[]
        for level,alpha in enumerate(self.base):
            self._new_level(alpha)
            for g in self.strong:
                if all(g[b]==b for b in self.base[:level]):self._add_generator(level,g)
        strong_set=set(self.strong)
        while self._active:
            level=-self._active[0]
            batch=self._pending[level][0];row,column=batch.take()
            if batch.row==batch.end:
                self._pending[level].popleft()
                if not self._pending[level]:heappop(self._active)
            point=self.orbit_points[level][row];g=self.level_generators[level][column]
            representative=self.transversals[level][point]
            schreier=multiply(multiply(representative,g),self.inverse_transversals[level][g[point]])
            if schreier==self.identity:continue
            word=None
            if self._words is not None:
                tracker=self._words
                word=(tracker.path(level,point)+[tracker.ids[g]]+
                      tracker.invert(tracker.path(level,g[point])))
            residue,failed=self.sift(schreier,level+1,word)
            if residue==self.identity:continue
            if failed==len(self.base):
                alpha=next(i for i,image in enumerate(residue) if i!=image)
                self.base.append(alpha);self._new_level(alpha)
            if residue in strong_set:
                raise RuntimeError('stabilizer-chain extension did not make progress')
            if self._words is not None:self._words.append(residue,word)
            strong_set.add(residue);self.strong.append(residue)
            # The residue is already in every group through the originating
            # level. Only intermediate stabilizers require a new generator.
            for index in range(level+1,failed+1):self._add_generator(index,residue)
        self.order=prod(map(len,self.transversals))

    def _new_level(self,alpha):
        self.transversals.append({alpha:self.identity})
        self.inverse_transversals.append({alpha:self.identity})
        self.orbit_points.append([alpha]);self.level_generators.append([])
        if self._words is not None:self._words.paths.append({alpha:None})
        self._pending.append(deque())

    def _enqueue(self,level,start,end,column,end_column):
        if start==end or column==end_column:return
        if not self._pending[level]:heappush(self._active,-level)
        self._pending[level].append(_PairBatch(start,end,column,end_column))

    def _add_generator(self,level,g):
        orbit=self.transversals[level];points=self.orbit_points[level]
        old_count=len(points);generators=self.level_generators[level]
        column=len(generators);generators.append(g)
        def discover(point,step):
            image=step[point]
            if image not in orbit:
                representative=multiply(orbit[point],step)
                orbit[image]=representative
                self.inverse_transversals[level][image]=inverse(representative)
                points.append(image)
                if self._words is not None:
                    paths=self._words.paths[level]
                    paths[image]=(paths[point],self._words.ids[step])
        for row in range(old_count):discover(points[row],g)
        row=old_count
        while row<len(points):
            for step in generators:discover(points[row],step)
            row+=1
        self._enqueue(level,0,old_count,column,column+1)
        self._enqueue(level,old_count,len(points),0,len(generators))

    def sift(self,g,start=0,word=None):
        for level in range(start,len(self.base)):
            if g==self.identity:break
            representative=self.inverse_transversals[level].get(g[self.base[level]])
            if representative is None:return g,level
            if word is not None:
                word.extend(self._words.invert(self._words.path(level,g[self.base[level]])))
            g=multiply(g,representative)
        return g,len(self.base)

    def contains(self,g):return self.sift(g)[0]==self.identity

    def generate(self):
        # Mixed-radix traversal caches partial compositions. Trivial orbit
        # factors are omitted, giving O(degree) amortized work per output.
        choices=[list(level.values()) for level in self.transversals if len(level)>1]
        if not choices:
            yield self.identity
            return
        indices=[0]*len(choices);partial=[self.identity]*(len(choices)+1)
        for level,values in enumerate(choices):partial[level+1]=multiply(values[0],partial[level])
        while True:
            yield partial[-1]
            level=len(choices)-1
            while level>=0:
                indices[level]+=1
                if indices[level]<len(choices[level]):break
                indices[level]=0;level-=1
            if level<0:return
            for index in range(level,len(choices)):
                partial[index+1]=multiply(choices[index][indices[index]],partial[index])
