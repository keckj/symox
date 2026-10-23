"""Bijection between labelled trees and Prüfer words."""
import heapq
import operator


class Prufer:
    def __init__(self,*args,**kwargs):
        if kwargs or not 1<=len(args)<=2:
            raise TypeError('expected a Prüfer sequence or edges and node count')
        values=list(args[0])
        if values and isinstance(values[0],(list,tuple)):
            n=operator.index(args[1]) if len(args)==2 else max(max(e) for e in values)+1
            self._code=self.to_prufer(values,n)
        else:
            self._code=[operator.index(v) for v in values]
            n=len(self._code)+2
            if any(not 0<=v<n for v in self._code):
                raise ValueError('Prüfer entries must be node labels')
        self.nodes=n
        # The reference prints its constructor arguments: Prufer([1, 2, 0]).
        self._args=(list(self._code) if not values or not isinstance(values[0],(list,tuple)) else [list(e) for e in values],)+tuple(args[1:])

    def __repr__(self):
        return 'Prufer(%s)'%', '.join(map(str,self._args))

    __str__=__repr__

    @property
    def prufer_repr(self):return list(self._code)
    @property
    def tree_repr(self):return self.to_tree(self._code)
    @property
    def size(self):return self.nodes**(self.nodes-2)
    @property
    def rank(self):return self.prufer_rank()

    def prufer_rank(self):
        result=0
        for v in self._code:
            result=result*self.nodes+v
        return result

    @staticmethod
    def to_tree(prufer):
        code=[operator.index(v) for v in prufer];n=len(code)+2
        if any(not 0<=v<n for v in code):
            raise ValueError('invalid Prüfer label')
        degrees=[1]*n
        for v in code:degrees[v]+=1
        leaves=[v for v in range(n) if degrees[v]==1];heapq.heapify(leaves)
        edges=[]
        for v in code:
            leaf=heapq.heappop(leaves);edges.append(sorted([leaf,v]))
            degrees[v]-=1
            if degrees[v]==1:heapq.heappush(leaves,v)
        edges.append(sorted(leaves))
        return edges

    @staticmethod
    def to_prufer(tree,n):
        n=operator.index(n)
        if n<2:
            raise ValueError('a tree needs at least two vertices')
        edges=[tuple(operator.index(v) for v in e) for e in tree]
        if len(edges)!=n-1 or any(len(e)!=2 or e[0]==e[1] or any(not 0<=v<n for v in e) for e in edges):
            raise ValueError('invalid tree edges')
        adjacent=[set() for _ in range(n)]
        for a,b in edges:adjacent[a].add(b);adjacent[b].add(a)
        if sum(map(len,adjacent))!=2*(n-1):raise ValueError('duplicate edges')
        leaves=[v for v in range(n) if len(adjacent[v])==1];heapq.heapify(leaves)
        code=[]
        for _ in range(n-2):
            if not leaves:raise ValueError('edges do not form a tree')
            leaf=heapq.heappop(leaves);parent=next(iter(adjacent[leaf]));code.append(parent)
            adjacent[leaf].clear();adjacent[parent].remove(leaf)
            if len(adjacent[parent])==1:heapq.heappush(leaves,parent)
        if len(leaves)!=2 or leaves[1] not in adjacent[leaves[0]]:
            raise ValueError('edges do not form a connected tree')
        return code

    @staticmethod
    def edges(*runs):
        edges={tuple(sorted((operator.index(a),operator.index(b)))) for run in runs for a,b in zip(run,run[1:])}
        if not edges:raise ValueError('empty tree')
        labels={v for e in edges for v in e};lo=min(labels);hi=max(labels)
        if labels!=set(range(lo,hi+1)):raise ValueError('missing node labels')
        return sorted([a-lo,b-lo] for a,b in edges),hi-lo+1

    @classmethod
    def unrank(cls,rank,n):
        n=operator.index(n);rank=operator.index(rank)
        if n<2 or not 0<=rank<n**(n-2):raise ValueError('tree rank out of range')
        code=[0]*(n-2)
        for i in range(n-3,-1,-1):rank,code[i]=divmod(rank,n)
        return cls(code)

    def next(self,delta=1):return self.unrank((self.rank+operator.index(delta))%self.size,self.nodes)
    def prev(self,delta=1):return self.next(-operator.index(delta))
