"""Iterator and container helpers, independent of external symbolic engines."""
from collections import defaultdict
from collections.abc import Mapping
from functools import reduce
import heapq
import itertools
import operator


def sift(seq, keyfunc, binary=False):
    if binary:
        yes,no=[],[]
        for value in seq:
            key=keyfunc(value)
            if key not in (True,False):
                raise ValueError('Binary sift keys must be True or False')
            (yes if key else no).append(value)
        return yes,no
    result=defaultdict(list)
    for value in seq:
        result[keyfunc(value)].append(value)
    return result


def _sequence(value):
    from .. import _core
    if isinstance(value,(str,bytes,Mapping)):
        return False
    if isinstance(value,_core.Expr):
        return getattr(value.func,'__name__','')=='Tuple'
    return hasattr(value,'__getitem__') and hasattr(value,'__iter__')


def flatten(iterable, levels=None, cls=None):
    if levels is not None:
        levels=operator.index(levels)
        if levels<0:
            raise ValueError('Flattening levels must be nonnegative')
        if levels==0:
            return iterable
    result=[]
    active=set()
    def visit(sequence,remaining):
        identity=id(sequence)
        if identity in active:
            raise ValueError('Cyclic flatten input')
        active.add(identity)
        try:
            for value in sequence:
                reducible=(isinstance(value,cls) if cls is not None else
                           _sequence(value) or isinstance(value,(set,frozenset)))
                if reducible and remaining!=0:
                    children=(value.args if hasattr(value,'args') and
                              not getattr(value,'is_NDimArray',False) else value)
                    visit(children,None if remaining is None else remaining-1)
                else:
                    result.append(value)
        finally:
            active.remove(identity)
    visit(iterable,levels)
    return result


def unflatten(iter, n=2):
    try:
        n=operator.index(n)
    except TypeError:
        raise ValueError('Group size must be a positive integer') from None
    values=list(iter)
    if n<=0 or len(values)%n:
        raise ValueError('Input length must be a multiple of the positive group size')
    return [tuple(values[i:i+n]) for i in range(0,len(values),n)]


def reshape(seq, how):
    def width(template):
        total=0
        for item in template:
            try:
                size=operator.index(item)
            except TypeError:
                size=width(item)
            if size<0:
                raise ValueError('Reshape template sizes must be nonnegative')
            total+=size
        return total
    size=width(how)
    values=list(seq)
    if not size or len(values)%size:
        raise ValueError('Input length must be a multiple of the template size')
    cursor=iter(values)
    def build(template):
        result=[]
        for item in template:
            try:
                count=operator.index(item)
            except TypeError:
                result.append(build(item))
            else:
                result.extend(itertools.islice(cursor,count))
        return type(template)(result)
    result=[build(how) for _ in range(len(values)//size)]
    return tuple(result) if isinstance(seq,tuple) else result


def group(seq, multiple=True):
    groups=[list(values) for _,values in itertools.groupby(seq)]
    return groups if multiple else [(values[0],len(values)) for values in groups]


def subsets(seq, k=None, repetition=False):
    operation=itertools.combinations_with_replacement if repetition else itertools.combinations
    if k is not None:
        return operation(seq,operator.index(k))
    seq=tuple(seq)
    return itertools.chain.from_iterable(operation(seq,size) for size in range(len(seq)+1))


def variations(seq, n, repetition=False):
    n=operator.index(n)
    return itertools.product(seq,repeat=n) if repetition else itertools.permutations(seq,n)


cartes=itertools.product


def prefixes(seq):
    for i in range(1,len(seq)+1):
        yield seq[:i]


def postfixes(seq):
    for i in range(1,len(seq)+1):
        yield seq[-i:]


def rotations(s, dir=1):
    seq=list(s)
    if dir not in (1,-1):
        raise ValueError('Rotation direction must be 1 or -1')
    for i in range(len(seq)):
        offset=(dir*i)%len(seq)
        yield seq[offset:]+seq[:offset]


def has_dups(seq):
    seen=set()
    unhashable=[]
    previous_values=[]
    for value in seq:
        try:
            if value in seen or any(value==previous for previous in unhashable):
                return True
            seen.add(value)
        except TypeError:
            if any(value==previous for previous in previous_values):
                return True
            unhashable.append(value)
        previous_values.append(value)
    return False


def has_variety(seq):
    iterator=iter(seq)
    sentinel=object()
    first=next(iterator,sentinel)
    return False if first is sentinel else any(value!=first for value in iterator)


def take(iter, n):
    return [value for _,value in zip(range(operator.index(n)),iter)]


def prod(a, start=1):
    return reduce(operator.mul,a,start)


def dict_merge(*dicts):
    merged={}
    for mapping in dicts:
        merged.update(mapping)
    return merged


def topological_sort(graph, key=None):
    vertices,edges=graph
    vertices=list(dict.fromkeys(vertices))
    positions={value:i for i,value in enumerate(vertices)}
    adjacency={value:set() for value in vertices}
    degree=dict.fromkeys(vertices,0)
    for source,target in edges:
        if source not in adjacency or target not in adjacency:
            raise ValueError('Graph edges must reference declared vertices')
        if target not in adjacency[source]:
            adjacency[source].add(target)
            degree[target]+=1
    key=(lambda value:value) if key is None else key
    ready=[(key(value),positions[value],value) for value in vertices if degree[value]==0]
    heapq.heapify(ready)
    result=[]
    while ready:
        _,_,value=heapq.heappop(ready)
        result.append(value)
        for target in adjacency[value]:
            degree[target]-=1
            if degree[target]==0:
                heapq.heappush(ready,(key(target),positions[target],target))
    if len(result)!=len(vertices):
        raise ValueError('cycle detected')
    return result


__all__=['sift','flatten','unflatten','reshape','group','subsets','variations',
    'cartes','prefixes','postfixes','rotations','has_dups','has_variety','take',
    'prod','dict_merge','topological_sort']

from .sorting import ordered

# Native container algorithms and SymPy-compatible iterator interfaces.
from ._iterable_extensions import (
    is_palindromic, iproduct, multiset, ibin,
    filter_symbols, numbered_symbols, capture, common_prefix,
    common_suffix, strongly_connected_components, connected_components, rotate_left,
    rotate_right, least_rotation, multiset_combinations, multiset_permutations,
    multiset_partitions, partitions, ordered_partitions, binary_partitions,
    uniq, generate_bell, generate_involutions, multiset_derangements,
    generate_derangements, necklaces, bracelets, generate_oriented_forest,
    minlex, runs, sequence_partitions, sequence_partitions_empty,
    kbins, permute_signs, signed_permutations, roundrobin,
    NotIterable, iterable, is_sequence, postorder_traversal,
    default_sort_key, random_derangement, interactive_traversal,
)
__all__ += [
    'is_palindromic', 'iproduct', 'multiset', 'ibin',
    'filter_symbols', 'numbered_symbols', 'capture', 'common_prefix',
    'common_suffix', 'strongly_connected_components', 'connected_components', 'rotate_left',
    'rotate_right', 'least_rotation', 'multiset_combinations', 'multiset_permutations',
    'multiset_partitions', 'partitions', 'ordered_partitions', 'binary_partitions',
    'uniq', 'generate_bell', 'generate_involutions', 'multiset_derangements',
    'generate_derangements', 'necklaces', 'bracelets', 'generate_oriented_forest',
    'minlex', 'runs', 'sequence_partitions', 'sequence_partitions_empty',
    'kbins', 'permute_signs', 'signed_permutations', 'roundrobin',
    'NotIterable', 'iterable', 'is_sequence', 'postorder_traversal',
    'default_sort_key', 'random_derangement', 'interactive_traversal',
]
