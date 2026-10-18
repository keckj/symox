"""Exact signed tensor canonicalization by bounded permutation orbits.

The reference uses specialized stabilizer algorithms. This implementation uses
explicit finite orbits, rejects orbits exceeding 100000 states, and detects
opposite signs of the same tensor exactly.
"""
from collections import deque
from operator import index
from .permutations import Permutation
from .perm_groups import PermutationGroup
MAX_ORBIT=100000


def _array(value):return value.array_form if isinstance(value,Permutation) else list(value)
def _identity(n):return Permutation(range(n))
def _group(gens):
    values=[g if isinstance(g,Permutation) else Permutation(g) for g in gens]
    return PermutationGroup(values)
def _rmul(a,b):return [a[i] for i in b]


def dummy_sgs(dummies,sym,n):
    n=index(n);dummies=list(map(index,dummies))
    if len(dummies)%2 or len(set(dummies))!=len(dummies) or any(i<0 or i>=n for i in dummies):raise ValueError('dummy indices require disjoint pairs inside the tensor slots')
    if sym not in (0,1,None):raise ValueError('metric symmetry must be 0, 1 or None')
    pairs=list(zip(dummies[::2],dummies[1::2]));result=[]
    if sym is not None:
        for a,b in pairs:
            values=list(range(n+2));values[a],values[b]=values[b],values[a]
            if sym==1:values[-2:]=[n+1,n]
            result.append(values)
    for (a,b),(c,d) in zip(pairs,pairs[1:]):
        values=list(range(n+2));values[a],values[b],values[c],values[d]=c,d,a,b;result.append(values)
    return result


def perm_af_direct_product(gens1,gens2,signed=True):
    first,second=[_array(g) for g in gens1],[_array(g) for g in gens2]
    tail=2 if signed else 0;n1=len(first[0])-tail;n2=len(second[0])-tail;size=n1+n2+tail;result=[]
    for gens,offset,width in [(first,0,n1),(second,n1,n2)]:
        for generator in gens:
            values=list(range(size))
            for i in range(width):values[offset+i]=offset+generator[i]
            if signed and generator[-1]!=width+1:values[-2:]=[size-1,size-2]
            result.append(values)
    return result


def bsgs_direct_product(base1,gens1,base2,gens2,signed=True):
    shift=len(_array(gens1[0]))-(2 if signed else 0)
    arrays=perm_af_direct_product(gens1,gens2,signed)
    unique=list(dict.fromkeys(tuple(a) for a in arrays if a!=list(range(len(a)))))
    return list(base1)+[shift+i for i in base2],[Permutation(a) for a in unique] or [_identity(len(arrays[0]))]


def get_symmetric_group_sgs(n,antisym=False):
    n=index(n)
    if n<1:raise ValueError('tensor symmetry rank must be positive')
    generators=[]
    for i in range(n-1):
        values=list(range(n+2));values[i],values[i+1]=values[i+1],values[i]
        if antisym:values[-2:]=[n+1,n]
        generators.append(Permutation(values))
    return list(range(n-1)),generators or [_identity(n+2)]


def get_transversals(base,gens):
    elements=list(_group(gens)._enumerate());result=[]
    for point in base:
        representatives={}
        for element in elements:representatives.setdefault(element(point),element.array_form)
        result.append(representatives)
        elements=[p for p in elements if p(point)==point]
    return result


def get_minimal_bsgs(base,gens):
    elements=list(_group(gens)._enumerate());size=elements[0].size;result=[];strong=set()
    while len(elements)>1:
        point=next(i for i in range(size) if any(p(i)!=i for p in elements))
        result.append(point)
        selected=[];generated={_identity(size)}
        for element in elements:
            if element not in generated:
                selected.append(element);generated=set(PermutationGroup(selected)._enumerate())
                if len(generated)==len(elements):break
        strong.update(selected);elements=[p for p in elements if p(point)==point]
    strong.discard(_identity(size))
    return result,sorted(strong,key=lambda p:p.array_form) or [_identity(size)]


def transversal2coset(size,base,transversal):
    lookup=dict(zip(base,transversal));identity=list(range(size))
    result=[sorted(lookup[i].values()) if i in lookup else [identity[:]] for i in range(size)]
    while result and result[-1]==[identity]:result.pop()
    return result


def tensor_gens(base,gens,list_free_indices,sym=0):
    if sym not in (0,1,None):raise ValueError('tensor exchange symmetry must be 0, 1 or None')
    size=len(_array(gens[0]));rank=size-2;count=len(list_free_indices);total=count*rank+2;result=[]
    if count<1:return 2,[],[_identity(2)]
    group=_group(gens)
    for block,free in enumerate(list_free_indices):
        if any(i<0 or i>=rank for i in free):raise ValueError('fixed tensor slot is out of range')
        candidates=[p for p in group._enumerate() if all(p(i)==i for i in free)] if free else group.generators
        for permutation in candidates:
            if permutation==group.identity:continue
            values=list(range(total))
            for i in range(rank):values[block*rank+i]=block*rank+permutation(i)
            if permutation(size-1)!=size-1:values[-2:]=[total-1,total-2]
            result.append(Permutation(values))
    movable=[block for block,free in enumerate(list_free_indices) if not free]
    if sym is not None:
        for a,b in zip(movable,movable[1:]):
            values=list(range(total))
            for i in range(rank):values[a*rank+i],values[b*rank+i]=b*rank+i,a*rank+i
            if sym==1:values[-2:]=[total-1,total-2]
            result.append(Permutation(values))
    if not result:return total,[],[_identity(total)]
    base_,strong=get_minimal_bsgs([],result)
    return total,base_,strong


def gens_products(*v):
    if not v:return 2,[],[_identity(2)]
    size,base,gens=tensor_gens(*v[0])
    for specification in v[1:]:
        size2,base2,gens2=tensor_gens(*specification)
        base,gens=bsgs_direct_product(base,gens,base2,gens2)
        size+=size2-2
    return size,base,gens


def canonical_free(base,gens,g,num_free):
    values=_array(g);num_free=index(num_free)
    if not 0<=num_free<=len(values)-2:raise ValueError('invalid free index count')
    if not num_free or not base:return values
    for transversal in get_transversals(base,gens):
        possibilities=[_rmul(values,permutation) for permutation in transversal.values()]
        values=min(possibilities,key=lambda row:tuple(row.index(i) for i in range(num_free)))
    return values


def _canonical_orbit(g,slot_gens,dummy_gens):
    initial=tuple(_array(g));size=len(initial);slots=[_array(p) for p in slot_gens];dummies=[_array(p) for p in dummy_gens]
    if any(len(p)!=size for p in slots+dummies):raise ValueError('tensor and symmetry permutation sizes differ')
    queue=deque([initial]);seen={initial};signs={initial[:-2]:initial[-1]};minimum=initial
    while queue:
        current=queue.popleft()
        candidates=[tuple(current[i] for i in permutation) for permutation in slots]
        candidates.extend(tuple(permutation[i] for i in current) for permutation in dummies)
        for candidate in candidates:
            if candidate in seen:continue
            if candidate[:-2] in signs and signs[candidate[:-2]]!=candidate[-1]:return 0
            if len(seen)>=MAX_ORBIT:raise NotImplementedError('tensor canonicalization orbit exceeds 100000 states')
            signs[candidate[:-2]]=candidate[-1];seen.add(candidate);queue.append(candidate);minimum=min(minimum,candidate)
    return list(minimum)


def double_coset_can_rep(dummies,sym,b_S,sgens,S_transversals,g):
    if sym in (0,1,None):dummies,sym=[dummies],[sym]
    if len(dummies)!=len(sym):raise ValueError('one metric symmetry is required for each dummy type')
    size=len(_array(g));dg=[generator for group,metric in zip(dummies,sym) for generator in dummy_sgs(group,metric,size-2)]
    return _canonical_orbit(g,sgens,dg)


def canonicalize(g,dummies,msym,*v):
    specifications=[]
    for base,gens,count,exchange in v:
        count=index(count)
        if count<0:raise ValueError('tensor multiplicities must be nonnegative')
        specifications.append((base,gens,[[] for _ in range(count)],exchange))
    size,base,gens=gens_products(*specifications)
    if size!=len(_array(g)):raise ValueError('symmetry ranks do not match the tensor permutation')
    return double_coset_can_rep(dummies,msym,base,gens,[],g)


riemann_bsgs=([0,2],[Permutation([1,0,2,3,5,4]),Permutation([0,1,3,2,5,4]),Permutation([2,3,0,1,4,5])])
__all__=['dummy_sgs','perm_af_direct_product','bsgs_direct_product','get_symmetric_group_sgs','get_transversals',
 'get_minimal_bsgs','transversal2coset','tensor_gens','gens_products','canonical_free','double_coset_can_rep','canonicalize']
