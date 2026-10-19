"""Exact signed-orbit tensor canonicalization and metric elimination.

The public canon_bp interface uses bounded orbit enumeration here, rather
than pretending to implement the large-group Butler--Portugal machinery.
"""
from ._indexing import TensorIndex,TensorSymmetry,_IndexStructure,_index_key,_name_key,MAX_ORBIT
from ._expressions import Tensor,TensExpr,TensAdd,TensMul,_split,_product,_normalize


def _state(factors):
    factors=tuple(_normalize(factors))
    return factors,tuple((factor.component,factor.indices) for factor in factors)


def _lex(factors):
    structure=_IndexStructure([i for f in factors for i in f.indices])
    dummies={_name_key(structure.indices[upper]) for upper,lower in structure.dum}
    return tuple(((factor.component.name,repr(factor.component.args[1:])),
                  tuple((_name_key(index) in dummies,*_index_key(index)) for index in factor.indices))
                 for factor in factors)


def canon_bp(p):
    from .. import S
    if not isinstance(p,TensExpr):return p
    if p.is_canon_bp:return p
    if isinstance(p,TensAdd):return TensAdd(*(canon_bp(arg) for arg in p.args))
    if not isinstance(p,(Tensor,TensMul)):
        raise NotImplementedError('Canonicalization requires tensors, their sums or their products')
    coefficient,factors=_split(p)
    factors,key=_state(factors);seen={key:1};queue=[(factors,1)]
    best=factors;best_key=_lex(factors);best_sign=1
    for state,sign in queue:
        moves=[]
        for position,factor in enumerate(state):
            rank=len(factor.indices)
            for generator in factor.component.symmetry.generators:
                images=generator.array_form
                if images==list(range(rank+2)):continue
                indices=tuple(factor.indices[i] for i in images[:rank])
                replacement=Tensor._raw(factor.component,indices,type(factor))
                moves.append((state[:position]+(replacement,)+state[position+1:],-1 if images[-1]!=rank+1 else 1))
        for position in range(len(state)-1):
            comm=state[position].component.commutes_with(state[position+1].component)
            if comm is not None:
                moves.append((state[:position]+(state[position+1],state[position])+state[position+2:],-1 if comm else 1))
        structure=_IndexStructure([i for f in state for i in f.indices])
        for upper,lower in structure.dum:
            index=structure.indices[upper];symmetry=index.tensor_index_type.args[4]
            if symmetry not in (1,-1):continue
            changed=[Tensor._raw(f.component,tuple(-i if _name_key(i)==_name_key(index) else i for i in f.indices),type(f)) for f in state]
            moves.append((changed,symmetry))
        for changed,scale in moves:
            changed,key=_state(changed);new_sign=sign*scale
            if key in seen:
                if seen[key]!=new_sign:return S.Zero
                continue
            if len(seen)>=MAX_ORBIT:raise NotImplementedError('Tensor symmetry orbit exceeds 100000 states')
            seen[key]=new_sign;queue.append((changed,new_sign))
            lex=_lex(changed)
            if lex<best_key:best,best_key,best_sign=changed,lex,new_sign
    return _product(coefficient*best_sign,list(best),canonical=True)


def contract_metric(t,g):
    if not isinstance(t,TensExpr):return t
    if isinstance(t,TensAdd):return TensAdd(*(contract_metric(term,g) for term in t.args))
    if g.rank!=2 or g.index_types[0]!=g.index_types[1]:raise ValueError('Metric contraction requires a rank-two head of one index type')
    if g.symmetry not in (TensorSymmetry.fully_symmetric(2),TensorSymmetry.fully_symmetric(-2)):
        raise NotImplementedError('Automatic metric contraction requires a symmetric or antisymmetric metric')
    coefficient,factors=_split(t)
    antisymmetric=g.symmetry==TensorSymmetry.fully_symmetric(-2)
    while True:
        changed=False
        for position,factor in enumerate(factors):
            if factor.component!=g:continue
            left,right=factor.indices
            if _name_key(left)==_name_key(right):
                coefficient*=g.index_types[0].dim
                if antisymmetric and left.is_up:coefficient=-coefficient
                factors.pop(position);changed=True;break
            for slot,index in enumerate(factor.indices):
                for other_pos,other in enumerate(factors):
                    if other_pos==position:continue
                    for other_slot,other_index in enumerate(other.indices):
                        if other_index!=-index:continue
                        replacement=list(other.indices);replacement[other_slot]=factor.indices[1-slot]
                        factors[other_pos]=Tensor._raw(other.component,replacement,type(other))
                        if antisymmetric and slot==0 and index.is_up:coefficient=-coefficient
                        factors.pop(position);changed=True;break
                    if changed:break
                if changed:break
            if changed:break
        if not changed:break
    return _product(coefficient,factors)


def substitute_indices(t,*pairs):
    if not isinstance(t,TensExpr):return t
    if isinstance(t,TensAdd):return TensAdd(*(substitute_indices(term,*pairs) for term in t.args))
    if not isinstance(t,(Tensor,TensMul)):
        raise NotImplementedError('Index substitution on selected components is not implemented')
    mapping={}
    for old,new in pairs:
        if not isinstance(old,TensorIndex) or not isinstance(new,TensorIndex):raise TypeError('Index substitutions require typed indices')
        if old.tensor_index_type!=new.tensor_index_type:raise ValueError('An index substitution must preserve its type')
        if old in mapping and mapping[old]!=new:raise ValueError('Conflicting index substitutions')
        mapping[old]=new;mapping[-old]=-new
    indices=t.get_indices();structure=_IndexStructure(indices)
    target_names={_name_key(i) for i in mapping.values()};old_names={_name_key(i) for i in mapping}
    occupied={_name_key(i) for i in indices}|target_names;rename={}
    for upper,lower in structure.dum:
        index=indices[upper]
        if _name_key(index) in target_names and _name_key(index) not in old_names:
            counter=0
            while (index.tensor_index_type,'_tensor_subst_'+str(counter)) in occupied:counter+=1
            new=TensorIndex('_tensor_subst_'+str(counter),index.tensor_index_type);occupied.add(_name_key(new))
            rename[index]=new;rename[-index]=-new
    indices=[rename.get(index,index) for index in indices]
    return t._set_indices(*(mapping.get(index,index) for index in indices))


def riemann_cyclic_replace(t_r):
    from .. import Rational
    if not isinstance(t_r,Tensor) or t_r.component.rank!=4:raise ValueError('Riemann replacement requires a rank-four tensor')
    a,b,c,d=t_r.indices;head=t_r.component
    return TensAdd(TensMul(Rational(2,3),head(a,b,c,d)),
                   TensMul(Rational(-1,3),head(a,d,b,c)),
                   TensMul(Rational(1,3),head(a,c,b,d)))


def riemann_cyclic(t2):
    if not isinstance(t2,TensExpr):return t2
    if isinstance(t2,TensAdd):return canon_bp(TensAdd(*(riemann_cyclic(term) for term in t2.args)))
    coefficient,factors=_split(t2)
    replaced=[riemann_cyclic_replace(factor) if factor.component.symmetry==TensorSymmetry.riemann() else factor for factor in factors]
    return canon_bp(TensMul(coefficient,*replaced))
