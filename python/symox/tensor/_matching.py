"""Typed tensor wildcard matching with bounded factor backtracking."""
import itertools

from ._indexing import _Object,TensorIndex,TensorIndexType,TensorSymmetry
from ._expressions import TensorHead,Tensor,TensExpr,TensMul,TensAdd,_split,_product


class WildTensorIndex(TensorIndex):
    __slots__=()
    def __new__(cls,name,tensor_index_type,is_up=True,ignore_updown=False):
        result=super().__new__(cls,name,tensor_index_type,is_up)
        result._args+= (bool(ignore_updown),);return result
    @property
    def ignore_updown(self):return self.args[3]
    def __neg__(self):return type(self)(self.name,self.tensor_index_type,not self.is_up,self.ignore_updown)
    def matches(self,expr,repl_dict=None,old=False):
        if not isinstance(expr,TensorIndex) or expr.tensor_index_type!=self.tensor_index_type:return None
        if not self.ignore_updown and self.is_up!=expr.is_up:return None
        result=dict(repl_dict or {})
        if self in result and result[self]!=expr or -self in result and result[-self]!=-expr:return None
        result[self]=expr;return result
    def to_sympy(self):
        from sympy.tensor.tensor import WildTensorIndex
        return WildTensorIndex(self.name,self.tensor_index_type.to_sympy(),self.is_up,self.ignore_updown)


class WildTensorHead(TensorHead):
    __slots__=()
    def __new__(cls,name,index_types=None,symmetry=None,comm=0,unordered_indices=False):
        result=super().__new__(cls,name,[] if index_types is None else index_types,symmetry,comm)
        if any(not g.is_Identity for g in result.symmetry.generators):raise NotImplementedError('Wildcard tensor slot symmetries are not supported')
        result._args+=(bool(unordered_indices),);return result
    @property
    def unordered_indices(self):return self.args[4]
    def __call__(self,*indices,**kwargs):return WildTensor(self,indices,**kwargs)
    def to_sympy(self):
        from sympy.tensor.tensor import WildTensorHead
        return WildTensorHead(self.name,[t.to_sympy() for t in self.index_types],comm=self.comm,unordered_indices=self.unordered_indices)


class _WildTensExpr(_Object):
    __slots__=()
    def __new__(cls,expression,index_order=None):
        order=expression.get_free_indices() if index_order is None and isinstance(expression,TensExpr) else (() if index_order is None else index_order)
        result=object.__new__(cls);result._args=(expression,tuple(order));return result
    @property
    def expr(self):return self.args[0]
    def __call__(self,*indices):
        if not isinstance(self.expr,TensExpr):return self.expr
        if len(indices)!=len(self.args[1]):raise ValueError('A captured tensor requires its original number of free indices')
        return self.expr.substitute_indices(*zip(self.args[1],indices))
    def __repr__(self):return '_WildTensExpr(%s)' % self.expr
    def to_sympy(self):
        from sympy.tensor.tensor import _WildTensExpr
        return _WildTensExpr(self.expr.to_sympy())


class WildTensor(Tensor):
    __slots__=()
    def __new__(cls,tensor_head,indices,**kwargs):
        if type(tensor_head) is TensorHead:return Tensor(tensor_head,indices,**kwargs)
        if isinstance(tensor_head,_WildTensExpr):return tensor_head(*indices)
        indices=tuple(indices)
        if any(not isinstance(index,TensorIndex) for index in indices):raise TypeError('Wildcard tensors require typed indices')
        if tensor_head.rank and tensor_head.index_types!=[i.tensor_index_type for i in indices]:raise ValueError('Wildcard index types do not match')
        head=WildTensorHead(tensor_head.name,[i.tensor_index_type for i in indices],comm=tensor_head.comm,unordered_indices=tensor_head.unordered_indices)
        return super().__new__(cls,head,indices,**kwargs)
    def matches(self,expr,repl_dict=None,old=False):
        if not isinstance(expr,TensExpr) and expr!=1:return None
        result=dict(repl_dict or {})
        if not self.indices:
            if self in result and result[self]!=expr:return None
            result[self]=expr;return result
        if not isinstance(expr,TensExpr):return None
        free=expr.get_free_indices()
        if len(free)!=len(self.indices):return None
        if len(free)>8 and self.component.unordered_indices:raise NotImplementedError('Unordered wildcard matching is limited to eight free indices')
        orders=itertools.permutations(free) if self.component.unordered_indices else [free]
        for order in orders:
            candidate=result.copy()
            for pattern,index in zip(self.indices,order):
                candidate=pattern.matches(index,candidate,old)
                if candidate is None:break
            if candidate is not None:
                value=_WildTensExpr(expr,order)
                if self.component in candidate and candidate[self.component]!=value:continue
                candidate[self.component]=value;return candidate
        return None
    def __repr__(self):return super().__repr__() if self.indices else self.component.name
    def to_sympy(self):return self.component.to_sympy()(*(index.to_sympy() for index in self.indices))


def _match_tensor(pattern,expr,result,old):
    if not isinstance(expr,Tensor) or pattern.component!=expr.component:return None
    if len(pattern.indices)!=len(expr.indices):return None
    for left,right in zip(pattern.indices,expr.indices):
        result=left.matches(right,result,old)
        if result is None:return None
    return result


def matches(pattern,expr,repl_dict=None,old=False):
    result=dict(repl_dict or {})
    if isinstance(pattern,WildTensor):return pattern.matches(expr,result,old)
    if pattern==expr:return result
    if isinstance(pattern,Tensor):return _match_tensor(pattern,expr,result,old)
    if isinstance(pattern,TensAdd):
        if not isinstance(expr,TensAdd) or len(pattern.args)!=len(expr.args):return None
        if len(expr.args)>8:raise NotImplementedError('Tensor sum matching is limited to eight terms')
        for order in itertools.permutations(expr.args):
            candidate=result.copy()
            for left,right in zip(pattern.args,order):
                candidate=matches(left,right,candidate,old) if isinstance(left,TensExpr) else candidate if left==right else None
                if candidate is None:break
            if candidate is not None:return candidate
        return None
    if not isinstance(pattern,TensMul) or not isinstance(expr,TensExpr):return None
    pcoeff,pfactors=_split(pattern);ecoeff,efactors=_split(expr)
    flexible=[factor for factor in pfactors if isinstance(factor,WildTensor) and not factor.indices]
    if len(flexible)>1:raise NotImplementedError('Only one unrestricted tensor wildcard is supported per product')
    scalar_match=ecoeff.match(pcoeff)
    from .. import S,Wild
    capture_coefficient=S.One
    if scalar_match is None:
        if not flexible or pcoeff.has(Wild):return None
        capture_coefficient=ecoeff/pcoeff
    else:result.update(scalar_match)
    if len(efactors)>8:raise NotImplementedError('Tensor product matching is limited to eight factors')
    fixed=[factor for factor in pfactors if factor not in flexible]
    if not flexible and len(fixed)!=len(efactors):return None
    if len(fixed)>len(efactors):return None
    commuting=all(a.component.commutes_with(b.component)==0 for a in efactors for b in efactors)
    if commuting:orders=itertools.permutations(range(len(efactors)),len(fixed))
    elif flexible:
        position=pfactors.index(flexible[0]);suffix=len(pfactors)-position-1
        orders=[tuple(range(position))+tuple(range(len(efactors)-suffix,len(efactors)))]
    else:orders=[tuple(range(len(fixed)))]
    for order in orders:
        candidate=result.copy()
        for left,position in zip(fixed,order):
            candidate=left.matches(efactors[position],candidate,old)
            if candidate is None:break
        if candidate is None:continue
        if flexible:
            remaining=[factor for i,factor in enumerate(efactors) if i not in order]
            candidate=flexible[0].matches(_product(capture_coefficient,remaining),candidate,old)
        if candidate is not None:return candidate
    return None
