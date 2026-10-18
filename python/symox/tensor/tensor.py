"""Typed indexed tensor algebra over native scalar and array algorithms."""
import warnings

from ._indexing import TensorIndexType,TensorIndex,TensorSymmetry,_Object,_IndexStructure
from ._expressions import TensorHead,Tensor,TensExpr,TensAdd,TensMul,TensorElement,tensor_heads
from ._canonical import canon_bp,contract_metric,substitute_indices,riemann_cyclic,riemann_cyclic_replace
from ._matching import WildTensorIndex,WildTensorHead,WildTensor


def tensor_indices(s,typ):
    from ._indexing import tensor_indices as construct
    return construct(s,typ)


def tensor_mul(*a):return TensMul(*a)
def get_indices(t):return t.get_indices() if isinstance(t,TensExpr) else ()
def get_free_indices(t):return t.get_free_indices() if isinstance(t,TensExpr) else ()
def get_dummy_indices(t):
    if not isinstance(t,TensExpr):return ()
    free=t.get_free_indices()
    return [i for i in t.get_indices() if i not in free]
def get_index_structure(t):
    if isinstance(t,TensAdd):raise NotImplementedError('A tensor sum has no single monomial index structure')
    return t._index_structure if isinstance(t,TensExpr) else _IndexStructure([])
def get_coeff(t):
    if isinstance(t,TensAdd):raise ValueError('A tensor sum has no unique scalar coefficient')
    return t.coeff if isinstance(t,TensExpr) else t


def components_canon_args(components):
    grouped=[]
    for head in components:
        if grouped and grouped[-1][0]==head:grouped[-1][1]+=1
        else:grouped.append([head,1])
    return [(head.symmetry.base,head.symmetry.generators,count,head.comm if head.comm in (0,1) else None)
            for head,count in grouped]


def perm2tensor(t,g,is_canon_bp=False):
    if not isinstance(t,TensExpr):return t
    if not isinstance(t,(Tensor,TensMul)):raise NotImplementedError('Index permutations require a tensor monomial')
    structure=get_index_structure(t).perm2tensor(g,is_canon_bp)
    result=t._set_indices(*structure.indices,is_canon_bp=is_canon_bp)
    images=g.array_form if hasattr(g,'array_form') else list(g)
    return -result if images[-1]!=len(images)-1 else result


def get_lines(ex,index_type):
    if not isinstance(ex,(Tensor,TensMul)):raise ValueError('Matrix lines require a tensor monomial')
    arguments=ex.args if isinstance(ex,TensMul) else (ex,)
    slots={}
    for position,tensor in enumerate(arguments):
        if not isinstance(tensor,Tensor):continue
        positions=[i for i,t in enumerate(tensor.component.index_types) if t==index_type]
        if len(positions)>2:raise ValueError('Matrix lines require at most two slots of the selected type')
        if len(positions)==2:slots[position]=positions
    edges={};incoming={};traces=[]
    for up,down,upper_arg,lower_arg in ex.dum_in_args:
        if upper_arg not in slots or lower_arg not in slots:continue
        if up not in slots[upper_arg] or down not in slots[lower_arg]:continue
        if upper_arg==lower_arg:traces.append([upper_arg]);continue
        first,second=slots[upper_arg].index(up),slots[lower_arg].index(down)
        if first==second:raise NotImplementedError('Matrix-line transposition flags are not supported')
        left,right=(upper_arg,lower_arg) if first==1 else (lower_arg,upper_arg)
        if left in edges and edges[left]!=right or right in incoming and incoming[right]!=left:
            raise ValueError('The contraction graph is not a matrix line')
        edges[left]=right;incoming[right]=left
    lines=[];used={vertex for trace in traces for vertex in trace}
    vertices=set(edges)|set(incoming)
    for start in sorted(vertices):
        if start in incoming or start in used:continue
        line=[start];used.add(start)
        while line[-1] in edges:
            following=edges[line[-1]]
            if following in used:raise ValueError('Invalid matrix-line contraction graph')
            line.append(following);used.add(following)
        lines.append(line)
    for start in sorted(vertices-used):
        if start in used:continue
        trace=[start];used.add(start)
        while edges.get(trace[-1])!=start:
            following=edges.get(trace[-1])
            if following is None or following in used:raise ValueError('Invalid traced matrix line')
            trace.append(following);used.add(following)
        traces.append(trace)
    return lines,traces,[i for i in range(len(arguments)) if i not in used]


def get_postprocessor(cls):
    name=getattr(cls,'__name__','')
    if name not in ('Add','Mul'):raise ValueError('Tensor postprocessors apply to Add or Mul')
    constructor=TensAdd if name=='Add' else TensMul
    def postprocess(expr):
        return constructor(*expr.args) if any(isinstance(arg,TensExpr) for arg in expr.args) else expr
    return postprocess


def deprecate_call():warnings.warn('Calling tensor expressions is deprecated; use substitute_indices()',DeprecationWarning,stacklevel=2)
def deprecate_data():warnings.warn('Tensor .data is deprecated; use replace_with_arrays()',DeprecationWarning,stacklevel=2)
def deprecate_fun_eval():warnings.warn('Tensor fun_eval() is deprecated; use substitute_indices()',DeprecationWarning,stacklevel=2)


def tensorsymmetry(*args):
    warnings.warn('tensorsymmetry() is deprecated; use TensorSymmetry',DeprecationWarning,stacklevel=2)
    from ..combinatorics import Permutation
    if not args:return TensorSymmetry.no_symmetry(0)
    if len(args)==2 and args[1] and isinstance(args[1][0],Permutation):return TensorSymmetry(*args)
    if len(args)==1 and len(args[0])==2 and isinstance(args[0][1],(tuple,list)):
        return TensorSymmetry(args[0])
    groups=[]
    for shape in args:
        shape=list(shape)
        if shape==[2,2]:group=TensorSymmetry.riemann()
        elif len(shape)==1:group=TensorSymmetry.fully_symmetric(-shape[0])
        elif all(value==1 for value in shape):group=TensorSymmetry.fully_symmetric(len(shape))
        else:raise NotImplementedError('This Young diagram symmetry is not supported')
        groups.append(group)
    size=sum(group.rank for group in groups);base=[];generators=[];offset=0
    for group in groups:
        base.extend(offset+i for i in group.base)
        for generator in group.generators:
            local=generator.array_form;images=list(range(size+2))
            for i in range(group.rank):images[offset+i]=offset+local[i]
            if local[-1]!=len(local)-1:images[-2:]=reversed(images[-2:])
            generators.append(Permutation(images))
        offset+=group.rank
    return TensorSymmetry(base,generators)


def tensorhead(name,typ,sym=None,comm=0):
    warnings.warn('tensorhead() is deprecated; use tensor_heads()',DeprecationWarning,stacklevel=2)
    if sym is None:symmetry=None
    elif isinstance(sym,TensorSymmetry):symmetry=sym
    else:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',DeprecationWarning);symmetry=tensorsymmetry(*sym)
    return tensor_heads(name,typ,symmetry,comm)


class TensorType(_Object):
    __slots__=()
    is_commutative=False
    def __new__(cls,index_types,symmetry,**kwargs):
        warnings.warn('TensorType is deprecated; use TensorHead or tensor_heads()',DeprecationWarning,stacklevel=2)
        checked=TensorHead('_type',index_types,symmetry)
        result=object.__new__(cls);result._args=(tuple(checked.index_types),symmetry);return result
    @property
    def index_types(self):return list(self.args[0])
    @property
    def symmetry(self):return self.args[1]
    @property
    def types(self):return sorted(set(self.index_types))
    def __call__(self,s,comm=0):return tensor_heads(s,self.index_types,self.symmetry,comm)
    def to_sympy(self):
        from sympy.tensor.tensor import TensorType
        return TensorType([t.to_sympy() for t in self.index_types],self.symmetry.to_sympy())


__all__=['TensAdd','TensExpr','TensMul','Tensor','TensorElement','TensorHead','TensorIndex',
         'TensorIndexType','TensorSymmetry','TensorType','WildTensor','WildTensorHead','WildTensorIndex',
         'canon_bp','components_canon_args','contract_metric','deprecate_call','deprecate_data',
         'deprecate_fun_eval','get_coeff','get_dummy_indices','get_free_indices','get_index_structure',
         'get_indices','get_lines','get_postprocessor','perm2tensor','riemann_cyclic',
         'riemann_cyclic_replace','substitute_indices','tensor_heads','tensor_indices','tensor_mul',
         'tensorhead','tensorsymmetry']
