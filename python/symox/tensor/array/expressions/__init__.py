"""Symbolic, shaped array expressions with exact indexed evaluation.

These nodes preserve tensor operations without allocating their entries.
Explicit evaluation uses the owned native-scalar array implementation.
"""
from itertools import product
from math import prod
import operator
import symox as S
from .. import Array,NDimArray
from ...._kind import matrix_tag
from .... import _core


_NO_SHAPE=object()
def _shape(value):
    shape=getattr(value,'shape',_NO_SHAPE)
    if shape is _NO_SHAPE:return ()
    # native expressions are already sympified (sympify returns them unchanged)
    return tuple(n if isinstance(n,_core.Expr) else S.sympify(n) for n in shape)
def _maybe_shape(value):
    # None for an object whose shape is undefined (IndexedBase without
    # shape), as SymPy's get_shape; the array nodes then have shape None.
    try:
        if getattr(value,'shape',_NO_SHAPE) is None:return None
    except Exception:return ()
    return _shape(value)
def _entry(value,indices):
    if not indices:return S.sympify(value)
    return S.sympify(value[tuple(indices)])
def _dimensions(shape):
    shape=tuple(S.sympify(n) for n in shape)
    if any(n.is_integer is False or n.is_nonnegative is False for n in shape):raise ValueError('array dimensions must be nonnegative integers')
    return shape


class _ArrayExpr:
    is_commutative=False
    is_Array=True
    @property
    def func(self):return type(self)
    @property
    def free_symbols(self):
        return set().union(*(getattr(arg,'free_symbols',set()) for arg in self.args))
    def rank(self):return len(self.shape)
    def __getitem__(self,index):
        indices=index if isinstance(index,tuple) else (index,)
        if len(indices)!=len(self.shape):raise IndexError('wrong number of array indices')
        normalized=[]
        for i,n in zip(indices,self.shape):
            i=S.sympify(i)
            if i.is_integer is False:raise IndexError('array indices must be integers')
            n=S.sympify(n)
            if i.is_number and n.is_number:
                if i<0:i+=n
                if not 0<=i<n:raise IndexError('array index out of range')
            normalized.append(i)
        return self._entry(tuple(normalized))
    def as_explicit(self):
        shape=tuple(operator.index(n) for n in self.shape)
        if prod(shape)>100000:raise ValueError('explicit array exceeds 100000 entries')
        return Array([self[index] for index in product(*(range(n) for n in shape))],shape)
    def __add__(self,other):return ArrayAdd(self,other)
    __radd__=__add__
    def __sub__(self,other):return ArrayAdd(self,ArrayTensorProduct(-1,other))
    def __mul__(self,other):
        if _shape(other):raise TypeError('use ArrayTensorProduct for array multiplication')
        return ArrayTensorProduct(other,self)
    __rmul__=__mul__
    def __neg__(self):return ArrayTensorProduct(-1,self)
    def __eq__(self,other):return type(self) is type(other) and self.args==other.args
    def __hash__(self):return hash((type(self),self.args))
    def dummy_eq(self,other,symbol=None):
        def eq(a,b):
            if isinstance(a,tuple):return isinstance(b,tuple) and len(a)==len(b) and all(eq(u,v) for u,v in zip(a,b))
            if hasattr(a,'dummy_eq'):
                try:return a.dummy_eq(b) if isinstance(a,_ArrayExpr) else a.dummy_eq(b,symbol)
                except (TypeError,AttributeError):return False
            return a==b
        return type(self) is type(other) and eq(self.args,other.args)
    def __repr__(self):return type(self).__name__+'('+', '.join(map(str,self.args))+')'
    def doit(self,**hints):return self
    def subs(self,*args,**kwargs):
        mapping=args[0] if len(args)==1 and isinstance(args[0],dict) else {args[0]:args[1]}
        if self in mapping:return mapping[self]
        def replace(value):
            if isinstance(value,tuple):return tuple(replace(v) for v in value)
            return value.subs(mapping,**kwargs) if hasattr(value,'subs') else value
        return type(self)(*(replace(arg) for arg in self.args))


class ArraySymbol(_ArrayExpr):
    def __init__(self,name,shape):self.name=str(name);self.shape=_dimensions(shape);self.args=(self.name,self.shape)
    @property
    def free_symbols(self):return {self}
    def __repr__(self):return self.name
    def _entry(self,indices):return ArrayElement(self,indices)


from ....codegen.ast import _IndexedArithmetic

class ArrayElement(_IndexedArithmetic):
    def __init__(self,name,indices):
        # SymPy's ArrayElement.__new__ and _check_shape
        if isinstance(name,str):name=S.Symbol(name)
        if not isinstance(indices,(list,tuple)):indices=(indices,)
        self.name,self.indices=name,tuple(S.sympify(i) for i in indices)
        self._check_shape(name,self.indices)
    @classmethod
    def _check_shape(cls,name,indices):
        indices=tuple(indices)
        if hasattr(name,'shape'):
            if len(indices)!=len(name.shape):raise IndexError('number of indices does not match shape of the array')
            if any((i>=s)==True for i,s in zip(indices,name.shape)):raise ValueError('shape is out of bounds')
        if any((i<0)==True for i in indices):raise ValueError('shape contains negative values')
    args=property(lambda self:(self.name,self.indices))
    @property
    def free_symbols(self):return {self.name}|set().union(*(i.free_symbols for i in self.indices))
    def _symox_(self):
        return _array_element_node(self.name)(S.Symbol(str(getattr(self.name,'name',self.name))),*self.indices)
    def diff(self,*symbols,**assumptions):
        # The native node differentiates (SymPy ArrayElement._eval_derivative).
        return S.diff(self._symox_(),*(s._symox_() if isinstance(s,ArrayElement) else s for s in symbols),**assumptions)
    def __repr__(self):return str(getattr(self.name,'name',self.name))+'['+', '.join(map(str,self.indices))+']'
    def __eq__(self,other):return isinstance(other,ArrayElement) and self.args==other.args
    def __hash__(self):return hash((type(self),self.args))


_ARRAY_ELEMENT_NODES={}
def _array_element_node(parent):
    # The scalar form of an ArrayElement: an applied function printed as
    # SymPy prints ArrayElement (``A[0]``, LaTeX ``{{A}_{0}}``). One class
    # per parent array, so that elements of distinct arrays sharing a name
    # stay distinct (SymPy compares ArrayElement.name, the array itself) and
    # the parent is recovered by convert_indexed_to_array.
    try:key=(type(parent),parent);cls=_ARRAY_ELEMENT_NODES.get(key)
    except TypeError:key=('id',id(parent));cls=_ARRAY_ELEMENT_NODES.get(key)
    if cls is None:
        class ArrayElement(S.Function):
            _symox_user_hooks=True
            # Differentiated natively as SymPy's ArrayElement (_diff_wrt).
            _symox_array_element=True
            _symox_element_parent=parent
            def _sympystr(self,printer):return printer._print(self.args[0])+'['+', '.join(printer._print(i) for i in self.args[1:])+']'
            def _latex(self,printer):return '{{%s}_{%s}}'%(printer._print(self.args[0]),', '.join(printer._print(i) for i in self.args[1:]))
        cls=_ARRAY_ELEMENT_NODES[key]=ArrayElement
    return cls


class _ConstantArray(_ArrayExpr):
    value=None
    def __init__(self,*shape):
        if len(shape)==1 and isinstance(shape[0],(list,tuple)):shape=shape[0]
        self.shape=_dimensions(shape);self.args=self.shape
    def _entry(self,indices):return S.Integer(self.value)
    @property
    def free_symbols(self):return set().union(*(n.free_symbols for n in self.shape))
class ZeroArray(_ConstantArray):value=0
class OneArray(_ConstantArray):value=1


class ArrayTensorProduct(_ArrayExpr):
    def __init__(self,*args):
        flattened=[]
        for arg in args:flattened.extend(arg.args if isinstance(arg,ArrayTensorProduct) else [arg])
        self.args=tuple(flattened);shapes=[_maybe_shape(a) for a in self.args]
        self.shape=None if any(sh is None for sh in shapes) else tuple(n for sh in shapes for n in sh)
    def doit(self,**hints):
        from ._canon import _doit
        return _doit(self,hints.get('deep',True))
    def _entry(self,indices):
        result=S.Integer(1);offset=0
        for arg in self.args:
            rank=len(_shape(arg));result*=_entry(arg,indices[offset:offset+rank]);offset+=rank
        return result


class ArrayAdd(_ArrayExpr):
    def __init__(self,*args):
        if not args:raise ValueError('an array sum requires arguments')
        self.shape=_shape(args[0])
        if any(_shape(a)!=self.shape for a in args):raise ValueError('array sum shapes must match')
        self.args=tuple(args)
    def doit(self,**hints):
        from ._canon import _doit
        return _doit(self,hints.get('deep',True))
    def _entry(self,indices):return S.Add(*(_entry(a,indices) for a in self.args))


class PermuteDims(_ArrayExpr):
    def __init__(self,expr,permutation=None,index_order_old=None,index_order_new=None):
        if index_order_old is not None:permutation=[list(index_order_old).index(i) for i in index_order_new]
        permutation=tuple(permutation.array_form if hasattr(permutation,'array_form') else permutation)
        if sorted(permutation)!=list(range(len(_shape(expr)))):raise ValueError('permutation must contain each axis exactly once')
        from ....combinatorics import Permutation
        self.expr,self.permutation=expr,Permutation(list(permutation))
        self.args=(expr,self.permutation);self.shape=tuple(expr.shape[i] for i in permutation)
    def doit(self,**hints):
        from ._canon import _doit
        return _doit(self,hints.get('deep',True))
    def nest_permutation(self):
        from ._canon import nest_permutation
        return nest_permutation(self)
    def _entry(self,indices):
        original=[None]*len(indices)
        for axis,old in enumerate(self.permutation.array_form):original[old]=indices[axis]
        return _entry(self.expr,original)


def _af_invert(form):
    inverse=[0]*len(form)
    for i,j in enumerate(form):inverse[j]=i
    return inverse


def _permute_dims(expr,permutation):
    """``PermuteDims(expr, permutation)`` in SymPy's canonical form.

    As SymPy's ``PermuteDims._canonicalize``: nested permutations are
    composed, the arguments of a permuted tensor product are sorted by the
    positions their axes take, and an identity permutation is dropped.
    """
    from ....combinatorics import Permutation
    form=list(permutation.array_form if hasattr(permutation,'array_form') else permutation)
    if sorted(form)!=list(range(len(_shape(expr)))):raise ValueError('Permutation size must be the length of the shape of expr')
    if isinstance(expr,PermuteDims):
        inner=expr.permutation.array_form
        form=[inner[i] for i in form];expr=expr.expr
    if isinstance(expr,ArrayContraction):
        expr,form=_permute_contraction(expr,form)
    if isinstance(expr,ArrayTensorProduct):
        image=_af_invert(form);args=list(expr.args);offsets=[0]
        for arg in args:offsets.append(offsets[-1]+len(_shape(arg)))
        components=[image[offsets[i]:offsets[i+1]] for i in range(len(args))]
        order=sorted(range(len(args)),key=lambda i:sorted(components[i]))
        form=_af_invert([j for i in order for j in components[i]])
        expr=ArrayTensorProduct(*(args[i] for i in order)) if len(args)>1 else args[0]
    if isinstance(expr,ZeroArray):return ZeroArray(*[expr.shape[i] for i in form])
    if form==sorted(form):return expr
    return PermuteDims(expr,Permutation(form))


def _sort_key(value):
    if isinstance(value,ArraySymbol):
        return ((5,0,'ArraySymbol'),(1,(value.name,)+tuple(str(n) for n in value.shape)))
    try:return S.default_sort_key(value)
    except Exception:return (str(type(value).__name__),str(value))


def _array_contraction(expr,*groups):
    """``ArrayContraction(expr, *groups)`` in part of SymPy's canonical form.

    As SymPy's ``ArrayContraction._canonicalize`` for a tensor product: the
    fully contracted arguments move first, sorted by their sort key, and a
    single-axis contraction of a dimension one axis is dropped. Other nested
    forms are kept as given (an equivalent, not canonical, expression).
    """
    groups=[tuple(sorted(operator.index(i) for i in g)) for g in groups]
    if not groups:return expr
    if isinstance(expr,ArrayTensorProduct) and not any(isinstance(a,ArrayAdd) for a in expr.args):
        args=expr.args;offsets=[0]
        for a in args:offsets.append(offsets[-1]+len(_shape(a)))
        flat={j for g in groups for j in g}
        full=[all(j in flat for j in range(offsets[i],offsets[i+1])) for i in range(len(args))]
        order=sorted(range(len(args)),key=lambda i:(0,_sort_key(args[i])) if full[i] else (1,))
        form=_af_invert([j for i in order for j in range(offsets[i],offsets[i+1])])
        groups=[tuple(sorted(form[j] for j in g)) for g in groups]
        expr=ArrayTensorProduct(*(args[i] for i in order))
    shape=_shape(expr)
    groups=[g for g in groups if len(g)>1 or shape[g[0]]!=1]
    if not groups:return expr
    return ArrayContraction(expr,*groups)


def _is_symbolic_array(value):
    """Whether SymPy's array operations build an array expression for it."""
    return isinstance(value,_ArrayExpr) or matrix_tag(value)=='MatrixSymbol'


class _AxisOperation(_ArrayExpr):
    def __init__(self,expr,*groups):
        self.expr=expr;shape=_maybe_shape(expr);used=set();normalized=[]
        if shape is None:
            self.groups=tuple(tuple(operator.index(i) for i in g) for g in groups);self.remaining=None
            self.args=(expr,*self.groups);self.shape=None;return
        for group in groups:
            group=tuple(operator.index(i) for i in group)
            if not group or len(set(group))!=len(group) or used.intersection(group) or any(i<0 or i>=len(shape) for i in group):raise ValueError('invalid or overlapping array axes')
            if len({shape[i] for i in group})!=1:raise ValueError('contracted or diagonal axes must have equal dimensions')
            used.update(group);normalized.append(group)
        self.groups=tuple(normalized);self.remaining=tuple(i for i in range(len(shape)) if i not in used)
        self.args=(expr,*self.groups);self.shape=tuple(shape[i] for i in self.remaining)

class ArrayContraction(_AxisOperation):
    def __init__(self,expr,*groups):
        # SymPy sorts each contraction and the contractions by first axis
        groups=sorted((tuple(sorted(operator.index(i) for i in g)) for g in groups),key=lambda g:min(g) if g else -1)
        super().__init__(expr,*groups)
    def doit(self,**hints):
        from ._canon import _doit
        return _doit(self,hints.get('deep',True))
    @property
    def contraction_indices(self):return self.groups
    def sort_args_by_name(self):
        # SymPy's ArrayContraction.sort_args_by_name
        expr=self.expr
        if not isinstance(expr,ArrayTensorProduct):return self
        args=expr.args;offsets=[0]
        for a in args:offsets.append(offsets[-1]+len(_shape(a)))
        owner=[(k,j-offsets[k]) for k in range(len(args)) for j in range(offsets[k],offsets[k+1])]
        order=sorted(range(len(args)),key=lambda i:_sort_key(args[i]))
        new_position={old:new for new,old in enumerate(order)}
        sorted_args=[args[i] for i in order];new_offsets=[0]
        for a in sorted_args:new_offsets.append(new_offsets[-1]+len(_shape(a)))
        groups=[tuple(new_offsets[new_position[owner[j][0]]]+owner[j][1] for j in g) for g in self.groups]
        return _array_contraction(ArrayTensorProduct(*sorted_args),*groups)
    def _entry(self,indices):
        original=[None]*len(_shape(self.expr))
        for axis,index in zip(self.remaining,indices):original[axis]=index
        dummies=[S.Dummy('i',integer=True) for _ in self.groups]
        limits=[]
        for group,index in zip(self.groups,dummies):
            for axis in group:original[axis]=index
            limits.append((index,0,self.expr.shape[group[0]]-1))
        if all(self.expr.shape[group[0]].is_number if hasattr(self.expr.shape[group[0]],'is_number') else True for group in self.groups):
            terms=[]
            for values in product(*(range(int(self.expr.shape[g[0]])) for g in self.groups)):
                at=list(original)
                for group,index in zip(self.groups,values):
                    for axis in group:at[axis]=index
                terms.append(_entry(self.expr,at))
            return S.Add(*terms)
        expression=_entry(self.expr,original)
        return S.Sum(expression,*limits) if limits else expression

class ArrayDiagonal(_AxisOperation):
    def __init__(self,expr,*groups):
        groups=[tuple(sorted(operator.index(i) for i in g)) for g in groups]
        super().__init__(expr,*groups)
        if any(len(g)<2 for g in self.groups):raise ValueError('diagonals require at least two axes')
        if self.shape is not None:self.shape+=tuple(expr.shape[g[0]] for g in self.groups)
    def doit(self,**hints):
        from ._canon import _doit
        return _doit(self,hints.get('deep',True))
    @property
    def diagonal_indices(self):return self.groups
    def _entry(self,indices):
        original=[None]*len(_shape(self.expr))
        for axis,index in zip(self.remaining,indices):original[axis]=index
        for group,index in zip(self.groups,indices[len(self.remaining):]):
            for axis in group:original[axis]=index
        return _entry(self.expr,original)


class Reshape(_ArrayExpr):
    def __init__(self,expr,shape):
        self.expr=expr;self.shape=_dimensions(shape);self.args=(expr,self.shape)
        if S.simplify(prod(_shape(expr))-prod(self.shape))!=0:raise ValueError('reshape must preserve the number of entries')
    def _entry(self,indices):
        linear=S.Integer(0)
        for index,n in zip(indices,self.shape):linear=linear*n+index
        original=[]
        for n in reversed(_shape(self.expr)):
            original.append(S.Mod(linear,n));linear=S.floor(linear/n)
        return _entry(self.expr,tuple(reversed(original)))


class ArrayElementwiseApplyFunc(_ArrayExpr):
    def __init__(self,function,expr):
        if not isinstance(function,S.Lambda):
            d=S.Dummy('d');function=S.Lambda(d,function(d))
        self.function,self.expr=function,expr;self.args=(function,expr);self.shape=_maybe_shape(expr)
    def _entry(self,indices):return self.function(_entry(self.expr,indices))
    def _get_function_fdiff(self):
        # SymPy's ArrayElementwiseApplyFunc._get_function_fdiff
        d=S.Dummy('d');fdiff=self.function(d).diff(d)
        return type(fdiff) if isinstance(fdiff,S.Function) else S.Lambda(d,fdiff)


def array_derive(expr,x):
    from ._native import array_derive
    return array_derive(expr,x)


def convert_array_to_indexed(expr,indices):
    from ._native import convert_array_to_indexed
    return convert_array_to_indexed(expr,indices)


def convert_matrix_to_array(expr):
    from ._canon import convert_matrix_to_array
    return convert_matrix_to_array(expr)


def _as_matrix_factor(value):
    if isinstance(value,_ArrayExpr):value=convert_array_to_matrix(value)
    if isinstance(value,_ArrayExpr) or len(_shape(value))!=2 or not hasattr(value,'T'):return None
    if isinstance(value,NDimArray):return None
    return value


def _matrix_chains(expr):
    """Matrix form of a contraction of a tensor product of matrices.

    Each connected set of contracted matrix axes is a chain (a matrix
    product whose two ends are free axes) or a cycle (a trace). The result
    is None when the free axes of the chains are not consecutive, which
    needs an axis permutation that a matrix product cannot express.
    """
    from functools import reduce
    inner,groups=(expr.expr,expr.groups) if isinstance(expr,ArrayContraction) else (expr,())
    factors=inner.args if isinstance(inner,ArrayTensorProduct) else (inner,)
    if any(len(g)!=2 for g in groups):return None
    owner={};matrices={};scalars=[];axis=0
    for k,factor in enumerate(factors):
        rank=len(_shape(factor))
        if rank==0:
            scalars.append(S.sympify(factor));continue
        matrix=_as_matrix_factor(factor)
        if matrix is None:return None
        matrices[k]=matrix;owner[axis]=(k,0);owner[axis+1]=(k,1);axis+=2
    position={v:a for a,v in owner.items()}
    partner={}
    for a,b in groups:partner[a]=b;partner[b]=a
    free=[a for a in range(axis) if a not in partner]
    def walk(k,side):
        chain=[];start=(k,side)
        while True:
            chain.append(matrices[k] if side==0 else matrices[k].T);seen.add(k)
            out=position[(k,1-side)]
            if out not in partner:return chain,out
            k,side=owner[partner[out]]
            if (k,side)==start:return chain,None
    seen=set();chains=[];ends=[];pending=list(free)
    while pending:
        start=pending.pop(0);k,side=owner[start]
        chain,end=walk(k,side)
        if end is None or end not in pending:return None
        pending.remove(end);ends.extend((start,end))
        chains.append(reduce(lambda a,b:a*b,chain))
    for k in matrices:
        if k not in seen:
            chain,_=walk(k,0)
            # The trace of a cycle is evaluated (Trace(A.T) is Trace(A), a
            # scalar factor leaves the trace), as SymPy's recognized traces.
            scalars.append(S.Trace(reduce(lambda a,b:a*b,chain)).doit(deep=False))
    scalar=S.Mul(*scalars)
    if not chains:return scalar
    if scalar!=1:chains[0]=scalar*chains[0]
    if len(chains)==1:return chains[0]
    result=ArrayTensorProduct(*chains)
    # axis t of the product is free axis ends[t] of the expression
    form=_af_invert([free.index(e) for e in ends])
    return _permute_dims(result,form)


def convert_array_to_matrix(expr):
    if isinstance(expr,(ArrayContraction,ArrayTensorProduct)):
        chained=_matrix_chains(expr)
        if chained is not None:return chained
    if len(_shape(expr))!=2:return expr
    if isinstance(expr,ArraySymbol):return S.MatrixSymbol(expr.name,*expr.shape)
    if isinstance(expr,ArrayAdd):return S.MatAdd(*(convert_array_to_matrix(a) for a in expr.args))
    if isinstance(expr,PermuteDims) and expr.permutation.array_form==[1,0]:return convert_array_to_matrix(expr.expr).T
    if isinstance(expr,ArrayContraction) and isinstance(expr.expr,ArrayTensorProduct):
        factors=expr.expr.args;matrices=[a for a in factors if len(_shape(a))==2]
        if expr.groups==tuple((2*i+1,2*i+2) for i in range(len(matrices)-1)):
            return S.MatMul(*(convert_array_to_matrix(a) if _shape(a) else a for a in factors))
    if isinstance(expr,_ArrayExpr):
        i,j=S.symbols('_i _j',integer=True)
        return S.FunctionMatrix(*expr.shape,S.Lambda((i,j),expr[i,j]))
    return expr


def _permute_contraction(expr,form):
    """SymPy's ``PermuteDims._PermuteDims_denestarg_ArrayContraction``."""
    if not isinstance(expr.expr,ArrayTensorProduct):return expr,form
    args=expr.expr.args;ranks=[len(_shape(a)) for a in args];offsets=[0]
    for r in ranks:offsets.append(offsets[-1]+r)
    contracted={j for g in expr.groups for j in g};image=_af_invert(form)
    blocks_up=[];counter=0;index_blocks_up=[]
    for i in range(len(args)):
        current=[];up=[]
        for j in range(offsets[i],offsets[i+1]):
            if j in contracted:continue
            current.append(image[counter]);up.append(image[counter]);counter+=1
        blocks_up.append(current);index_blocks_up.append(up)
    order=sorted(range(len(args)),key=lambda i:index_blocks_up[i])
    new_blocks=[list(range(offsets[i],offsets[i+1])) for i in order]
    new_form=_af_invert([j for b in new_blocks for j in b])
    groups=[tuple(new_form[j] for j in g) for g in expr.groups]
    new_args=[args[i] for i in order]
    product_=ArrayTensorProduct(*new_args) if len(new_args)>1 else new_args[0]
    return _array_contraction(product_,*groups),_af_invert([j for i in order for j in blocks_up[i]])


def convert_indexed_to_array(expr,first_indices=None):
    from ._native import convert_indexed_to_array
    return convert_indexed_to_array(expr,first_indices)


__all__=['ArrayAdd','ArrayContraction','ArrayDiagonal','ArrayElement','ArrayElementwiseApplyFunc','ArraySymbol',
 'ArrayTensorProduct','OneArray','PermuteDims','Reshape','ZeroArray','array_derive','convert_array_to_indexed',
 'convert_array_to_matrix','convert_indexed_to_array','convert_matrix_to_array']
