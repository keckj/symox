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


def _shape(value):return tuple(S.sympify(n) for n in value.shape) if hasattr(value,'shape') else ()
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
    def _entry(self,indices):return ArrayElement(self,indices)


from ....codegen.ast import _IndexedArithmetic

class ArrayElement(_IndexedArithmetic):
    def __init__(self,name,indices):
        if not hasattr(name,'shape'):raise TypeError('an array element requires a shaped parent')
        self.name,self.indices=name,tuple(S.sympify(i) for i in indices)
        if len(self.indices)!=len(name.shape):raise IndexError('wrong number of array indices')
    args=property(lambda self:(self.name,self.indices))
    @property
    def free_symbols(self):return {self.name}|set().union(*(i.free_symbols for i in self.indices))
    def _symox_(self):return S.Function('ArrayElement')(S.Symbol(self.name.name),*self.indices)
    def __repr__(self):return str(self.name.name)+'['+', '.join(map(str,self.indices))+']'
    def __eq__(self,other):return isinstance(other,ArrayElement) and self.args==other.args
    def __hash__(self):return hash((type(self),self.args))


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
        self.args=tuple(flattened);self.shape=tuple(n for a in self.args for n in _shape(a))
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
    def _entry(self,indices):return S.Add(*(_entry(a,indices) for a in self.args))


class PermuteDims(_ArrayExpr):
    def __init__(self,expr,permutation=None,index_order_old=None,index_order_new=None):
        if index_order_old is not None:permutation=[list(index_order_old).index(i) for i in index_order_new]
        permutation=tuple(permutation.array_form if hasattr(permutation,'array_form') else permutation)
        if sorted(permutation)!=list(range(len(_shape(expr)))):raise ValueError('permutation must contain each axis exactly once')
        self.expr,self.permutation=expr,permutation;self.args=(expr,permutation);self.shape=tuple(expr.shape[i] for i in permutation)
    def _entry(self,indices):
        original=[None]*len(indices)
        for axis,old in enumerate(self.permutation):original[old]=indices[axis]
        return _entry(self.expr,original)


class _AxisOperation(_ArrayExpr):
    def __init__(self,expr,*groups):
        self.expr=expr;shape=_shape(expr);used=set();normalized=[]
        for group in groups:
            group=tuple(operator.index(i) for i in group)
            if not group or len(set(group))!=len(group) or used.intersection(group) or any(i<0 or i>=len(shape) for i in group):raise ValueError('invalid or overlapping array axes')
            if len({shape[i] for i in group})!=1:raise ValueError('contracted or diagonal axes must have equal dimensions')
            used.update(group);normalized.append(group)
        self.groups=tuple(normalized);self.remaining=tuple(i for i in range(len(shape)) if i not in used)
        self.args=(expr,*self.groups);self.shape=tuple(shape[i] for i in self.remaining)

class ArrayContraction(_AxisOperation):
    @property
    def contraction_indices(self):return self.groups
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
        super().__init__(expr,*groups)
        if any(len(g)<2 for g in self.groups):raise ValueError('diagonals require at least two axes')
        self.shape+=tuple(expr.shape[g[0]] for g in self.groups)
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
    def __init__(self,function,expr):self.function,self.expr=function,expr;self.args=(function,expr);self.shape=_shape(expr)
    def _entry(self,indices):return self.function(_entry(self.expr,indices))


class _ArrayDerivative(_ArrayExpr):
    def __init__(self,expr,variable):self.expr,self.variable=expr,variable;self.args=(expr,variable);self.shape=_shape(variable)+_shape(expr)
    def _entry(self,indices):
        rank=len(_shape(self.variable));deriving,output=indices[:rank],indices[rank:]
        if isinstance(self.expr,ArraySymbol) and isinstance(self.variable,ArraySymbol):
            return S.Mul(*(S.KroneckerDelta(i,j) for i,j in zip(deriving,output))) if self.expr==self.variable else S.Integer(0)
        return _entry(self.expr,output).diff(_entry(self.variable,deriving) if rank else self.variable)


def array_derive(expr,x):
    if not _shape(x):
        return ArrayElementwiseApplyFunc(lambda value:value.diff(x),expr) if _shape(expr) else S.diff(expr,x)
    return _ArrayDerivative(expr,x)


def convert_array_to_indexed(expr,indices):
    if len(indices)!=len(_shape(expr)):raise ValueError('indexed conversion requires one index for each free axis')
    return _entry(expr,indices)


def convert_matrix_to_array(expr):
    kind=matrix_tag(expr)
    if kind=='MatrixSymbol':return ArraySymbol(expr.name,expr.shape)
    if kind=='MatAdd':return ArrayAdd(*(convert_matrix_to_array(a) for a in expr.args))
    if kind=='Transpose':return PermuteDims(convert_matrix_to_array(expr.args[0]),(1,0))
    if kind=='MatMul':
        factors=[convert_matrix_to_array(a) if hasattr(a,'shape') else a for a in expr.args]
        product_=ArrayTensorProduct(*factors);matrices=[a for a in factors if _shape(a)]
        groups=[(2*i+1,2*i+2) for i in range(len(matrices)-1)]
        return ArrayContraction(product_,*groups) if groups else product_
    if kind=='KroneckerProduct':
        product_=ArrayTensorProduct(*(convert_matrix_to_array(a) for a in expr.args));n=len(expr.args)
        order=tuple(range(0,2*n,2))+tuple(range(1,2*n,2))
        return Reshape(PermuteDims(product_,order),expr.shape)
    return expr


def convert_array_to_matrix(expr):
    if len(_shape(expr))!=2:return expr
    if isinstance(expr,ArraySymbol):return S.MatrixSymbol(expr.name,*expr.shape)
    if isinstance(expr,ArrayAdd):return S.MatAdd(*(convert_array_to_matrix(a) for a in expr.args))
    if isinstance(expr,PermuteDims) and expr.permutation==(1,0):return convert_array_to_matrix(expr.expr).T
    if isinstance(expr,ArrayContraction) and isinstance(expr.expr,ArrayTensorProduct):
        factors=expr.expr.args;matrices=[a for a in factors if len(_shape(a))==2]
        if expr.groups==tuple((2*i+1,2*i+2) for i in range(len(matrices)-1)):
            return S.MatMul(*(convert_array_to_matrix(a) if _shape(a) else a for a in factors))
    if isinstance(expr,_ArrayExpr):
        i,j=S.symbols('_i _j',integer=True)
        return S.FunctionMatrix(*expr.shape,S.Lambda((i,j),expr[i,j]))
    return expr


def convert_indexed_to_array(expr,first_indices=None):
    from ....codegen.ast import Indexed
    if isinstance(expr,(Indexed,ArrayElement)):
        result=expr.name if isinstance(expr,ArrayElement) else ArraySymbol(str(expr.base.label),expr.base.shape)
        indices=tuple(expr.indices)
        if len(set(indices))!=len(indices):
            groups=[tuple(j for j,i in enumerate(indices) if i==index) for index in dict.fromkeys(indices) if indices.count(index)>1]
            result=ArrayDiagonal(result,*groups)
            remaining=[indices[j] for j in range(len(indices)) if all(j not in group for group in groups)]+[indices[g[0]] for g in groups]
        else:remaining=list(indices)
        if first_indices is not None:
            order=list(first_indices)+[i for i in remaining if i not in first_indices]
            if order!=remaining:result=PermuteDims(result,tuple(remaining.index(i) for i in order))
        return result
    try:return convert_matrix_to_array(S.MatrixExpr.from_index_summation(expr,first_index=first_indices[0] if first_indices else None))
    except (TypeError,ValueError,NotImplementedError):
        raise NotImplementedError('indexed conversion requires an Indexed array or a recognized matrix summation') from None


__all__=['ArrayAdd','ArrayContraction','ArrayDiagonal','ArrayElement','ArrayElementwiseApplyFunc','ArraySymbol',
 'ArrayTensorProduct','OneArray','PermuteDims','Reshape','ZeroArray','array_derive','convert_array_to_indexed',
 'convert_array_to_matrix','convert_indexed_to_array','convert_matrix_to_array']
