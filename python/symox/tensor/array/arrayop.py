"""Concrete tensor operations using native scalar arithmetic."""
import itertools
import math
import operator

from ... import _core
from ._array import (Array, NDimArray, _array_class, _scalar, _shape_tuple,
                     _ravel, _unravel, MAX_ENTRIES)


def _array(value):
    return value if isinstance(value,NDimArray) else Array(value)


def _tensor_argument(value):
    if isinstance(value,NDimArray) or hasattr(value,'shape'):
        return _array(value)
    if isinstance(value,(tuple,list)) or (isinstance(value,_core.Expr) and
            getattr(value.func,'__name__','')=='Tuple'):
        return _array(value)
    return _scalar(value)


def _expressions():
    from . import expressions
    return expressions


def tensorproduct(*arguments):
    from ... import S
    if len(arguments)>1 and any(_expressions()._is_symbolic_array(arg) for arg in arguments):
        # as SymPy, a matrix symbol or array expression keeps the product held
        return _expressions().ArrayTensorProduct(*arguments)
    if len(arguments)==1 and _expressions()._is_symbolic_array(arguments[0]):
        return arguments[0]
    values=[_tensor_argument(arg) for arg in arguments]
    arrays=[value for value in values if isinstance(value,NDimArray)]
    if not arrays:
        result=S.One
        for value in values:
            result *= value
        return result
    shape=_shape_tuple(tuple(dimension for array in arrays for dimension in array.shape))
    sparse=all(array._sparse for array in arrays)
    size=math.prod(shape)
    if size==0:
        return _array_class(sparse)({},shape)
    # Sparse multiplication is valid only while zero annihilates every
    # factor. For example 0*oo is nan and must materialize those positions.
    if sparse:
        safe=True
        data={0:S.One}
        for value in values:
            items=list(value._items()) if isinstance(value,NDimArray) else [(0,value)]
            if any(S.Zero*v!=S.Zero or v*S.Zero!=S.Zero for _,v in items):
                safe=False
                break
            if any(S.Zero*v!=S.Zero or v*S.Zero!=S.Zero for v in data.values()):
                safe=False
                break
            if len(data)*len(items)>MAX_ENTRIES:
                raise ValueError('Tensor product exceeds 100000 stored intermediate entries')
            width=len(value) if isinstance(value,NDimArray) else 1
            data={i*width+j:old*new for i,old in data.items() for j,new in items}
        if safe:
            return _array_class(True)(data,shape)
    if size>MAX_ENTRIES:
        raise ValueError('Tensor product materialization exceeds 100000 entries')
    factors=[list(value._flat_values()) if isinstance(value,NDimArray) else [value] for value in values]
    def products():
        for factors_at_index in itertools.product(*factors):
            result=S.One
            for factor in factors_at_index:
                result *= factor
            yield result
    return _array_class(sparse)(products(),shape)


def _groups(array, groups):
    used=set()
    normalized=[]
    for group in groups:
        try:
            group=tuple(operator.index(axis) for axis in group)
        except TypeError:
            raise ValueError('Tensor axes must be groups of integers') from None
        if not group:
            raise ValueError('A tensor axis group must be nonempty')
        if any(not 0<=axis<array.rank() for axis in group):
            raise ValueError('Tensor axis is out of range')
        if len(set(group))!=len(group) or used.intersection(group):
            raise ValueError('A tensor axis may appear in only one group')
        if len({array.shape[axis] for axis in group})!=1:
            raise ValueError('Axes in a tensor group must have equal dimensions')
        used.update(group)
        normalized.append(group)
    remaining=tuple(axis for axis in range(array.rank()) if axis not in used)
    return tuple(normalized), remaining


def tensorcontraction(array, *contraction_axes):
    from ... import S
    if _expressions()._is_symbolic_array(array):
        if not contraction_axes:
            return array
        return _expressions()._array_contraction(array,*contraction_axes)
    array=_array(array)
    if not contraction_axes:
        return array
    groups,remaining=_groups(array,contraction_axes)
    shape=tuple(array.shape[axis] for axis in remaining)
    data={}
    for index,value in array._items():
        position=_unravel(index,array.shape)
        if all(all(position[axis]==position[group[0]] for axis in group[1:]) for group in groups):
            new_index=_ravel(tuple(position[axis] for axis in remaining),shape)
            data[new_index]=data.get(new_index,S.Zero)+value
    if not shape:
        return data.get(0,S.Zero)
    return type(array)(data,shape)


def tensordiagonal(array, *diagonal_axes):
    if _expressions()._is_symbolic_array(array):
        if not diagonal_axes:
            return array
        return _expressions().ArrayDiagonal(array,*diagonal_axes)
    array=_array(array)
    if not diagonal_axes:
        return array
    groups,remaining=_groups(array,diagonal_axes)
    if any(len(group)<2 for group in groups):
        raise ValueError('At least two axes are required to diagonalize')
    shape=tuple(array.shape[axis] for axis in remaining)+tuple(array.shape[group[0]] for group in groups)
    data={}
    for index,value in array._items():
        position=_unravel(index,array.shape)
        if all(all(position[axis]==position[group[0]] for axis in group[1:]) for group in groups):
            new_position=tuple(position[axis] for axis in remaining)+tuple(position[group[0]] for group in groups)
            data[_ravel(new_position,shape)]=value
    return type(array)(data,shape)


def permutedims(expr, perm=None, index_order_old=None, index_order_new=None):
    if _expressions()._is_symbolic_array(expr):
        if index_order_old is not None or index_order_new is not None:
            if perm is not None or index_order_old is None or index_order_new is None:
                raise ValueError('Specify a permutation or both index orders')
            old,new=tuple(index_order_old),tuple(index_order_new)
            perm=tuple(old.index(index) for index in new)
        if perm is None:
            raise ValueError('Permutation is not defined')
        return _expressions()._permute_dims(expr,perm)
    expr=_array(expr)
    if index_order_old is not None or index_order_new is not None:
        if perm is not None or index_order_old is None or index_order_new is None:
            raise ValueError('Specify a permutation or both index orders')
        old,new=tuple(index_order_old),tuple(index_order_new)
        if len(old)!=expr.rank() or len(set(old))!=len(old) or set(old)!=set(new) or len(old)!=len(new):
            raise ValueError('Index orders must be bijections of the array axes')
        perm=tuple(old.index(index) for index in new)
    if hasattr(perm,'array_form'):
        perm=perm.array_form
    if perm is None:
        raise ValueError('Permutation is not defined')
    perm=tuple(operator.index(axis) for axis in perm)
    if sorted(perm)!=list(range(expr.rank())):
        raise ValueError('Permutation must contain each array axis once')
    shape=tuple(expr.shape[axis] for axis in perm)
    data={}
    for index,value in expr._items():
        position=_unravel(index,expr.shape)
        data[_ravel(tuple(position[axis] for axis in perm),shape)]=value
    return type(expr)(data,shape)


def derive_by_array(expr, dx):
    from ... import S
    expr=_tensor_argument(expr)
    dx=_tensor_argument(dx)
    if not isinstance(dx,NDimArray):
        dx=_validate_variable(dx)
        return (_derivative(expr,dx) if not isinstance(expr,NDimArray)
                else expr.applyfunc(lambda value:_derivative(value,dx)).as_immutable())
    dx._materialize_limit()
    entries=list(dx._flat_values())
    for variable in entries:
        # Even an all-zero source must reject an invalid differentiation variable.
        _validate_variable(variable)
    if not isinstance(expr,NDimArray):
        return Array((_derivative(expr,variable) for variable in entries),dx.shape)
    shape=_shape_tuple(dx.shape+expr.shape)
    sparse=expr._sparse
    if not sparse and math.prod(shape)>MAX_ENTRIES:
        raise ValueError('Array derivative exceeds 100000 entries')
    data={}
    for i,variable in enumerate(entries):
        for j,value in expr._items():
            derivative=_derivative(value,variable)
            if derivative!=S.Zero:
                data[i*len(expr)+j]=derivative
                if len(data)>MAX_ENTRIES:
                    raise ValueError('Array derivative exceeds 100000 stored entries')
    # SymPy builds the result with the type of the immutable source.
    return type(expr.as_immutable())(data,shape)


def _validate_variable(variable):
    variable=_scalar(variable)
    if (getattr(variable,'is_Symbol',False) or getattr(variable,'is_Derivative',False)
            or getattr(variable,'is_UndefinedFunction',False)
            or getattr(variable,'is_Indexed',False)
            or getattr(variable.func,'__name__','') in ('MatrixElement','Indexed')):
        return variable
    raise ValueError('Cannot differentiate with respect to %s' % variable)


def _derivative(expression,variable):
    if getattr(variable.func,'__name__','')=='MatrixElement' and not getattr(variable,'is_UndefinedFunction',False):
        expression=_expand_matrix_scalars(expression)
    return expression.diff(variable)


def _expand_matrix_scalars(expression):
    """Expose finite scalar matrix functions before entry differentiation."""
    if not isinstance(expression,_core.Expr):
        return expression
    for _ in range(16):
        pending=[expression]
        seen=set()
        replacements={}
        while pending:
            value=pending.pop()
            if value in seen:
                continue
            seen.add(value)
            arguments=value.args
            name=getattr(value.func,'__name__','')
            if getattr(value,'is_UndefinedFunction',False):
                name=''
            try:
                if name in ('Trace','Determinant','DotProduct'):
                    explicit=[arg.as_explicit() if hasattr(arg,'as_explicit') else arg for arg in arguments]
                    expanded=value.func(*explicit).doit()
                elif name=='MatrixElement':
                    expanded=arguments[0][arguments[1],arguments[2]]
                else:
                    expanded=value
            except (NotImplementedError,ValueError):
                expanded=value
            if expanded!=value:
                replacements[value]=expanded
            pending.extend(arg for arg in arguments if isinstance(arg,_core.Expr))
        if not replacements:
            return expression
        expression=expression.xreplace(replacements)
    return expression


class Flatten:
    """SymPy's ``Flatten``: a lazy iterator over the flat entries of an
    array, a matrix or another iterable.

    >>> from symox.tensor.array.arrayop import Flatten
    >>> from symox.tensor.array import Array
    >>> A = Array(range(6)).reshape(2, 3)
    >>> Flatten(A)
    Flatten([[0, 1, 2], [3, 4, 5]])
    >>> [i for i in Flatten(A)]
    [0, 1, 2, 3, 4, 5]
    """
    def __init__(self, iterable):
        from collections.abc import Iterable
        from ... import _core
        is_matrix=isinstance(iterable,(_core.Matrix,)) or getattr(iterable,'is_Matrix',False)
        if not isinstance(iterable,Iterable) and not is_matrix:
            raise NotImplementedError("Data type not yet supported")
        if isinstance(iterable,list):
            iterable=NDimArray(iterable)
        self._iter=iterable
        self._idx=0

    def __iter__(self):
        return self

    def __next__(self):
        from ... import S
        if len(self._iter)>self._idx:
            source=self._iter
            if isinstance(source,NDimArray):
                result=(source._sparse_array.get(self._idx,S.Zero) if source._sparse
                        else source._array[self._idx])
            elif hasattr(source,'__next__'):
                result=next(source)
            else:
                result=source[self._idx]
        else:
            raise StopIteration
        self._idx+=1
        return result

    def next(self):
        return self.__next__()

    def _sympystr(self, printer):
        # SymPy's StrPrinter prints an NDimArray as its nested list.
        inner=str(self._iter) if isinstance(self._iter,NDimArray) else printer._print(self._iter)
        return type(self).__name__+'('+inner+')'

    # SymPy's Printable (imported lazily: the printing package needs the
    # complete symox namespace).
    def __str__(self):
        from ...printing.str import sstr
        return sstr(self, order=None)

    __repr__ = __str__


__all__=['tensorproduct','tensorcontraction','tensordiagonal','permutedims','derive_by_array','Flatten']
