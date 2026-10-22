"""Evaluate indexed tensors using native dense arrays and exact metrics."""
import itertools
import math

from ._indexing import TensorIndex,_IndexStructure,_name_key
from ._expressions import TensorHead,Tensor,TensExpr,TensAdd,TensMul,TensorElement,_split


def _array(value):
    from .array import NDimArray,Array
    return value if isinstance(value,NDimArray) else Array(value)


def _metric(typ,replacements):
    from .. import Matrix,cancel
    if typ not in replacements:raise ValueError('Missing covariant metric components for index type '+typ.name)
    value=_array(replacements[typ])
    if value.rank()!=2 or value.shape[0]!=value.shape[1]:raise ValueError('Metric components must form a square matrix')
    value._materialize_limit()
    if typ.args[4]!=1:raise NotImplementedError('Automatic component raising/lowering requires a symmetric metric')
    if any(cancel(value[i,j]-value[j,i])!=0 for i in range(value.shape[0]) for j in range(i)):
        raise ValueError('A symmetric tensor metric requires symmetric component data')
    return value.tomatrix()


def _change_axis(array,axis,typ,want_up,replacements):
    from .. import S
    from .array import Array
    array._materialize_limit()
    metric=_metric(typ,replacements)
    if want_up:metric=metric.inv()
    if array.shape[axis]!=metric.rows:raise ValueError('Metric dimension does not match the array axis')
    values=[]
    for point in itertools.product(*(range(size) for size in array.shape)):
        row=point[axis]
        values.append(sum((metric[row,k]*array[point[:axis]+(k,)+point[axis+1:]]
                           for k in range(metric.cols)),S.Zero))
    return Array(values,array.shape)


def _align(indices,array,wanted,replacements):
    from .array import permutedims
    if len(indices)!=len(wanted) or {_name_key(i) for i in indices}!={_name_key(i) for i in wanted}:
        raise ValueError('Requested component axes must be the tensor free indices')
    if not indices:return [],array
    for axis,index in enumerate(indices):
        target=next(i for i in wanted if _name_key(i)==_name_key(index))
        if index.is_up!=target.is_up:array=_change_axis(array,axis,index.tensor_index_type,target.is_up,replacements)
    order=tuple(next(position for position,index in enumerate(indices) if _name_key(index)==_name_key(target)) for target in wanted)
    return list(wanted),permutedims(array,order)


def _tensor_data(tensor,replacements):
    from .. import eye,cancel
    from .array import Array,tensorcontraction,permutedims
    head=tensor.component;source=None;source_indices=None
    if head.comm==1:raise NotImplementedError('Anticommuting components require a Grassmann coefficient algebra')
    for key,value in replacements.items():
        if isinstance(key,Tensor) and key.component==head:
            if key.dum:continue
            source=_array(value);source_indices=list(key.indices);break
        if isinstance(key,TensorHead) and key==head:
            source=_array(value);source_indices=[TensorIndex('_array_axis_'+str(i),typ) for i,typ in enumerate(head.index_types)];break
    if source is None:
        typ=head.index_types[0] if head.rank==2 and head.index_types[0]==head.index_types[1] else None
        if typ is not None and head==typ.metric:
            source=_array(_metric(typ,replacements));source_indices=[TensorIndex('_metric_axis_'+str(i),typ,False) for i in range(2)]
        elif typ is not None and head==typ.delta:
            dimension=int(typ.dim) if typ.dim.is_Integer else _metric(typ,replacements).rows
            source=Array(eye(dimension));source_indices=[TensorIndex('_delta0',typ),TensorIndex('_delta1',typ,False)]
        else:raise ValueError('Missing array components for tensor head '+head.name)
    if source.rank()!=head.rank:raise ValueError('Array rank does not match the tensor head')
    for generator in head.symmetry.generators:
        order=tuple(generator.array_form[:-2])
        if order==tuple(range(head.rank)):continue
        if any(source_indices[i].is_up!=source_indices[j].is_up for i,j in enumerate(order)):continue
        permuted=permutedims(source,order)
        sign=-1 if generator.array_form[-1]!=head.rank+1 else 1
        stored={index for index,_ in source._items()}|{index for index,_ in permuted._items()}
        if any(cancel(source._get(index)-sign*permuted._get(index))!=0 for index in stored):
            raise ValueError('Component data violate the declared tensor slot symmetry')
    for axis,(typ,index,given) in enumerate(zip(head.index_types,tensor.indices,source_indices)):
        if typ.dim.is_Integer and source.shape[axis]!=int(typ.dim):raise ValueError('Array dimension does not match its index type')
        if index.is_up!=given.is_up:source=_change_axis(source,axis,typ,index.is_up,replacements)
    structure=_IndexStructure(tensor.indices)
    if structure.dum:source=tensorcontraction(source,*structure.dum)
    return structure.get_free_indices(),source


def _extract(expression,replacements):
    from .array import tensorproduct,tensorcontraction
    if isinstance(expression,Tensor):return _tensor_data(expression,replacements)
    if isinstance(expression,TensorElement):
        indices,array=_extract(expression.expr,replacements)
        positions=expression.index_map
        if any(not value.is_Integer for value in positions.values()):raise NotImplementedError('Selected tensor components require concrete integer positions')
        selectors=tuple(int(positions[index]) if index in positions else slice(None) for index in indices)
        return [i for i in indices if i not in positions],array[selectors]
    if isinstance(expression,TensMul):
        indices=[];arrays=[]
        for factor in expression._factors:
            free,array=_extract(factor,replacements);indices.extend(free);arrays.append(array)
        source=tensorproduct(*arrays)
        structure=_IndexStructure(indices)
        if structure.dum:source=tensorcontraction(source,*structure.dum)
        return structure.get_free_indices(),source*expression.coeff
    if isinstance(expression,TensAdd):
        wanted=expression.get_free_indices();result=None
        for term in expression.args:
            if isinstance(term,TensExpr):indices,array=_extract(term,replacements)
            else:indices,array=[],term
            _,array=_align(indices,array,wanted,replacements)
            result=array if result is None else result+array
        return wanted,result
    raise NotImplementedError('No component evaluator for this tensor expression')


def _remap_indices(expression,indices):
    """TensExpr.replace_with_arrays: symbols ``i`` and ``-i`` in the requested
    order stand for the free tensor indices of the same name and variance.
    The reference rewrites the caller's list in place."""
    if not indices:return indices
    remap=None
    for position,index in enumerate(indices):
        if getattr(index,'is_Symbol',False) or getattr(index,'is_Mul',False):
            if remap is None:
                from .. import Symbol
                remap={}
                for k in expression.get_free_indices():
                    name=k.args[0] if hasattr(k.args[0],'is_Symbol') else Symbol(k.args[0])
                    remap[name if k.is_up else -name]=k
            indices[position]=remap[index] if index in remap else -remap[-index]
    return indices


def replace_with_arrays(expression,replacement_dict,indices=None):
    replacements=dict(replacement_dict)
    indices=_remap_indices(expression,indices)
    free,array=_extract(expression,replacements)
    wanted=free if indices is None else list(indices)
    if any(not isinstance(index,TensorIndex) for index in wanted):raise TypeError('Component axes require typed tensor indices')
    _,array=_align(free,array,wanted,replacements)
    from .array import NDimArray
    return array[()] if isinstance(array,NDimArray) and array.rank()==0 else array
