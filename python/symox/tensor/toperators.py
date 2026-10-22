"""Tensor partial derivatives with explicit component differentiation."""
from ._expressions import TensExpr,Tensor
from ._indexing import _IndexStructure


def _contract_indices_for_derivative(expr,variables):
    """PartialDerivative._contract_indices_for_derivative: an index of the
    expression contracted with a variable slot is a dummy of the derivative and
    takes the next dummy name of its type, as the reference TensMul
    contraction with ``replace_indices=True``."""
    from ._expressions import TensMul,_normalize
    from ._indexing import TensorIndex,_name_key
    if not isinstance(expr,(Tensor,TensMul)) or not all(isinstance(v,Tensor) for v in variables):
        return expr,variables
    own=list(expr.get_indices())
    flipped=[-i for v in variables for i in v.get_indices()]
    try:structure=_IndexStructure(own+flipped)
    except ValueError:
        seen=set()
        for index in own+flipped:
            if index in seen:raise ValueError('Repeated index: %s' % index)
            seen.add(index)
        raise
    if not any(upper>=len(own) or lower>=len(own) for upper,lower in structure.dum):
        return expr,variables
    indices=structure.indices
    occupied={_name_key(index) for index,_ in structure.free};counters={};replacements={}
    for upper,lower in structure.dum:
        index=indices[upper];typ=index.tensor_index_type
        counter=counters.get(typ,0)
        while (typ,typ.dummy_name+'_'+str(counter)) in occupied:counter+=1
        new=TensorIndex(typ.dummy_name+'_'+str(counter),typ);counters[typ]=counter+1
        replacements[_name_key(index)]=new
    def rename(index):
        new=replacements.get(_name_key(index))
        return index if new is None else new if index.is_up else -new
    expr=expr._set_indices(*(rename(i) for i in own))
    variables=tuple(Tensor._raw(v.component,tuple(rename(i) for i in v.get_indices()),type(v)) for v in variables)
    return expr,variables


class PartialDerivative(TensExpr):
    __slots__=()
    def __new__(cls,expr,*variables):
        if isinstance(expr,PartialDerivative):variables=expr.variables+variables;expr=expr.expr
        if not variables:raise ValueError('a partial derivative requires a variable')
        expr,variables=_contract_indices_for_derivative(expr,variables)
        result=super().__new__(cls,expr,*variables)
        _IndexStructure(result.get_indices())
        return result
    expr=property(lambda self:self.args[0])
    variables=property(lambda self:self.args[1:])
    @property
    def coeff(self):
        from symox import S
        return S.One
    nocoeff=property(lambda self:self)
    def get_indices(self):
        indices=list(self.expr.get_indices()) if isinstance(self.expr,TensExpr) else []
        for variable in self.variables:
            if isinstance(variable,TensExpr):indices.extend(-index for index in variable.get_free_indices())
        return indices
    def __repr__(self):return 'PartialDerivative('+', '.join(map(str,self.args))+')'
    def _extract_data(self,replacements):
        from ._components import _extract
        from .array import Array,NDimArray,derive_by_array,permutedims,tensorcontraction
        if isinstance(self.expr,TensExpr):indices,array=_extract(self.expr,replacements)
        else:indices,array=[],self.expr
        for variable in self.variables:
            if isinstance(variable,TensExpr):
                var_indices,values=_extract(variable,replacements)
                old_rank=len(indices);rank=len(var_indices)
                array=derive_by_array(array,values)
                array=permutedims(array,tuple(range(rank,rank+old_rank))+tuple(range(rank)))
                indices.extend(-index for index in var_indices)
            else:
                array=derive_by_array(array,variable) if isinstance(array,NDimArray) else array.diff(variable)
        structure=_IndexStructure(indices)
        if structure.dum:array=tensorcontraction(array,*structure.dum)
        return structure.get_free_indices(),array
    def replace_with_arrays(self,replacement_dict,indices=None):
        from ._components import _align,_remap_indices
        from .array import NDimArray
        indices=_remap_indices(self,indices)
        free,array=self._extract_data(dict(replacement_dict))
        _,array=_align(free,array,free if indices is None else list(indices),replacement_dict)
        return array[()] if isinstance(array,NDimArray) and not array.rank() else array

__all__=['PartialDerivative']
