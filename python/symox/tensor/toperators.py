"""Tensor partial derivatives with explicit component differentiation."""
from ._expressions import TensExpr,Tensor
from ._indexing import _IndexStructure


class PartialDerivative(TensExpr):
    __slots__=()
    def __new__(cls,expr,*variables):
        if isinstance(expr,PartialDerivative):variables=expr.variables+variables;expr=expr.expr
        if not variables:raise ValueError('a partial derivative requires a variable')
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
        from ._components import _align
        from .array import NDimArray
        free,array=self._extract_data(dict(replacement_dict))
        _,array=_align(free,array,free if indices is None else list(indices),replacement_dict)
        return array[()] if isinstance(array,NDimArray) and not array.rank() else array

__all__=['PartialDerivative']
