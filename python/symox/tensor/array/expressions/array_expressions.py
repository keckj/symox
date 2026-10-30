"""``sympy.tensor.array.expressions.array_expressions`` import path.

Re-exports the array-expression classes of ``symox.tensor.array.expressions``
and SymPy's canonicalizing constructors, whose ``_canonicalize`` rules run in
Rust (``_canon``). The ``_EditArrayContraction`` machinery is not ported
(docs/algorithmic-parity.md, "arrexprrust").
"""
from . import (ArrayAdd, ArrayContraction, ArrayDiagonal, ArrayElement,
               ArrayElementwiseApplyFunc, ArraySymbol, ArrayTensorProduct,
               OneArray, PermuteDims, Reshape, ZeroArray, _ArrayExpr)
from ._canon import (_array_add, _array_contraction, _array_diagonal,
                     _array_tensor_product, _permute_dims, nest_permutation)



class _ArgE:
    """SymPy's ``_ArgE`` (not ported: the mutable contraction editor is a gap)."""
    def __init__(self, *args, **kwargs):
        raise NotImplementedError("_ArgE: the _EditArrayContraction editor is not ported")


class _EditArrayContraction:
    """SymPy's ``_EditArrayContraction`` (not ported, see the module docstring)."""
    def __init__(self, *args, **kwargs):
        raise NotImplementedError("_EditArrayContraction is not ported")


__all__ = ['ArrayAdd', 'ArrayContraction', 'ArrayDiagonal', 'ArrayElement',
           'ArrayElementwiseApplyFunc', 'ArraySymbol', 'ArrayTensorProduct',
           'OneArray', 'PermuteDims', 'Reshape', 'ZeroArray', 'nest_permutation']
