"""``sympy.tensor.array.expressions.array_expressions`` import path.

Re-exports the array-expression classes of ``symox.tensor.array.expressions``
(no computation here). SymPy's private canonicalizing constructors and the
``_EditArrayContraction`` machinery are not ported (docs/algorithmic-parity.md,
"arrexprrust").
"""
from . import (ArrayAdd, ArrayContraction, ArrayDiagonal, ArrayElement,
               ArrayElementwiseApplyFunc, ArraySymbol, ArrayTensorProduct,
               OneArray, PermuteDims, Reshape, ZeroArray, _ArrayExpr)

__all__ = ['ArrayAdd', 'ArrayContraction', 'ArrayDiagonal', 'ArrayElement',
           'ArrayElementwiseApplyFunc', 'ArraySymbol', 'ArrayTensorProduct',
           'OneArray', 'PermuteDims', 'Reshape', 'ZeroArray']
