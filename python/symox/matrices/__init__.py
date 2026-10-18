"""Dense, sparse and symbolic matrix APIs."""
from .. import (Matrix, MutableDenseMatrix, ImmutableMatrix, ImmutableDenseMatrix,
                SparseMatrix, MutableSparseMatrix, ImmutableSparseMatrix,
                eye, zeros, ones, diag, hessian, wronskian, jacobian)
from .._matrixexpr import *
from .._matrix_permutations import PermutationMatrix, MatrixPermute

# Compatible imports reuse the existing native objects.
from symox._api_exports import export as _export
_export(globals(), 'symox', 'DeferredVector GramSchmidt MatrixBase MutableMatrix NonSquareMatrixError Permanent ShapeError banded blockcut casoratian det diagonalize_vector jordan_cell list2numpy matrix2numpy matrix_multiply_elementwise matrix_symbols per randMatrix rot_axis1 rot_axis2 rot_axis3 rot_ccw_axis1 rot_ccw_axis2 rot_ccw_axis3 rot_givens symarray trace')
from .kind import MatrixKind
from .expressions.sets import MatrixSet


def dotprodsimp(expr,withsimp=False):
    from symox import cancel,expand,sympify
    original=sympify(expr);result=cancel(expand(original))
    return (result,result!=original) if withsimp else result

__all__ += ['MatrixKind','MatrixSet','dotprodsimp']
