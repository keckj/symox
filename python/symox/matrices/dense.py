"""Dense matrix names (SymPy's ``sympy.matrices.dense``)."""
from .. import Matrix, MutableDenseMatrix, eye, zeros, ones, diag
from .. import (GramSchmidt, casoratian, hessian, jordan_cell, list2numpy,
                matrix2numpy, matrix_multiply_elementwise, randMatrix,
                rot_axis1, rot_axis2, rot_axis3, rot_ccw_axis1, rot_ccw_axis2,
                rot_ccw_axis3, rot_givens, symarray, wronskian)
from .. import DenseMatrix
MutableMatrix = MutableDenseMatrix
