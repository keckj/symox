"""Symbolic matrix expression constructors."""
from ..._matrixexpr import *
from ..._matrixexpr import __all__
from ..._matrix_permutations import PermutationMatrix, MatrixPermute

__all__ = [*__all__, 'PermutationMatrix', 'MatrixPermute']

# Compatible imports reuse the existing native objects.
from symox._api_exports import export as _export
_export(globals(), 'symox', 'Permanent blockcut det diagonalize_vector matrix_symbols per trace')
from .companion import CompanionMatrix
from .sets import MatrixSet
from .kronecker import combine_kronecker
__all__ += ['CompanionMatrix','MatrixSet','combine_kronecker']
