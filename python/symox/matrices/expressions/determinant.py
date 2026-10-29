from ..._matrixexpr import Determinant
# SymPy's module also defines det, per and Permanent; symox's are the
# top-level ones.
from ... import det, per, Permanent

__all__ = ['Determinant', 'det', 'Permanent', 'per']
