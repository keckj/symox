from ... import Piecewise,piecewise_fold,piecewise_exclusive
__all__ = ['Piecewise','piecewise_fold','piecewise_exclusive']

from ..._matrixexpr_methods import _ExprCondPairView as ExprCondPair
__all__.append('ExprCondPair')

from ... import S as _S
# the reference alias of nan for the holes of a Piecewise
Undefined = _S.NaN
__all__.append('Undefined')
