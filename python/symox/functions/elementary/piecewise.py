from ... import Piecewise,piecewise_fold,piecewise_exclusive
__all__ = ['Piecewise','piecewise_fold','piecewise_exclusive']

from ..._matrixexpr_methods import _ExprCondPairView as ExprCondPair
__all__.append('ExprCondPair')

from ... import S as _S
# the reference alias of nan for the holes of a Piecewise
Undefined = _S.NaN
__all__.append('Undefined')

from ..._public_calculus import (piecewise_simplify, piecewise_simplify_arguments,
                                 _piecewise_simplify_eq_and, _piecewise_simplify_equal_to_next_segment,
                                 _blessed)
__all__ += ['piecewise_simplify', 'piecewise_simplify_arguments']
