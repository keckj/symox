"""Symbolic probability over native SymOx expressions.

Distribution formulas, finite probability spaces, moments, joint models and
sampling are implemented locally. General integrals and sums may remain exact
unevaluated symbolic expressions. See docs/statistics.md for the supported scope.
"""
from ._univariate import *
from ._queries import *
from ._joint import *
from ._matrix import *
from ._processes import *
from ._ensembles import *
from ._symbolic import *
from . import _univariate, _queries, _joint, _matrix, _processes, _ensembles, _symbolic

__all__ = list(dict.fromkeys(_univariate.__all__ + _queries.__all__ + _joint.__all__ +
    _matrix.__all__ + _processes.__all__ + _ensembles.__all__ + _symbolic.__all__))
_symbolic.install()
