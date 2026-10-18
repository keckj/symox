"""Native finite, formal and Fourier series."""
from .formal import FormalPowerSeries, fps
from .fourier import FourierSeries, FiniteFourierSeries, fourier_series
from .sequences import (SeqFormula, SeqPer, SeqAdd, SeqMul, EmptySequence, sequence,
                        SeqExpr,SeqExprOp,RecursiveSeq)
from .limits import Limit,limit
from .order import O,Order
from .residues import residue
from .approximants import approximants
from .limitseq import difference_delta,limit_seq
from .gruntz import gruntz


def series(expr, x=None, x0=0, n=6, dir='+'):
    from .. import _series_call
    return _series_call(expr, x, x0, n, dir)


__all__ = ['series', 'fps', 'FormalPowerSeries', 'fourier_series', 'FourierSeries',
    'FiniteFourierSeries','SeqFormula','SeqPer','SeqAdd','SeqMul','EmptySequence','sequence']
__all__ += ['SeqExpr','SeqExprOp','RecursiveSeq']
__all__ += ['Limit','limit','O','Order','residue','approximants','difference_delta','limit_seq','gruntz']
