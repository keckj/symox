"""Native finite, formal and Fourier series."""
# The formal series classes are host expression classes: they are created
# once symox itself has finished importing (``symox/__init__``), not here.
_FORMAL = ('FormalPowerSeries', 'FiniteFormalPowerSeries', 'FormalPowerSeriesProduct',
           'FormalPowerSeriesCompose', 'FormalPowerSeriesInverse')


def fps(f, x=None, x0=0, dir=1, hyper=True, order=4, rational=True, full=False):
    from .formal import fps as _fps
    return _fps(f, x, x0, dir, hyper, order, rational, full)


def __getattr__(name):
    if name in _FORMAL:
        from . import formal
        return getattr(formal, name)
    raise AttributeError(name)
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
