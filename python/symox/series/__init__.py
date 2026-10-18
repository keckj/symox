"""Native finite, formal and Fourier series."""
from .formal import FormalPowerSeries, fps
from .fourier import FourierSeries, FiniteFourierSeries, fourier_series
from .sequences import SeqFormula, SeqPer, SeqAdd, SeqMul, EmptySequence, sequence


def series(expr, x=None, x0=0, n=6, dir='+'):
    from .. import _series_call
    return _series_call(expr, x, x0, n, dir)


__all__ = ['series', 'fps', 'FormalPowerSeries', 'fourier_series', 'FourierSeries',
    'FiniteFourierSeries','SeqFormula','SeqPer','SeqAdd','SeqMul','EmptySequence','sequence']
