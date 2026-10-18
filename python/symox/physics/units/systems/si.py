"""International System of Units."""
from .. import SI
from .._core import DimensionSystem
dimsys_SI = SI.get_dimension_system()
__all__ = ['SI','dimsys_SI']
