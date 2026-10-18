"""Base and derived physical dimensions."""
from .. import _dimensions,Dimension
for _name,_value in _dimensions.items():globals()[_name] = _value
angle = Dimension('angle')
__all__ = list(_dimensions)+['angle']
