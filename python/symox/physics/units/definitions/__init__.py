"""Physical quantities and constants from the reference unit catalog."""
from .. import _quantities,ALIASES
for _name,(_kind,_key) in ALIASES.items():
    if _kind == 'quantity':globals()[_name] = _quantities[_key]
__all__ = [n for n,(kind,key) in ALIASES.items() if kind == 'quantity']
