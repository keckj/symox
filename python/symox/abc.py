"""Predefined Latin and Greek symbols for interactive symbolic calculations."""
from string import ascii_letters as _letters
from . import Symbol as _Symbol

_greek = ('alpha beta gamma delta epsilon zeta eta theta iota kappa lamda mu '
          'nu xi omicron pi rho sigma tau upsilon phi chi psi omega').split()
__all__ = list(_letters) + _greek
for _name in __all__:
    globals()[_name] = _Symbol(_name)

# sympify's locals convention treats an empty string as a request for a Symbol.
_clash1 = dict.fromkeys(('O', 'Q', 'N', 'I', 'E', 'S'), '')
_clash2 = dict.fromkeys(('beta', 'zeta', 'gamma', 'pi'), '')
_clash = {**_clash1, **_clash2}
