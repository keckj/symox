"""Bounded identity reuse for canonical exact monomial sums."""
from collections import OrderedDict
from . import _core

_values=OrderedDict()
_arguments=0


def add(args):
    global _arguments
    if not 2<=len(args)<=1024:
        return None
    key=_core._exact_expr_key(args)
    if key is None:
        return None
    if key in _values:
        _values.move_to_end(key)
        return _values[key][1]
    result=_core._monomial_add(args)
    if result is None:
        return None
    while len(_values)>=32 or _arguments+len(args)>4096:
        owners,_=_values.popitem(last=False)[1]
        _arguments-=len(owners)
    # Retain all argument owners: their Python object IDs cannot be recycled
    # while the entry exists. Binder presentations, floats and foreign operand
    # protocols never enter this cache.
    _values[key]=(args,result)
    _arguments+=len(args)
    return result
