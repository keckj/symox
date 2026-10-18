"""Deterministic ordering with lazy, successive tie-breaking keys."""
from collections import defaultdict
from collections.abc import Mapping


def _nodes(value):
    if getattr(value,'is_Derivative',False):
        return _nodes(value.expr)+sum(count for _,count in value.variable_count)
    if isinstance(value,Mapping):
        return 1+sum(_nodes(k)+_nodes(v) for k,v in value.items())
    if hasattr(value,'args'):
        return 1+sum(_nodes(child) for child in value.args)
    if not isinstance(value,(str,bytes)):
        try:
            return 1+sum(_nodes(child) for child in value)
        except TypeError:
            pass
    return 1


def ordered(seq, keys=None, default=True, warn=False):
    from .. import default_sort_key
    keys=list(keys) if isinstance(keys,(tuple,list)) else [keys] if keys else []
    groups=defaultdict(list)
    if keys:
        key=keys.pop(0)
        for value in seq:
            groups[key(value)].append(value)
    else:
        if not default:
            raise ValueError('if default=False then keys must be provided')
        groups[None].extend(seq)
    for _,values in sorted(groups.items()):
        if len(values)>1:
            if keys:
                yield from ordered(values,keys,default,warn)
                continue
            if default:
                yield from ordered(values,(_nodes,default_sort_key),False,warn)
                continue
            if warn and any(value!=values[0] for value in values[1:]):
                raise ValueError('not enough keys to break ties: %s' % values)
        yield from values


__all__=['ordered']
