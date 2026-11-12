"""Deterministic ordering with lazy, successive tie-breaking keys."""
from collections import defaultdict


def _node_count(e, memo=None):
    # this not only counts nodes, it affirms that the
    # args are Basic (i.e. have an args property).
    # The count is that of the expanded tree (SymPy), but a shared subgraph
    # is counted once and reused: ``memo`` maps each node (by equality, or
    # by identity if unhashable) to its count for the duration of the
    # outermost call. Native wrappers are fresh objects at each ``args``.
    if getattr(e, 'is_Float', False):
        return 0.5
    if isinstance(e, str):
        # The reference stores a MatrixSymbol name as a ``Str`` atom.
        return 1
    if memo is None:
        memo = {}
    try:
        key = (e, type(e))
        hit = memo.get(key)
    except TypeError:
        key = (id(e), e)  # the node is held so that its id is not reused
        hit = memo.get(key[0])
        key = key[0]
        if hit is None:
            memo.setdefault('_held', []).append(e)
    if hit is not None:
        return hit
    args = getattr(e, 'args', None)
    if args is None:
        args = _immutable_matrix_args(e)
    count = 1 + sum(_node_count(a, memo) for a in args)
    memo[key] = count
    return count


def _immutable_matrix_args(e):
    # Reference immutable matrices are Basic with args (rows, cols, entries):
    # a Tuple of the flat entries, or a Dict of the nonzero entries if sparse.
    from .. import _core, Dict, Integer, Tuple
    rows, cols = e.shape
    if isinstance(e, _core.Matrix):
        entries = Tuple(*e.flat())
    else:
        entries = Dict(e.todok())
    return (Integer(rows), Integer(cols), entries)


def _nodes(e):
    """
    A helper for ordered() which returns the node count of ``e`` which
    for Basic objects is the number of Basic nodes in the expression tree
    but for other objects is 1 (unless the object is an iterable or dict
    for which the sum of nodes is returned).
    """
    from .._basic import Basic
    from ._iterable_extensions import iterable
    if isinstance(e, Basic):
        if getattr(e, 'is_Derivative', False):
            # Native variable counts may be plain Python integers.
            return _nodes(e.expr) + sum(i[1] if isinstance(i[1], int) or i[1].is_Number
                else _nodes(i[1]) for i in e.variable_count)
        return _node_count(e)
    elif iterable(e):
        return 1 + sum(_nodes(ei) for ei in e)
    elif isinstance(e, dict):
        return 1 + sum(_nodes(k) + _nodes(v) for k, v in e.items())
    else:
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
