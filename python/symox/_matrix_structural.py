"""Argument normalization shared by native dense and sparse transformations.

Permutation algebra and elementary scalar arithmetic are native. This module
only reads Python protocols and normalizes SymPy's legacy method signatures.
"""
import itertools
import operator

_MAX_ITEMS = 100_000


def _index(value):
    try:
        return operator.index(value)
    except TypeError:
        if getattr(value, 'is_Integer', False):
            return int(value)
        raise TypeError('Matrix transformation indices must be integers') from None


def _bounded_list(values):
    result = list(itertools.islice(values, _MAX_ITEMS + 1))
    if len(result) > _MAX_ITEMS:
        raise NotImplementedError('Matrix structural transformation resource limit exceeded')
    return result


def normalize_permutation(perm, size, orientation='rows', direction='forward'):
    orientation = {'columns': 'cols'}.get(orientation, orientation)
    direction = {'forwards': 'forward', 'backwards': 'backward'}.get(direction, direction)
    if orientation not in ('rows', 'cols'):
        raise TypeError("orientation must be 'rows' or 'cols'")
    if direction not in ('forward', 'backward'):
        raise TypeError("direction must be 'forward' or 'backward'")
    # Interoperability with explicitly supplied permutation objects. Never
    # construct/import a foreign symbolic object to perform the transformation.
    if hasattr(type(perm), 'array_form') or 'array_form' in getattr(perm, '__dict__', {}):
        declared = getattr(perm, 'size', None)
        if declared is not None and _index(declared) > _MAX_ITEMS:
            raise NotImplementedError('Matrix structural transformation resource limit exceeded')
        perm = perm.array_form
    try:
        values = _bounded_list(iter(perm))
    except TypeError:
        raise ValueError('perm must be a list, a list of lists, or a permutation object') from None
    from .utilities.iterables import iterable
    cycles = bool(values) and iterable(values[0], exclude=(str, bytes))
    if cycles:
        normalized = []
        total = 0
        for cycle in values:
            row = [_index(i) for i in _bounded_list(iter(cycle))]
            total += len(row)
            if total > _MAX_ITEMS:
                raise NotImplementedError('Matrix structural transformation resource limit exceeded')
            normalized.append(row)
    else:
        normalized = [[_index(i) for i in values]]
    return normalized, cycles, direction == 'backward', orientation == 'rows'


def normalize_elementary(op, size, target=None, coefficient=None, first=None, second=None):
    if op not in ('n->kn', 'n<->m', 'n->n+km'):
        raise ValueError('Unknown elementary operation: %s' % op)
    if op == 'n->kn':
        target = target if target is not None else first
        if target is None or coefficient is None:
            raise ValueError('Scaling requires an index and a coefficient')
        target = _index(target)
        source, kind = target, 0
    elif op == 'n<->m':
        candidates = {x for x in (target, coefficient, first, second) if x is not None}
        if len(candidates) > 2:
            candidates = {x for x in (target, first, second) if x is not None}
        if len(candidates) != 2:
            raise ValueError('Swapping requires two distinct indices')
        target, source = sorted(_index(x) for x in candidates)
        coefficient, kind = 1, 1
    else:
        target = target if target is not None else first
        source = second if second is not None else first
        if target is None or source is None or coefficient is None:
            raise ValueError('Adding a multiple requires two indices and a coefficient')
        target, source = _index(target), _index(source)
        if target == source:
            raise ValueError('Elementary operation indices must be distinct')
        kind = 2
    if not (0 <= target < size and 0 <= source < size):
        raise ValueError('Elementary operation index out of range')
    return kind, target, source, coefficient
