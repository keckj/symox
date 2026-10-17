"""Sparse adapters for the native structural transformation kernels.

Only finite permutation support and stored entries are visited. A transform
that fills implicit zeros must pass the native allocation bound first.
"""
from ._core import Matrix, sympify
from ._matrix_structural import normalize_permutation, normalize_elementary

_ZERO = sympify(0)


def permute(matrix, perm, orientation='rows', direction='forward'):
    rows = orientation == 'rows'
    size = matrix.rows if rows else matrix.cols
    spec, cycles, backward, rows = normalize_permutation(perm, size, orientation, direction)
    size = matrix.rows if rows else matrix.cols
    mapping = Matrix._structural_plan(spec, size, cycles, backward)
    Matrix._structural_limit(len(matrix._data))
    data = {(mapping.get(i, i), j) if rows else (i, mapping.get(j, j)): value
            for (i, j), value in matrix._data.items()}
    return matrix._new(matrix.rows, matrix.cols, data)


def elementary(matrix, rows, op, target, coefficient, first, second):
    size, width = (matrix.rows, matrix.cols) if rows else (matrix.cols, matrix.rows)
    kind, target, source, coefficient = normalize_elementary(op, size, target, coefficient, first, second)
    Matrix._structural_limit(len(matrix._data))
    if kind == 1:
        return permute(matrix, [[target, source]], 'rows' if rows else 'cols')
    if width == 0:
        return matrix._new(matrix.rows, matrix.cols, matrix._data)
    coefficient = sympify(coefficient)
    zero = Matrix._structural_values(kind == 2, coefficient, [_ZERO], [_ZERO])[0]
    at = (lambda key: key[0]) if rows else (lambda key: key[1])
    pos = (lambda key: key[1]) if rows else (lambda key: key[0])
    if zero != 0:
        # Multiplying implicit zeros by infinity really does produce NaN.
        # Do not expand an unbounded sparse row merely to discover this.
        Matrix._structural_limit(width + len(matrix._data))
        positions = range(width)
    else:
        positions = sorted({pos(key) for key in matrix._data if at(key) == target or (kind == 2 and at(key) == source)})
        Matrix._structural_limit(len(positions) + len(matrix._data))
    key = (lambda index, p: (index, p)) if rows else (lambda index, p: (p, index))
    targets = [matrix._data.get(key(target, p), _ZERO) for p in positions]
    sources = targets if kind == 0 else [matrix._data.get(key(source, p), _ZERO) for p in positions]
    values = Matrix._structural_values(kind == 2, coefficient, targets, sources)
    data = dict(matrix._data)
    for p, value in zip(positions, values):
        coordinate = key(target, p)
        if value != 0:
            data[coordinate] = value
        else:
            data.pop(coordinate, None)
    return matrix._new(matrix.rows, matrix.cols, data)


def callback(matrix, rows, index, source, function):
    from ._sparse import _index
    matrix._mutable()
    size, width = (matrix.rows, matrix.cols) if rows else (matrix.cols, matrix.rows)
    if not width:
        return
    index = _index(index, size)
    source = None if source is None else _index(source, size)
    Matrix._structural_limit(width + len(matrix._data))
    key = (lambda i, p: (i, p)) if rows else (lambda i, p: (p, i))
    for p in range(width):
        a = matrix[key(index, p)]
        b = p if source is None else matrix[key(source, p)]
        matrix[key(index, p)] = function(a, b)
