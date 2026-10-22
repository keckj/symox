"""SymPy's ``MatrixBase.extract`` index normalization and the structural
predicates ``is_symmetric``/``is_anti_symmetric`` with their ``simplify``
argument."""
from inspect import isfunction


def _is_sequence(value):
    return hasattr(value, '__getitem__') and hasattr(value, '__iter__') and not isinstance(value, str)


def _a2idx(j, n):
    # sympy.matrices.utilities.a2idx
    if not isinstance(j, int):
        try:
            j = j.__index__()
        except AttributeError:
            raise IndexError('Invalid index a[%r]' % (j,))
    if n is not None:
        if j < 0:
            j += n
        if not (j >= 0 and j < n):
            raise IndexError('Index out of range: a[%s]' % (j,))
    return int(j)


def extract_indices(matrix, rowsList, colsList):
    """The row and column lists of ``extract``: boolean masks select the
    positions of their True entries and negative indices count from the
    end."""
    if not _is_sequence(rowsList) or not _is_sequence(colsList):
        raise TypeError('rowsList and colsList must be iterable')
    rowsList, colsList = list(rowsList), list(colsList)
    if rowsList and all(isinstance(i, bool) for i in rowsList):
        rowsList = [index for index, item in enumerate(rowsList) if item]
    if colsList and all(isinstance(i, bool) for i in colsList):
        colsList = [index for index, item in enumerate(colsList) if item]
    rows, cols = matrix.shape
    return [_a2idx(k, rows) for k in rowsList], [_a2idx(k, cols) for k in colsList]


def _simpfunc(simplify):
    if isfunction(simplify):
        return simplify
    if simplify:
        from . import simplify as _simplify
        return _simplify
    return lambda value: value


def _fuzzy_and(values):
    result = True
    for value in values:
        if value is None:
            result = None
        elif not value:
            return False
    return result


def is_symmetric(matrix, simplify=True):
    """SymPy's ``RepMatrix.is_symmetric``: the square matrices whose
    ``self - self.T`` has no nonzero entry once ``simplify`` (``simplify``
    itself by default, a given function, or none) is applied."""
    simpfunc = _simpfunc(simplify)
    rows, cols = matrix.shape
    if rows != cols:
        return False
    for i in range(rows):
        for j in range(i + 1, cols):
            difference = matrix[i, j] - matrix[j, i]
            if difference != 0 and simpfunc(difference) != 0:
                return False
    return True


def is_anti_symmetric(matrix, simplify=True):
    """SymPy's ``MatrixBase.is_anti_symmetric``: the fuzzy conjunction of
    ``simpfunc(a[i, j] + a[j, i]).is_zero`` over the nonzero entries."""
    simpfunc = _simpfunc(simplify)
    rows, cols = matrix.shape
    if rows != cols:
        return False
    return _fuzzy_and(simpfunc(matrix[i, j] + matrix[j, i]).is_zero
                      for i in range(rows) for j in range(cols) if matrix[i, j] != 0)


def is_diagonal(matrix):
    """SymPy's ``MatrixBase.is_diagonal``: the fuzzy conjunction of
    ``is_zero`` over the nonzero off-diagonal entries."""
    rows, cols = matrix.shape
    return _fuzzy_and(matrix[i, j].is_zero for i in range(rows) for j in range(cols)
                      if i != j and matrix[i, j] != 0)


def install():
    from . import _core
    native_extract = _core.Matrix.extract

    def extract(self, rowsList, colsList):
        """Return a submatrix by specifying a list of rows and columns
        (SymPy's ``MatrixBase.extract``)."""
        rowsList, colsList = extract_indices(self, rowsList, colsList)
        return native_extract(self, rowsList, colsList)

    extract.__doc__ = native_extract.__doc__ or extract.__doc__
    _core.Matrix.extract = extract
    _core.Matrix.is_symmetric = is_symmetric
    _core.Matrix.is_anti_symmetric = is_anti_symmetric
    _core.Matrix.is_diagonal = is_diagonal
