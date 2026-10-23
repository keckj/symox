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


def _is_Identity(M):
    # RepMatrix._eval_is_Identity
    if not M.is_square:
        return False
    if not all(M[i, i] == 1 for i in range(M.rows)):
        return False
    return sum(1 for v in M if not (v.is_Number and v.is_zero)) == M.rows


def _is_upper_hessenberg(M):
    return all(M[i, j].is_zero for i in range(M.rows) for j in range(M.cols)
               if i > j + 1)


def _is_lower_hessenberg(M):
    return all(M[i, j].is_zero for i in range(M.rows) for j in range(M.cols)
               if i + 1 < j)


def _hermitian_check(native):
    """``MatrixBase.is_hermitian``: ``fuzzy_and`` of
    ``simplify(M[i, j] - M[j, i].adjoint()).is_zero`` over the nonzero
    entries.  The native test decides most cases; it is undecided when its
    zero test cannot conclude, then SymPy's rule is applied."""
    def is_hermitian(M):
        r = native.__get__(M)
        if r is not None or not M.is_square:
            return r
        from . import simplify, adjoint
        return _fuzzy_and(simplify(M[i, j] - adjoint(M[j, i])).is_zero
                          for i in range(M.rows) for j in range(M.cols)
                          if M[i, j] != 0)
    is_hermitian.__doc__ = native.__doc__
    return property(is_hermitian)


def _dirac_conjugate(M):
    """``MatrixBase.D``: the Dirac conjugate ``M.H * mgamma(0)``."""
    from .physics.matrices import mgamma
    if M.rows != 4:
        raise AttributeError
    return M.H * mgamma(0)


def _minor_index_check(native, square):
    """``_minor_submatrix``/``_minor``/``_cofactor`` argument checks of
    ``sympy/matrices/determinant.py``."""
    def wrapper(M, i, j, *args, **kwargs):
        if square and (not M.is_square or M.rows < 1):
            from . import NonSquareMatrixError
            raise NonSquareMatrixError()
        ii = i + M.rows if i < 0 else i
        jj = j + M.cols if j < 0 else j
        if not 0 <= ii < M.rows or not 0 <= jj < M.cols:
            raise ValueError("`i` and `j` must satisfy 0 <= i < ``M.rows`` "
                             "(%d)" "and 0 <= j < ``M.cols`` (%d)."
                             % (M.rows, M.cols))
        return native(M, i, j, *args, **kwargs)
    wrapper.__doc__ = native.__doc__
    wrapper.__name__ = native.__name__
    return wrapper


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
    for name, fn in (('is_Identity', _is_Identity),
                     ('is_upper_hessenberg', _is_upper_hessenberg),
                     ('is_lower_hessenberg', _is_lower_hessenberg),
                     ('D', _dirac_conjugate)):
        if not hasattr(_core.Matrix, name):
            setattr(_core.Matrix, name, property(fn))
    _core.Matrix.is_hermitian = _hermitian_check(_core.Matrix.is_hermitian)
    for name in ('minor_submatrix', 'minor', 'cofactor'):
        setattr(_core.Matrix, name,
                _minor_index_check(getattr(_core.Matrix, name), name != 'minor_submatrix'))
    _core.Matrix.is_anti_symmetric = is_anti_symmetric
    _core.Matrix.is_diagonal = is_diagonal
