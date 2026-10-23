"""Ports of SymPy 1.14's matrix definiteness predicates
and the mpmath eigenvalue path for numeric Float matrices together with the
``simplify``/``rational`` flags of ``eigenvals``
(``sympy/matrices/eigen.py``), ``is_echelon`` (``reductions.py``) and the
pivot search ``_find_reasonable_pivot`` (``determinant.py``)."""


def _load():
    # Deferred: this module is installed while ``symox`` is still importing.
    g = globals()
    if 'sqrt' not in g:
        from .core.logic import fuzzy_and, fuzzy_or
        from . import Float, Integer, S, sqrt
        g.update(fuzzy_and=fuzzy_and, fuzzy_or=fuzzy_or, Float=Float,
                 Integer=Integer, S=S, sqrt=sqrt)


def _iszero(x):
    return getattr(x, 'is_zero', None)


def _simplify(expr):
    from . import simplify
    return simplify(expr)


def _find_reasonable_pivot(col, iszerofunc=_iszero, simpfunc=_simplify):
    """ Find the lowest index of an item in ``col`` that is
    suitable for a pivot.  If ``col`` consists only of
    Floats, the pivot with the largest norm is returned.
    Otherwise, the first element where ``iszerofunc`` returns
    False is used.  If ``iszerofunc`` does not return false,
    items are simplified and retested until a suitable
    pivot is found.

    Returns a 4-tuple
        (pivot_offset, pivot_val, assumed_nonzero, newly_determined)
    where pivot_offset is the index of the pivot, pivot_val is
    the (possibly simplified) value of the pivot, assumed_nonzero
    is True if an assumption that the pivot was non-zero
    was made without being proved, and newly_determined are
    elements that were simplified during the process of pivot
    finding."""

    newly_determined = []
    col = list(col)
    # a column that contains a mix of floats and integers
    # but at least one float is considered a numerical
    # column, and so we do partial pivoting
    if all(isinstance(x, (Float, Integer)) for x in col) and any(
            isinstance(x, Float) for x in col):
        col_abs = [abs(x) for x in col]
        max_value = max(col_abs)
        if iszerofunc(max_value):
            # just because iszerofunc returned True, doesn't
            # mean the value is numerically zero.  Make sure
            # to replace all entries with numerical zeros
            if max_value != 0:
                newly_determined = [(i, 0) for i, x in enumerate(col) if x != 0]
            return (None, None, False, newly_determined)
        index = col_abs.index(max_value)
        return (index, col[index], False, newly_determined)

    # PASS 1 (iszerofunc directly)
    possible_zeros = []
    for i, x in enumerate(col):
        is_zero = iszerofunc(x)
        # is someone wrote a custom iszerofunc, it may return
        # BooleanFalse or BooleanTrue instead of True or False,
        # so use == for comparison instead of `is`
        if is_zero == False:
            # we found something that is definitely not zero
            return (i, x, False, newly_determined)
        possible_zeros.append(is_zero)

    # by this point, we've found no certain non-zeros
    if all(possible_zeros):
        # if everything is definitely zero, we have
        # no pivot
        return (None, None, False, newly_determined)

    # PASS 2 (iszerofunc after simplify)
    # we haven't found any for-sure non-zeros, so
    # go through the elements iszerofunc couldn't
    # make a determination about and opportunistically
    # simplify to see if we find something
    for i, x in enumerate(col):
        if possible_zeros[i] is not None:
            continue
        simped = simpfunc(x)
        is_zero = iszerofunc(simped)
        if is_zero in (True, False):
            newly_determined.append((i, simped))
        if is_zero == False:
            return (i, simped, False, newly_determined)
        possible_zeros[i] = is_zero

    # after simplifying, some things that were recognized
    # as zeros might be zeros
    if all(possible_zeros):
        # if everything is definitely zero, we have
        # no pivot
        return (None, None, False, newly_determined)

    # PASS 3 (.equals(0))
    # some expressions fail to simplify to zero, but
    # ``.equals(0)`` evaluates to True.  As a last-ditch
    # attempt, apply ``.equals`` to these expressions
    for i, x in enumerate(col):
        if possible_zeros[i] is not None:
            continue
        if x.equals(S.Zero):
            # ``.iszero`` may return False with
            # an implicit assumption (e.g., ``x.equals(0)``
            # when ``x`` is a symbol), so only treat it
            # as proved when ``.equals(0)`` returns True
            possible_zeros[i] = True
            newly_determined.append((i, S.Zero))

    if all(possible_zeros):
        return (None, None, False, newly_determined)

    # at this point there is nothing that could definitely
    # be a pivot.  To maintain compatibility with existing
    # behavior, we'll assume that an illdetermined thing is
    # non-zero.  We should probably raise a warning in this case
    i = possible_zeros.index(None)
    return (i, col[i], True, newly_determined)


def _has_positive_diagonals(M):
    return fuzzy_and(M[i, i].is_positive for i in range(M.rows))


def _has_nonnegative_diagonals(M):
    return fuzzy_and(M[i, i].is_nonnegative for i in range(M.rows))


def _diagonally_dominant(M, strict):
    # MatrixBase.is_weakly/strongly_diagonally_dominant
    from . import Abs
    if not M.is_square:
        return False
    rows, cols = M.shape

    def test_row(i):
        summation = S.Zero
        for j in range(cols):
            if i != j:
                summation += Abs(M[i, j])
        diff = Abs(M[i, i]) - summation
        return diff.is_positive if strict else diff.is_nonnegative
    return fuzzy_and(test_row(i) for i in range(rows))


def _is_weakly_diagonally_dominant(M):
    return _diagonally_dominant(M, False)


def _is_strongly_diagonally_dominant(M):
    return _diagonally_dominant(M, True)


def _fuzzy_positive_definite(M):
    positive_diagonals = _has_positive_diagonals(M)
    if positive_diagonals is False:
        return False
    if positive_diagonals and _is_strongly_diagonally_dominant(M):
        return True
    return None


def _fuzzy_positive_semidefinite(M):
    nonnegative_diagonals = _has_nonnegative_diagonals(M)
    if nonnegative_diagonals is False:
        return False
    if nonnegative_diagonals and _is_weakly_diagonally_dominant(M):
        return True
    return None


def _is_positive_definite(M):
    if not M.is_hermitian:
        if not M.is_square:
            return False
        M = M + M.H

    fuzzy = _fuzzy_positive_definite(M)
    if fuzzy is not None:
        return fuzzy

    return _is_positive_definite_GE(M)


def _is_positive_semidefinite(M):
    if not M.is_hermitian:
        if not M.is_square:
            return False
        M = M + M.H

    fuzzy = _fuzzy_positive_semidefinite(M)
    if fuzzy is not None:
        return fuzzy

    return _is_positive_semidefinite_cholesky(M)


def _is_negative_definite(M):
    return _is_positive_definite(-M)


def _is_negative_semidefinite(M):
    return _is_positive_semidefinite(-M)


def _is_indefinite(M):
    if M.is_hermitian:
        eigen = M.eigenvals()
        args1        = [x.is_positive for x in eigen.keys()]
        any_positive = fuzzy_or(args1)
        args2        = [x.is_negative for x in eigen.keys()]
        any_negative = fuzzy_or(args2)

        return fuzzy_and([any_positive, any_negative])

    elif M.is_square:
        return (M + M.H).is_indefinite

    return False


def _is_positive_definite_GE(M):
    """A division-free gaussian elimination method for testing
    positive-definiteness."""
    M = M.as_mutable()
    size = M.rows

    for i in range(size):
        is_positive = M[i, i].is_positive
        if is_positive is not True:
            return is_positive
        for j in range(i+1, size):
            M[j, i+1:] = M[i, i] * M[j, i+1:] - M[j, i] * M[i, i+1:]
    return True


def _is_positive_semidefinite_cholesky(M):
    """Uses Cholesky factorization with complete pivoting

    References
    ==========

    .. [1] http://eprints.ma.man.ac.uk/1199/1/covered/MIMS_ep2008_116.pdf

    .. [2] https://www.value-at-risk.net/cholesky-factorization/
    """
    M = M.as_mutable()
    for k in range(M.rows):
        diags = [M[i, i] for i in range(k, M.rows)]
        pivot, pivot_val, nonzero, _ = _find_reasonable_pivot(diags)

        if nonzero:
            return None

        if pivot is None:
            for i in range(k+1, M.rows):
                for j in range(k, M.cols):
                    iszero = M[i, j].is_zero
                    if iszero is None:
                        return None
                    elif iszero is False:
                        return False
            return True

        if M[k, k].is_negative or pivot_val.is_negative:
            return False
        elif not (M[k, k].is_nonnegative and pivot_val.is_nonnegative):
            return None

        if pivot > 0:
            M.col_swap(k, k+pivot)
            M.row_swap(k, k+pivot)

        M[k, k] = sqrt(M[k, k])
        M[k, k+1:] /= M[k, k]
        M[k+1:, k+1:] -= M[k, k+1:].H * M[k, k+1:]

    return M[-1, -1].is_nonnegative


def _is_echelon(M, iszerofunc=_iszero):
    """Returns `True` if the matrix is in echelon form. That is, all rows of
    zeros are at the bottom, and below each leading non-zero in a row are
    exclusively zeros."""

    if M.rows <= 0 or M.cols <= 0:
        return True

    zeros_below = all(iszerofunc(t) for t in M[1:, 0])

    if iszerofunc(M[0, 0]):
        return zeros_below and _is_echelon(M[:, 1:], iszerofunc)

    return zeros_below and _is_echelon(M[1:, 1:], iszerofunc)


DEFAULT_MAXPREC = 333


def _float_numeric(M):
    from . import Float
    return M.has(Float) and all(x.is_number for x in M)


def _eigenvals_eigenvects_mpmath(M):
    # sympy.matrices.eigen._eigenvals_eigenvects_mpmath
    import mpmath as mp
    from mpmath import workprec
    from mpmath.libmp.libmpf import prec_to_dps
    from . import Float, PrecisionExhausted
    norm2 = lambda v: mp.sqrt(sum(i**2 for i in v))

    def entry(x):
        re, im = x.as_real_imag()
        return mp.mpc(mp.mpmathify(str(re)), mp.mpmathify(str(im)))

    v1 = None
    prec = max(x._prec for x in M.atoms(Float))
    eps = 2**-prec
    while prec < DEFAULT_MAXPREC:
        with workprec(prec):
            N = M.evalf(n=prec_to_dps(prec))
            A = mp.matrix([[entry(N[i, j]) for j in range(N.cols)]
                           for i in range(N.rows)])
            E, ER = mp.eig(A)
            v2 = norm2([i for e in E for i in (mp.re(e), mp.im(e))])
            if v1 is not None and mp.fabs(v1 - v2) < eps:
                return E, ER
            v1 = v2
        prec *= 2
    raise PrecisionExhausted


def _mp_to_expr(x):
    import mpmath as mp
    from . import sympify
    if isinstance(x, mp.mpc) and x.imag == 0:
        x = x.real
    return sympify(x)


def _eigenvals_mpmath(M, multiple=False):
    from collections import Counter
    E, _ = _eigenvals_eigenvects_mpmath(M)
    result = [_mp_to_expr(x) for x in E]
    if multiple:
        return result
    return dict(Counter(result))


def _eigenvects_mpmath(M):
    from . import Matrix
    E, ER = _eigenvals_eigenvects_mpmath(M)
    result = []
    for i in range(M.rows):
        eigenval = _mp_to_expr(E[i])
        eigenvect = Matrix([[_mp_to_expr(ER[j, i])] for j in range(ER.rows)])
        result.append((eigenval, 1, [eigenvect]))
    return result


_native_eigenvals = None
_wrapped_eigenvects = None


def _eigenvals(M, error_when_incomplete=True, *, simplify=False,
               multiple=False, rational=False, **flags):
    """SymPy's ``_eigenvals`` dispatch: mpmath for numeric Float matrices,
    the ``rational`` conversion and the ``simplify`` post-processing around
    the native exact kernel."""
    from types import FunctionType
    from . import Float, NonSquareMatrixError, nsimplify
    if not M:
        return [] if multiple else {}
    if not M.is_square:
        raise NonSquareMatrixError("{} must be a square matrix.".format(M))
    if _float_numeric(M):
        return _eigenvals_mpmath(M, multiple=multiple)
    if rational:
        M = M.applyfunc(
            lambda x: nsimplify(x, rational=True) if x.has(Float) else x)
    try:
        ret = _native_eigenvals(M, error_when_incomplete, multiple, **flags)
    except NotImplementedError as exc:
        from . import _matrix_eigen_rootof as rootof
        if not rootof.is_unsupported(exc):
            raise
        ret = rootof.eigenvals(M, error_when_incomplete, multiple,
                               simplify, **flags)
    if not simplify:
        return ret
    if not isinstance(simplify, FunctionType):
        simplify = _simplify
    if multiple:
        return [simplify(value) for value in ret]
    return {simplify(key): value for key, value in ret.items()}


def _eigenvects(M, error_when_incomplete=True, iszerofunc=None, **flags):
    if _float_numeric(M):
        return _eigenvects_mpmath(M)
    try:
        return _wrapped_eigenvects(M, error_when_incomplete, iszerofunc, **flags)
    except NotImplementedError as exc:
        from . import _matrix_eigen_rootof as rootof
        from . import CRootOf, MatrixError
        if not rootof.is_unsupported(exc):
            raise
        if rootof._is_rational_matrix(M):
            return rootof.eigenvects(M, flags.get('simplify', True),
                                     flags.get('simplify', False))
        # _eigenvects_sympy
        eigs = rootof.eigenvals(M, error_when_incomplete, False)
        if any(x.has(CRootOf) for x in eigs):
            raise MatrixError(
                "Eigenvector computation is not implemented if the matrix have "
                "eigenvalues in CRootOf form")
        raise


def _jordan_form(M, calc_transform=True, *, chop=False, **kwargs):
    from . import Float
    has_floats = M.is_square and M.has(Float)
    if has_floats:
        # SymPy runs the whole computation on nsimplify(x, rational=True)
        # and evaluates the results at the largest precision of the entries
        from . import nsimplify
        from .printing._mpf_str import prec_to_dps
        precs = [x._prec for x in M if isinstance(x, Float)]
        max_dps = max(prec_to_dps(max(precs) if precs else 53), 15)
        R = M.applyfunc(lambda x: nsimplify(x, rational=True))
        res = _jordan_form(R, calc_transform, **kwargs)
        if calc_transform:
            # restore_floats returns a list here, as in SymPy
            return [m.evalf(n=max_dps, chop=chop) for m in res]
        return res.evalf(n=max_dps, chop=chop)
    try:
        return _native_jordan_form(M, calc_transform, **kwargs)
    except NotImplementedError as exc:
        from . import _matrix_eigen_rootof as rootof
        from . import MatrixError
        if not rootof.is_unsupported(exc):
            raise
        raise MatrixError(
            "Jordan normal form is not implemented if the matrix have "
            "eigenvalues in CRootOf form")


_native_jordan_form = None


def install():
    from . import _core
    from ._sparse import _SparseBase
    def wrap(fn):
        def run(M):
            _load()
            return fn(M)
        return run

    global _native_eigenvals, _wrapped_eigenvects
    if _native_eigenvals is None:
        _native_eigenvals = _core.Matrix.eigenvals
        _wrapped_eigenvects = _core.Matrix.eigenvects
        _core.Matrix.eigenvals = _eigenvals
        _core.Matrix.eigenvects = _eigenvects
    global _native_jordan_form
    if _native_jordan_form is None:
        _native_jordan_form = _core.Matrix.jordan_form
        _core.Matrix.jordan_form = _jordan_form
    for name, fn in (('is_positive_definite', _is_positive_definite),
                     ('is_positive_semidefinite', _is_positive_semidefinite),
                     ('is_negative_definite', _is_negative_definite),
                     ('is_negative_semidefinite', _is_negative_semidefinite),
                     ('is_indefinite', _is_indefinite),
                     ('is_echelon', _is_echelon),
                     ('is_weakly_diagonally_dominant', _is_weakly_diagonally_dominant),
                     ('is_strongly_diagonally_dominant', _is_strongly_diagonally_dominant)):
        if not hasattr(_core.Matrix, name):
            setattr(_core.Matrix, name, property(wrap(fn)))
        if not hasattr(_SparseBase, name):
            setattr(_SparseBase, name,
                    property(lambda self, _fn=wrap(fn): _fn(self.as_dense())))
