"""SymPy's ``MatrixBase.pow`` dispatch for symox matrices.

Integer exponents keep the native binary powering. Other exponents follow
SymPy 1.14 (``pow``, ``_matrix_pow_by_jordan_blocks``): diagonal matrices are
raised entrywise, then the Jordan form is used, and a singular matrix raised
to a power that may be a nonnegative integer stays an unevaluated ``MatPow``.
The Jordan-block helpers of ``exp``/``log`` are SymPy's
``_eval_matrix_exp_jblock``/``_eval_matrix_log_jblock``.
"""

_native_pow = None


def _errors():
    from . import MatrixError, NonInvertibleMatrixError, NonSquareMatrixError
    return MatrixError, NonInvertibleMatrixError, NonSquareMatrixError


def _new_like(self, value):
    if getattr(self, '_immutable', False) and hasattr(value, 'as_immutable'):
        return value.as_immutable()
    if not getattr(self, '_immutable', False) and hasattr(value, 'as_mutable'):
        return value.as_mutable()
    return value


def _matrix_pow_by_jordan_blocks(self, num):
    from . import Matrix, binomial, diag, KroneckerDelta, sympify
    _, NonInvertibleMatrixError, _ = _errors()
    num = sympify(num)

    def jordan_cell_power(jc, n):
        N = jc.shape[0]
        l = jc[0, 0]
        if l.is_zero:
            if N == 1 and n.is_nonnegative:
                jc[0, 0] = l**n
            elif not (n.is_integer and n.is_nonnegative):
                raise NonInvertibleMatrixError(
                    "Non-invertible matrix can only be raised to a nonnegative integer")
            else:
                for i in range(N):
                    jc[0, i] = KroneckerDelta(i, n)
        else:
            for i in range(N):
                bn = binomial(n, i)
                if isinstance(bn, binomial):
                    bn = bn._eval_expand_func()
                jc[0, i] = l**(n - i)*bn
        for i in range(N):
            for j in range(1, N - i):
                jc[j, i + j] = jc[j - 1, i + j - 1]

    P, J = self.jordan_form()
    jordan_cells = [Matrix(j).as_mutable() for j in J.get_diag_blocks()]
    for j in jordan_cells:
        jordan_cell_power(j, num)
    return _new_like(self, P.multiply(diag(*jordan_cells)).multiply(P.inv()))


def _eval_pow_by_recursion(self, num):
    if num == 1:
        return self
    if num % 2 == 1:
        a, b = self, _eval_pow_by_recursion(self, num - 1)
    else:
        a = b = _eval_pow_by_recursion(self, num // 2)
    return a.multiply(b)


def _is_integral_number(exp):
    return bool(exp.is_Number and exp % 1 == 0)


def pow(self, exp, method=None):
    """Return ``self**exp`` (SymPy's ``MatrixBase.pow``)."""
    from . import Integer, sympify
    _, NonInvertibleMatrixError, NonSquareMatrixError = _errors()
    if method is not None and method not in ['multiply', 'mulsimp', 'jordan', 'cayley']:
        raise TypeError('No such method')
    if self.rows != self.cols:
        raise NonSquareMatrixError()
    a = self
    exp = sympify(exp)
    integral = _is_integral_number(exp)
    if method is None and integral:
        # Native binary powering (SymPy's recursion/cayley heuristics give
        # the same exact values).
        return _native_pow(a, Integer(exp) if exp.is_Float else exp)
    if exp.is_zero:
        return _new_like(a, a.eye(a.rows))
    if exp == 1:
        return a
    if a.is_diagonal():
        out = a.as_mutable().copy()
        for i in range(a.rows):
            out[i, i] = a[i, i]**exp
        return _new_like(a, out)
    if integral:
        if a.rows == 1:
            return _new_like(a, a.__class__([[a[0]**exp]]))
        if exp < 0:
            exp = -exp
            a = a.inv()
    if method == 'jordan':
        return _matrix_pow_by_jordan_blocks(a, exp)
    if method in ('cayley', 'mulsimp', 'multiply'):
        if not integral:
            raise ValueError("%s method is only valid for integer powers" % method)
        if method == 'multiply':
            return _new_like(self, _eval_pow_by_recursion(a, int(exp)))
        return _new_like(self, _native_pow(a, Integer(exp)))
    try:
        return _matrix_pow_by_jordan_blocks(a, exp)
    except NonInvertibleMatrixError:
        if exp.is_integer is False or exp.is_nonnegative is False:
            raise
    from .matrices.expressions import MatPow
    return MatPow(a, exp)


def __pow__(self, exp, modulo=None):
    if modulo is not None:
        return _native_pow(self, exp, modulo)
    if isinstance(exp, int) and not isinstance(exp, bool):
        return _native_pow(self, exp)
    return pow(self, exp)


def _banded_like(self, size, bands):
    from . import Matrix
    from .matrices.sparsetools import banded
    return _new_like(self, Matrix(banded(size, bands)))


def _eval_matrix_exp_jblock(self):
    from . import exp, factorial
    size = self.rows
    exp_l = exp(self[0, 0])
    return _banded_like(self, size, {i: exp_l / factorial(i) for i in range(size)})


def _eval_matrix_log_jblock(self):
    from . import log
    MatrixError, _, _ = _errors()
    size = self.rows
    l = self[0, 0]
    if l.is_zero:
        raise MatrixError(
            'Could not take logarithm or reciprocal for the given '
            'eigenvalue {}'.format(l))
    bands = {0: log(l)}
    for i in range(1, size):
        bands[i] = -((-l) ** -i) / i
    return _banded_like(self, size, bands)


def install(Matrix):
    global _native_pow
    if _native_pow is None:
        _native_pow = Matrix.__pow__
    Matrix.__pow__ = __pow__
    Matrix.pow = pow
    Matrix._matrix_pow_by_jordan_blocks = _matrix_pow_by_jordan_blocks
    Matrix._eval_pow_by_recursion = _eval_pow_by_recursion
    Matrix._eval_matrix_exp_jblock = _eval_matrix_exp_jblock
    Matrix._eval_matrix_log_jblock = _eval_matrix_log_jblock
