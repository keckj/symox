"""SymPy 1.14's invertibility check of ``inverse_LU`` and the ``BLOCK``
inversion method (``sympy/matrices/inverse.py``) around the native kernels."""


def _iszero(x):
    return getattr(x, 'is_zero', None)


def _verify_invertible(M, iszerofunc=_iszero):
    from . import NonSquareMatrixError, NonInvertibleMatrixError
    if not M.is_square:
        raise NonSquareMatrixError("A Matrix must be square to invert.")
    d = M.det(method='berkowitz')
    zero = d.equals(0)
    if zero is None:
        ok = M.rref(simplify=True)[0]
        zero = any(iszerofunc(ok[j, j]) for j in range(ok.rows))
    if zero:
        raise NonInvertibleMatrixError("Matrix det == 0; not invertible.")
    return d


def _has_nc_entry(M):
    from . import _core
    return any(isinstance(x, _core.Expr) and _core._has_noncommutative_atoms(x)
               and not x.is_commutative for x in M)


def _like(M, rv):
    # The result keeps the class and the immutability of ``M``.
    from . import ImmutableMatrix
    if isinstance(M, ImmutableMatrix) and not isinstance(rv, ImmutableMatrix):
        rv = rv.as_immutable()
    if type(rv) is not type(M):
        try:
            rv = type(M)(rv)
        except Exception:
            pass
    return rv


_native_inverse_LU = None
_native_inv = None


def inverse_LU(M, iszerofunc=_iszero):
    from . import NonSquareMatrixError
    if not M.is_square:
        raise NonSquareMatrixError("A Matrix must be square to invert.")
    if M.free_symbols:
        _verify_invertible(M, iszerofunc=iszerofunc)
    if iszerofunc is _iszero:
        return _native_inverse_LU(M)
    # SymPy's LUsolve ignores a custom predicate; the native kernel applies
    # it to the pivots (stricter: it may raise, never a wrong value).
    return _native_inverse_LU(M, iszerofunc=iszerofunc)


def inverse_BLOCK(M, iszerofunc=_iszero):
    from . import NonInvertibleMatrixError
    from .matrices.expressions.blockmatrix import BlockMatrix
    i = M.shape[0]
    if i <= 20:
        return M.inv(method="LU", iszerofunc=_iszero)
    A = M[:i // 2, :i // 2]
    B = M[:i // 2, i // 2:]
    C = M[i // 2:, :i // 2]
    D = M[i // 2:, i // 2:]
    try:
        D_inv = inverse_BLOCK(D)
    except NonInvertibleMatrixError:
        return M.inv(method="LU", iszerofunc=_iszero)
    B_D_i = B*D_inv
    BDC = B_D_i*C
    A_n = A - BDC
    try:
        A_n = inverse_BLOCK(A_n)
    except NonInvertibleMatrixError:
        return M.inv(method="LU", iszerofunc=_iszero)
    B_n = -A_n*B_D_i
    dc = D_inv*C
    C_n = -dc*A_n
    D_n = D_inv + dc*-B_n
    return BlockMatrix([[A_n, B_n], [C_n, D_n]]).as_explicit()


def _inexact_domain(M):
    """``construct_domain`` of the entries when it is ``RR`` or ``CC``:
    returns ``(is_complex, prec)``, or None for any other domain."""
    from .core.evalf import pure_complex
    floats = complexes = False
    precs = []
    for e in M:
        if e.is_Rational:
            continue
        if e.is_Float:
            floats = True
            precs.append(e._prec)
            continue
        parts = pure_complex(e)
        if not parts or not all(v.is_Rational or v.is_Float for v in parts):
            return None
        complexes = True
        for v in parts:
            if v.is_Float:
                floats = True
                precs.append(v._prec)
    if not floats:
        return None
    return complexes, max(precs) if precs else 53


def _inv_inexact(M, complexes, prec):
    """SymPy's ``_inv_DM`` over ``RR``/``CC``: ``inv_den`` is the dense
    Gauss-Jordan ``ddm_irref`` of ``[M | I]`` with partial pivoting, in
    floating-point arithmetic at the domain precision."""
    import mpmath
    from mpmath.libmp.libmpf import prec_to_dps
    from . import Float, I, Integer, NonInvertibleMatrixError, Matrix
    n = M.rows
    dps = prec_to_dps(prec)
    with mpmath.workprec(prec):
        def conv(e):
            v = e.evalf(n=dps)
            if complexes:
                re_, im_ = v.as_real_imag()
                return mpmath.mpc(mpmath.mpf(re_._mpf_ if re_.is_Float else
                                             Float(re_, dps)._mpf_),
                                  mpmath.mpf(im_._mpf_ if im_.is_Float else
                                             Float(im_, dps)._mpf_))
            return mpmath.mpf(v._mpf_ if v.is_Float else Float(v, dps)._mpf_)
        zero = mpmath.mpc(0) if complexes else mpmath.mpf(0)
        one = mpmath.mpc(1) if complexes else mpmath.mpf(1)
        a = [[conv(M[i, j]) for j in range(n)] +
             [one if k == i else zero for k in range(n)] for i in range(n)]
        m, w = n, 2 * n
        i = 0
        pivots = []
        for j in range(w):
            ip = max(range(i, m), key=lambda ip: abs(a[ip][j]))
            a[i], a[ip] = a[ip], a[i]
            aij = a[i][j]
            if not aij:
                for ip in range(i + 1, m):
                    aij = a[ip][j]
                    if aij:
                        a[i], a[ip] = a[ip], a[i]
                        break
                else:
                    continue
            ai = a[i]
            aijinv = aij**-1
            for l in range(j, w):
                ai[l] *= aijinv
            for k, ak in enumerate(a):
                if k == i or not ak[j]:
                    continue
                akj = ak[j]
                ak[j] -= akj
                for l in range(j + 1, w):
                    ak[l] -= akj * ai[l]
            pivots.append(j)
            i += 1
            if i >= m:
                break
        if len(pivots) != n or (pivots and pivots[-1] >= n):
            raise NonInvertibleMatrixError("Matrix det == 0; not invertible.")

        def back(v):
            if not v:
                return Integer(0)
            if complexes:
                return Float(v.real, dps) + I*Float(v.imag, dps)
            return Float(v, dps)
        return Matrix([[back(a[r][n + c]) for c in range(n)] for r in range(n)])


def inv(M, method=None, iszerofunc=None, try_block_diag=False, **kwargs):
    if method is None and iszerofunc is None and not try_block_diag \
            and not kwargs and M.is_square and M.rows:
        dom = _inexact_domain(M)
        if dom is not None:
            return _like(M, _inv_inexact(M, *dom))
    if method in ('LU', 'BLOCK') and not try_block_diag and not kwargs:
        from . import NonSquareMatrixError
        if not M.is_square:
            raise NonSquareMatrixError("A Matrix must be square to invert.")
        z = iszerofunc if iszerofunc is not None else _iszero
        rv = (M.inverse_LU(iszerofunc=z) if method == 'LU'
              else inverse_BLOCK(M, iszerofunc=z))
        return _like(M, rv)
    if method is None and iszerofunc is None and not try_block_diag \
            and M.is_square and _has_nc_entry(M):
        # SymPy's default path first converts to a DomainMatrix, which
        # rejects noncommutative entries.
        from . import PolynomialError
        raise PolynomialError("non-commutative expressions are not supported")
    if iszerofunc is not None:
        kwargs['iszerofunc'] = iszerofunc
    if try_block_diag:
        kwargs['try_block_diag'] = try_block_diag
    return _native_inv(M, method=method, **kwargs)


def install():
    from . import _core
    global _native_inverse_LU, _native_inv
    if _native_inverse_LU is not None:
        return
    _native_inverse_LU = _core.Matrix.inverse_LU
    _native_inv = _core.Matrix.inv
    inverse_LU.__doc__ = _native_inverse_LU.__doc__
    inv.__doc__ = _native_inv.__doc__
    _core.Matrix.inverse_LU = inverse_LU
    _core.Matrix.inverse_BLOCK = inverse_BLOCK
    _core.Matrix.inv = inv
