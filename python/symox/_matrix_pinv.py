"""SymPy's ``_pinv``/``_pinv_diagonalization`` (sympy/matrices/inverse.py).

The diagonalization method uses the normalized eigenvectors of ``A.H*A``
(or ``A*A.H``) and their conjugate transpose, as SymPy does, instead of
inverting and cancelling the eigenvector matrix (which blew up on symbolic
entries). The rank-decomposition method stays native.
"""

_native_pinv = None


def _pinv_diagonalization(M):
    from . import MatrixError
    if M.is_zero_matrix:
        return M.H
    A = M
    AH = M.H
    try:
        if M.rows >= M.cols:
            P, D = AH.multiply(A).diagonalize(normalize=True)
            D_pinv = D.applyfunc(lambda x: 0 if x.is_zero else 1 / x)
            return P.multiply(D_pinv).multiply(P.H).multiply(AH)
        P, D = A.multiply(AH).diagonalize(normalize=True)
        D_pinv = D.applyfunc(lambda x: 0 if x.is_zero else 1 / x)
        return AH.multiply(P).multiply(D_pinv).multiply(P.H)
    except MatrixError:
        raise NotImplementedError(
            'pinv for rank-deficient matrices is not implemented yet.')


def pinv(self, method='RD'):
    if method == 'ED':
        result = _pinv_diagonalization(self)
        if getattr(self, '_immutable', False) and hasattr(result, 'as_immutable'):
            return result.as_immutable()
        return result
    return _native_pinv(self, method)


def install(Matrix):
    global _native_pinv
    if _native_pinv is None:
        _native_pinv = Matrix.pinv
    Matrix.pinv = pinv
