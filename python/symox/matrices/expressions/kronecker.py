from ..._matrixexpr import KroneckerProduct, kronecker_product

__all__ = ['KroneckerProduct', 'kronecker_product']


def combine_kronecker(expr):
    """Combine KroneckerProduct with expression (SymPy's rules, native).

    Products ``KroneckerProduct(A, B)*KroneckerProduct(C, D)`` of matching
    shapes become ``KroneckerProduct(A*C, B*D)`` and integer powers of a
    Kronecker product of square matrices are distributed. SymPy also rewrites
    ``KroneckerProduct(A, B) + KroneckerProduct(C, D)`` as
    ``KroneckerProduct(A + C, B + D)``, which is not an identity, and
    distributes non-integer powers; both are left unchanged here.
    """
    from symox import _core
    if isinstance(expr, _core._MatrixExpr):
        return _core._combine_kronecker(expr)
    doit = getattr(expr, 'doit', None)
    return doit() if doit is not None else expr
__all__.append('combine_kronecker')
