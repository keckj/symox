from ..._matrixexpr import KroneckerProduct, kronecker_product
from ..._kind import matrix_tag

__all__ = ['KroneckerProduct', 'kronecker_product']


def _kronecker_mat_mul(expr):
    from symox import MatMul
    factor, matrices = expr.as_coeff_matrices()
    matrices = list(matrices)
    i = 0
    while i < len(matrices) - 1:
        a, b = matrices[i:i+2]
        if isinstance(a, KroneckerProduct) and isinstance(b, KroneckerProduct) and a.has_matching_shape(b):
            matrices[i] = KroneckerProduct(*[x*y for x, y in zip(a.args, b.args)])
            matrices.pop(i+1)
        else:
            i += 1
    return factor*MatMul(*matrices)


def _kronecker_mat_pow(expr):
    from symox import MatPow
    base, exp = expr.args
    # (A x B)**k = A**k x B**k by the mixed-product property, for integer k.
    if isinstance(base, KroneckerProduct) and exp.is_integer and all(a.is_square for a in base.args):
        return KroneckerProduct(*[MatPow(a, exp) for a in base.args])
    return expr


def _combine_step(expr):
    from symox import _core
    if not isinstance(expr, _core._MatrixExpr):
        return expr
    args = expr.args
    new = [_combine_step(a) for a in args]
    if any(a is not b and a != b for a, b in zip(new, args)):
        expr = expr.func(*new)
        if not isinstance(expr, _core._MatrixExpr):
            return expr
    kind = matrix_tag(expr)
    if kind == 'MatMul':
        return _kronecker_mat_mul(expr)
    if kind == 'MatPow':
        return _kronecker_mat_pow(expr)
    return expr


def combine_kronecker(expr):
    """Combine KroneckerProduct with expression (SymPy's rules).

    Products ``KroneckerProduct(A, B)*KroneckerProduct(C, D)`` of matching
    shapes become ``KroneckerProduct(A*C, B*D)`` and integer powers of a
    Kronecker product of square matrices are distributed. SymPy also rewrites
    ``KroneckerProduct(A, B) + KroneckerProduct(C, D)`` as
    ``KroneckerProduct(A + C, B + D)``, which is not an identity; sums are
    left unchanged here.
    """
    while True:
        new = _combine_step(expr)
        if new == expr:
            break
        expr = new
    doit = getattr(expr, 'doit', None)
    return doit() if doit is not None else expr
__all__.append('combine_kronecker')
