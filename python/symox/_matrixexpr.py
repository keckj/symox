"""SymPy-style constructors for the native noncommutative matrix DAG."""
from . import _core


class _MatrixType(type):
    def __instancecheck__(cls, value):
        return isinstance(value, _core._MatrixExpr) and (cls._kinds is None or value.kind in cls._kinds)


class MatrixExpr(metaclass=_MatrixType):
    _kinds = None

    def __new__(cls, value):
        return _core._as_matrixexpr(value)


class MatrixSymbol(MatrixExpr):
    _kinds = ('MatrixSymbol',)

    def __new__(cls, name, n, m):
        return _core._matrixexpr_node('MatrixSymbol', (str(name), n, m))


class MatAdd(MatrixExpr):
    _kinds = ('MatAdd',)

    def __new__(cls, *args, evaluate=False, **options):
        if options:
            raise NotImplementedError('MatAdd options are not implemented: '+', '.join(options))
        return _core._matrixexpr_node('MatAdd', args, evaluate)


class MatMul(MatrixExpr):
    _kinds = ('MatMul',)

    def __new__(cls, *args, evaluate=False, **options):
        if options:
            raise NotImplementedError('MatMul options are not implemented: '+', '.join(options))
        return _core._matrixexpr_node('MatMul', args, evaluate)


class MatPow(MatrixExpr):
    _kinds = ('MatPow','Inverse')

    def __new__(cls, base, exp, evaluate=False):
        return _core._matrixexpr_node('MatPow',(base,exp),evaluate)


class Inverse(MatPow):
    _kinds = ('Inverse',)

    def __new__(cls, arg, exp=-1):
        if exp != -1:
            raise ValueError('Inverse exponent must be -1')
        return _core._matrixexpr_node('Inverse',(arg,))


class Transpose(MatrixExpr):
    _kinds = ('Transpose',)

    def __new__(cls, arg):
        return _core._matrixexpr_node('Transpose',(arg,))


class Adjoint(MatrixExpr):
    _kinds = ('Adjoint',)

    def __new__(cls, arg):
        return _core._matrixexpr_node('Adjoint',(arg,))


class MatrixConjugate(MatrixExpr):
    _kinds = ('MatrixConjugate',)

    def __new__(cls, arg):
        return _core._matrixexpr_node('MatrixConjugate',(arg,))


class Identity(MatrixExpr):
    _kinds = ('Identity',)

    def __new__(cls, n):
        return _core._matrixexpr_node('Identity',(n,))


class ZeroMatrix(MatrixExpr):
    _kinds = ('ZeroMatrix',)

    def __new__(cls, n, m):
        return _core._matrixexpr_node('ZeroMatrix',(n,m))


class _ElementType(type):
    def __instancecheck__(cls, value):
        return isinstance(value, _core.Expr) and getattr(value.func,'__name__',None)=='MatrixElement'


class MatrixElement(metaclass=_ElementType):
    def __new__(cls, parent, i, j):
        return _core._matrix_element(parent,i,j)


__all__ = ['MatrixExpr','MatrixSymbol','MatAdd','MatMul','MatPow','Transpose',
           'Adjoint','Inverse','Identity','ZeroMatrix','MatrixElement','MatrixConjugate']
