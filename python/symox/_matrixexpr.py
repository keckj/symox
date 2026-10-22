"""SymPy-style constructors for the native noncommutative matrix DAG."""
from . import _core

# The native string tag of a node (MatrixSymbol, MatMul, ...); the public
# ``kind`` becomes SymPy's MatrixKind once the kinds are installed.
if '_kind_tag' not in _core._MatrixExpr.__dict__:
    _core._MatrixExpr._kind_tag = _core._MatrixExpr.__dict__['kind']


class _MatrixType(type):
    def __getattr__(cls,name):
        from ._constructor_methods import instance_method
        return instance_method(_core._MatrixExpr,name)

    def __instancecheck__(cls, value):
        return isinstance(value, _core._MatrixExpr) and (cls._kinds is None or value._kind_tag in cls._kinds)

    def class_key(cls):
        name = cls.__name__
        if name in ('MatMul', 'MatAdd'):
            return (3, ('MatMul', 'MatAdd').index(name), name)
        return (5, 0, name)

    def fromiter(cls, args, **assumptions):
        return cls(*args, **assumptions)


class MatrixExpr(metaclass=_MatrixType):
    _kinds = None

    def __new__(cls, value):
        return _core._as_matrixexpr(value)

    @staticmethod
    def from_index_summation(expr, first_index=None, last_index=None, dimensions=None):
        from ._matrix_index_sums import from_index_summation
        return from_index_summation(expr, first_index, last_index, dimensions)


class MatrixSymbol(MatrixExpr):
    _kinds = ('MatrixSymbol',)

    def __new__(cls, name, n, m):
        return _core._matrixexpr_node('MatrixSymbol', (str(name), n, m))


def _normalize_matrix_function(function, arity, context):
    """Accept native Lambda values and known public symbolic function heads."""
    from . import Dummy, Function, Lambda
    if isinstance(function, str):
        raise NotImplementedError(f'{context} requires a native Lambda; function strings are unsupported')
    if isinstance(function, _core.Expr):
        return function
    if isinstance(function, Function):
        if function._nargs is not None and arity not in function._nargs:
            raise ValueError(f'{context} function must accept {arity} arguments')
    else:
        head = getattr(getattr(function, 'func', None), '__name__', None)
        if head in ('Lambda', 'IdentityFunction'):
            return _core.sympify(function)
        # Compare identities with the official public wrappers. Ordinary
        # Python callables are never invoked to manufacture symbolic bodies.
        import sys
        namespace = vars(sys.modules[__package__])
        if arity == 2:
            names = ('atan2', 'Max', 'Min', 'Mod', 'log', 'beta', 'binomial',
                     'polygamma', 'besselj', 'bessely', 'besseli', 'besselk',
                     'hankel1', 'hankel2', 'jn', 'yn')
        else:
            names = ('exp', 'sin', 'cos', 'tan', 'cot', 'sec', 'csc', 'asin',
                     'acos', 'atan', 'acot', 'sinh', 'cosh', 'tanh', 'coth',
                     'asinh', 'acosh', 'atanh', 'Abs', 'sign', 'floor',
                     'ceiling', 'gamma', 'loggamma', 'erf', 'erfc', 're',
                     'im', 'conjugate', 'arg', 'Heaviside', 'DiracDelta',
                     'Max', 'Min', 'log', 'factorial', 'factorial2',
                     'subfactorial', 'harmonic', 'bernoulli', 'euler',
                     'catalan', 'fibonacci', 'lucas', 'partition', 'digamma',
                     'trigamma', 'zeta', 'dirichlet_eta', 'erfi', 'erfinv',
                     'erfcinv', 'fresnels', 'fresnelc', 'Ei', 'li', 'Li',
                     'Si', 'Ci', 'Shi', 'Chi', 'E1', 'airyai', 'airybi',
                     'airyaiprime', 'airybiprime', 'elliptic_k', 'elliptic_e',
                     'LambertW')
        if not any(function is candidate or function is getattr(candidate, '__wrapped__', None) for name in names
                   if (candidate := namespace.get(name)) is not None):
            raise TypeError(f'{context} requires a native Lambda or a supported symbolic function head')
    parameters = tuple(Dummy(name) for name in (('i', 'j') if arity == 2 else ('d',)))
    body = function(*parameters)
    return Lambda(parameters if arity == 2 else parameters[0], body)


class FunctionMatrix(MatrixExpr):
    _kinds = ('FunctionMatrix',)

    def __new__(cls, rows, cols, lamda):
        lamda = _normalize_matrix_function(lamda, 2, 'FunctionMatrix')
        return _core._matrixexpr_node('FunctionMatrix', (rows, cols, lamda))


class ElementwiseApplyFunction(MatrixExpr):
    _kinds = ('ElementwiseApplyFunction',)

    def __new__(cls, function, expr):
        function = _normalize_matrix_function(function, 1, 'ElementwiseApplyFunction')
        return _core._matrixexpr_node('ElementwiseApplyFunction', (function, expr))


class MatAdd(MatrixExpr):
    """Represent a shape-checked sum of symbolic matrices."""
    _kinds = ('MatAdd',)

    def __new__(cls, *args, evaluate=False, **options):
        if options:
            raise NotImplementedError('MatAdd options are not implemented: '+', '.join(options))
        return _core._matrixexpr_node('MatAdd', args, evaluate)


class MatMul(MatrixExpr):
    """Represent an ordered product of matrix and scalar operands."""
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

    def __new__(cls, mat, exp=-1):
        if exp != -1:
            raise ValueError('Inverse exponent must be -1')
        return _core._matrixexpr_node('Inverse',(mat,))


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


class OneMatrix(MatrixExpr):
    _kinds = ('OneMatrix',)

    def __new__(cls, m, n, evaluate=False):
        return _core._matrixexpr_node('OneMatrix', (m, n), evaluate)


class HadamardProduct(MatrixExpr):
    _kinds = ('HadamardProduct',)

    def __new__(cls, *args, evaluate=False, check=None):
        if check is not None and check is not True:
            raise NotImplementedError('Disabling Hadamard shape validation is unsupported')
        return _core._matrixexpr_node('HadamardProduct', args, evaluate)


class HadamardPower(MatrixExpr):
    _kinds = ('HadamardPower',)

    def __new__(cls, base, exp):
        return _core._matrixexpr_node('HadamardPower', (base, exp))


class KroneckerProduct(MatrixExpr):
    _kinds = ('KroneckerProduct',)

    def __new__(cls, *args, check=True):
        if check is not True:
            raise NotImplementedError('Disabling Kronecker operand validation is unsupported')
        return _core._matrixexpr_node('KroneckerProduct', args)


class DiagMatrix(MatrixExpr):
    _kinds = ('DiagMatrix',)

    def __new__(cls, vector):
        return _core._matrixexpr_node('DiagMatrix', (vector,))


class DiagonalOf(MatrixExpr):
    _kinds = ('DiagonalOf',)

    def __new__(cls, matrix):
        return _core._matrixexpr_node('DiagonalOf', (matrix,))


class DiagonalMatrix(MatrixExpr):
    _kinds = ('DiagonalMatrix',)

    def __new__(cls, matrix):
        return _core._matrixexpr_node('DiagonalMatrix', (matrix,))


class BlockMatrix(MatrixExpr):
    _kinds = ('BlockMatrix', 'BlockDiagMatrix')

    def __new__(cls, blocks):
        return _core._matrixexpr_node('BlockMatrix', (blocks,))


def block_collapse(expr):
    return _core._block_collapse(expr)


class BlockDiagMatrix(BlockMatrix):
    _kinds = ('BlockDiagMatrix',)

    def __new__(cls, *matrices):
        return _core._matrixexpr_node('BlockDiagMatrix', matrices)


class _SliceTuple(tuple):
    """A semantic slice triple retaining its immutable native normalized range."""
    def __new__(cls, spec):
        value=tuple.__new__(cls,spec.args)
        object.__setattr__(value,'_symox_slice_spec',spec)
        return value

    def __setattr__(self,name,value):
        raise AttributeError('Normalized slice specifications are immutable')

    def __reduce__(self):
        return (_SliceTuple,(self._symox_slice_spec,))


def NormalizedSlice(start,stop,step,length,omitted=False):
    """Reconstruct the normalized slice metadata used by native srepr/pickle."""
    return _SliceTuple(_core._normalized_slice(start,stop,step,length,omitted))


class MatrixSlice(MatrixExpr):
    _kinds = ('MatrixSlice',)

    def __new__(cls, parent, rowslice, colslice):
        return _core._matrixexpr_node('MatrixSlice', (parent, rowslice, colslice))


def hadamard_product(*matrices):
    if not matrices:
        raise TypeError('Empty Hadamard product is undefined')
    if len(matrices) == 1:
        return matrices[0]
    return HadamardProduct(*matrices).doit()


def hadamard_power(base, exp):
    from . import sympify
    base, exp = sympify(base), sympify(exp)
    if exp == 1:
        return base
    if not getattr(base, 'is_Matrix', False):
        return base ** exp
    if getattr(exp, 'is_Matrix', False):
        raise ValueError('cannot raise expression to a matrix')
    return HadamardPower(base, exp)


def kronecker_product(*matrices):
    if not matrices:
        raise TypeError('Empty Kronecker product is undefined')
    if len(matrices) == 1:
        return matrices[0]
    return KroneckerProduct(*matrices).doit()


class _ElementType(type):
    def __instancecheck__(cls, value):
        return isinstance(value, _core.Expr) and getattr(value.func,'__name__',None)=='MatrixElement'


class MatrixElement(metaclass=_ElementType):
    def __new__(cls, parent, i, j):
        return _core._matrix_element(parent,i,j)


class _ScalarMatrixType(type):
    def __instancecheck__(cls, value):
        return (isinstance(value, _core.Expr)
                and not getattr(value, 'is_UndefinedFunction', False)
                and getattr(value.func, '__name__', None) == cls.__name__)


class Trace(metaclass=_ScalarMatrixType):
    def __new__(cls, mat, evaluate=False):
        return _core._matrix_trace(mat, evaluate)


class Determinant(metaclass=_ScalarMatrixType):
    def __new__(cls, mat, evaluate=False):
        return _core._matrix_determinant(mat, evaluate)


class DotProduct(metaclass=_ScalarMatrixType):
    def __new__(cls, arg1, arg2):
        return _core._matrix_dot_product(arg1, arg2)


__all__ = ['DotProduct','MatrixExpr','MatrixSymbol','FunctionMatrix','ElementwiseApplyFunction','MatAdd','MatMul','MatPow','Transpose',
           'Adjoint','Inverse','Identity','ZeroMatrix','MatrixElement','MatrixConjugate',
           'Trace','Determinant','OneMatrix','HadamardProduct','HadamardPower',
           'KroneckerProduct','DiagMatrix','DiagonalOf','DiagonalMatrix','BlockDiagMatrix','BlockMatrix','block_collapse','MatrixSlice','hadamard_product','hadamard_power','kronecker_product']
