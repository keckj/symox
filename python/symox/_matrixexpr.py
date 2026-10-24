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


class _ClassOrInstanceMethod:
    """The metaclass method on the class (``MatMul.class_key()``) and the
    native method on an instance, as the classes now derive from the
    native type."""

    def __init__(self, name):
        self.name = name

    def __get__(self, obj, cls=None):
        if obj is not None and self.name == 'fromiter':
            obj, cls = None, obj.func
        if obj is None:
            return getattr(_MatrixType, self.name).__get__(cls, type(cls))
        return getattr(_core._MatrixExpr, self.name).__get__(obj, type(obj))


class MatrixExpr(_core._MatrixExpr, metaclass=_MatrixType):
    # The native nodes are allocated as the public class of their kind
    # (see _register_classes), which adds no storage.
    __slots__ = ()
    _kinds = None
    class_key = _ClassOrInstanceMethod('class_key')
    fromiter = _ClassOrInstanceMethod('fromiter')

    def __new__(cls, value):
        return _core._as_matrixexpr(value)

    @staticmethod
    def from_index_summation(expr, first_index=None, last_index=None, dimensions=None):
        from ._matrix_index_sums import from_index_summation
        return from_index_summation(expr, first_index, last_index, dimensions)


class MatrixSymbol(MatrixExpr):
    __slots__ = ()
    _kinds = ('MatrixSymbol',)

    def __new__(cls, name, n, m):
        return _core._matrixexpr_node('MatrixSymbol', (str(name), n, m))


def _normalize_matrix_function(function, arity, context):
    """Accept native Lambda values and known public symbolic function heads."""
    from . import Dummy, Function, Lambda
    if isinstance(function, str) and context != 'FunctionMatrix':
        raise NotImplementedError(f'{context} requires a native Lambda; function strings are unsupported')
    if isinstance(function, str):
        # FunctionMatrix sympifies the string; only a Lambda or a function
        # head is accepted.
        from . import sympify
        function = sympify(function)
        if isinstance(function, _core.Expr) and getattr(function.func, '__name__', None) != 'Lambda':
            raise ValueError(f'{function} should be compatible with SymPy function classes.')
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
            names = ('atan2', 'Max', 'Min', 'Mod', 'log', 'beta', 'binomial', 'KroneckerDelta',
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
    __slots__ = ()
    _kinds = ('FunctionMatrix',)

    def __new__(cls, rows, cols, lamda):
        lamda = _normalize_matrix_function(lamda, 2, 'FunctionMatrix')
        return _core._matrixexpr_node('FunctionMatrix', (rows, cols, lamda))


class ElementwiseApplyFunction(MatrixExpr):
    __slots__ = ()
    _kinds = ('ElementwiseApplyFunction',)

    def __new__(cls, function, expr):
        function = _normalize_matrix_function(function, 1, 'ElementwiseApplyFunction')
        return _core._matrixexpr_node('ElementwiseApplyFunction', (function, expr))

    def _get_function_fdiff(self):
        # SymPy's ElementwiseApplyFunction._get_function_fdiff
        from . import Dummy, Function, Lambda
        d = Dummy("d")
        function = self.function(d)
        fdiff = function.diff(d)
        if isinstance(fdiff, Function):
            fdiff = fdiff.func
        else:
            fdiff = Lambda(d, fdiff)
        return fdiff


class MatAdd(MatrixExpr):
    """Represent a shape-checked sum of symbolic matrices."""
    __slots__ = ()
    _kinds = ('MatAdd',)

    def __new__(cls, *args, evaluate=False, **options):
        if options:
            raise NotImplementedError('MatAdd options are not implemented: '+', '.join(options))
        if not args:
            return GenericZeroMatrix()
        # MatAdd.__new__ drops its identity GenericZeroMatrix() first.
        args = tuple(a for a in args if not _is_generic(a, 'GenericZeroMatrix'))
        if not args:
            # validate(*args) with nothing left to unpack.
            raise ValueError('not enough values to unpack (expected 2, got 0)')
        return _core._matrixexpr_node('MatAdd', args, evaluate)


class MatMul(MatrixExpr):
    """Represent an ordered product of matrix and scalar operands."""
    __slots__ = ()
    _kinds = ('MatMul',)

    def __new__(cls, *args, evaluate=False, **options):
        if options:
            raise NotImplementedError('MatMul options are not implemented: '+', '.join(options))
        if not args:
            return GenericIdentity()
        # MatMul.__new__ drops its identity GenericIdentity() first.
        args = tuple(a for a in args if not _is_generic(a, 'GenericIdentity'))
        return _core._matrixexpr_node('MatMul', args, evaluate)


class MatPow(MatrixExpr):
    __slots__ = ()
    _kinds = ('MatPow','Inverse')

    def __new__(cls, base, exp, evaluate=False):
        return _core._matrixexpr_node('MatPow',(base,exp),evaluate)


class Inverse(MatPow):
    __slots__ = ()
    _kinds = ('Inverse',)

    def __new__(cls, mat, exp=-1):
        if exp != -1:
            raise ValueError('Inverse exponent must be -1')
        return _core._matrixexpr_node('Inverse',(mat,))


class Transpose(MatrixExpr):
    __slots__ = ()
    _kinds = ('Transpose',)

    def __new__(cls, arg):
        return _core._matrixexpr_node('Transpose',(arg,))


class Adjoint(MatrixExpr):
    __slots__ = ()
    _kinds = ('Adjoint',)

    def __new__(cls, arg):
        return _core._matrixexpr_node('Adjoint',(arg,))


class MatrixConjugate(MatrixExpr):
    __slots__ = ()
    _kinds = ('MatrixConjugate',)

    def __new__(cls, arg):
        return _core._matrixexpr_node('MatrixConjugate',(arg,))


class Identity(MatrixExpr):
    __slots__ = ()
    _kinds = ('Identity', 'GenericIdentity')

    def __new__(cls, n):
        return _core._matrixexpr_node('Identity',(n,))


class ZeroMatrix(MatrixExpr):
    __slots__ = ()
    _kinds = ('ZeroMatrix', 'GenericZeroMatrix')

    def __new__(cls, n, m):
        return _core._matrixexpr_node('ZeroMatrix',(n,m))


def _is_generic(value, tag):
    return isinstance(value, _core._MatrixExpr) and value._kind_tag == tag


def _shapeless(name):
    def getter(self):
        raise TypeError('%s does not have a specified shape' % name)
    return property(getter)


class GenericIdentity(Identity):
    """The shapeless identity of MatMul (``MatMul()``)."""
    __slots__ = ()
    _kinds = ('GenericIdentity',)
    is_Identity = True
    shape = rows = cols = _shapeless('GenericIdentity')
    func = property(lambda self: type(self))

    def __new__(cls):
        return _core._matrixexpr_node('GenericIdentity', ())

    # MatrixExpr.__mul__ is MatMul(self, other).doit(): MatMul drops the
    # GenericIdentity factors and a product left without matrices is its
    # scalar factor (``2*GenericIdentity()`` is 2, ``G*G`` is 1).
    def __mul__(self, other):
        if _is_generic(other, 'GenericIdentity'):
            return _core.sympify(1)
        return _core.sympify(other)

    __rmul__ = __mul__
    __matmul__ = __rmatmul__ = __mul__


class GenericZeroMatrix(ZeroMatrix):
    """The shapeless identity of MatAdd (``MatAdd()``)."""
    __slots__ = ()
    _kinds = ('GenericZeroMatrix',)
    is_ZeroMatrix = True
    shape = rows = cols = _shapeless('GenericZeroMatrix')
    func = property(lambda self: type(self))

    def __new__(cls):
        return _core._matrixexpr_node('GenericZeroMatrix', ())


class DFT(MatrixExpr):
    """The unitary discrete Fourier transform matrix (fourier.py)."""
    __slots__ = ()
    _kinds = ('DFT', 'IDFT')

    def __new__(cls, n):
        return _core._matrixexpr_node('DFT', (n,))

    n = property(lambda self: self.args[0])
    # Not exported at top level: resolve func through the class itself.
    func = property(lambda self: type(self))


class IDFT(DFT):
    """The inverse of ``DFT(n)``."""
    __slots__ = ()
    _kinds = ('IDFT',)

    def __new__(cls, n):
        return _core._matrixexpr_node('IDFT', (n,))


class OneMatrix(MatrixExpr):
    __slots__ = ()
    _kinds = ('OneMatrix',)

    def __new__(cls, m, n, evaluate=False):
        return _core._matrixexpr_node('OneMatrix', (m, n), evaluate)


class HadamardProduct(MatrixExpr):
    __slots__ = ()
    _kinds = ('HadamardProduct',)

    def __new__(cls, *args, evaluate=False, check=None):
        if check is not None and check is not True:
            raise NotImplementedError('Disabling Hadamard shape validation is unsupported')
        return _core._matrixexpr_node('HadamardProduct', args, evaluate)


class HadamardPower(MatrixExpr):
    __slots__ = ()
    _kinds = ('HadamardPower',)

    def __new__(cls, base, exp):
        return _core._matrixexpr_node('HadamardPower', (base, exp))


class KroneckerProduct(MatrixExpr):
    __slots__ = ()
    _kinds = ('KroneckerProduct',)

    def __new__(cls, *args, check=True):
        if check is not True:
            raise NotImplementedError('Disabling Kronecker operand validation is unsupported')
        return _core._matrixexpr_node('KroneckerProduct', args)


class DiagMatrix(MatrixExpr):
    __slots__ = ()
    _kinds = ('DiagMatrix',)

    def __new__(cls, vector):
        return _core._matrixexpr_node('DiagMatrix', (vector,))


class DiagonalOf(MatrixExpr):
    __slots__ = ()
    _kinds = ('DiagonalOf',)

    def __new__(cls, matrix):
        return _core._matrixexpr_node('DiagonalOf', (matrix,))


class DiagonalMatrix(MatrixExpr):
    __slots__ = ()
    _kinds = ('DiagonalMatrix',)

    def __new__(cls, matrix):
        return _core._matrixexpr_node('DiagonalMatrix', (matrix,))


def _check_block_regularity(rows):
    """SymPy's BlockMatrix regularity check, raising its ValueError texts."""
    is_mat = lambda i: getattr(i, 'is_Matrix', False)
    if not isinstance(rows, (list, tuple)) or not rows:
        return
    if len({is_mat(r) for r in rows}) != 1:
        raise ValueError('\nexpecting a sequence of 1 or more rows containing Matrices.')
    if is_mat(rows[0]):
        rows = [rows]
    if not all(isinstance(r, (list, tuple)) and all(is_mat(i) for i in r) for r in rows):
        return
    blocky = ok = len({len(r) for r in rows}) == 1
    if ok:
        for r in rows:
            ok = len({i.rows for i in r}) == 1
            if not ok:
                break
        blocky = ok
        if ok:
            for c in range(len(rows[0])):
                ok = len({rows[i][c].cols for i in range(len(rows))}) == 1
                if not ok:
                    break
    if not ok:
        ok = len({sum(i.cols for i in r) for r in rows}) == 1
        if blocky and ok:
            raise ValueError(
                '\nAlthough this matrix is comprised of blocks, the blocks do not fill\n'
                'the matrix in a size-symmetric fashion. To create a full matrix from\n'
                'these arguments, pass them directly to Matrix.')
        raise ValueError(
            "\nWhen there are not the same number of rows in each row's matrices or\n"
            'there are not the same number of total columns in each row, the matrix\n'
            'is not a block matrix. If this matrix is known to consist of blocks\n'
            'fully filling a 2-D space then see Matrix.irregular.')


class BlockMatrix(MatrixExpr):
    __slots__ = ()
    _kinds = ('BlockMatrix', 'BlockDiagMatrix')

    def __new__(cls, blocks):
        _check_block_regularity(blocks)
        return _core._matrixexpr_node('BlockMatrix', (blocks,))


def block_collapse(expr):
    return _core._block_collapse(expr)


class BlockDiagMatrix(BlockMatrix):
    __slots__ = ()
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
    __slots__ = ()
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


# SymPy's defining modules of the public classes.
_SYMPY_MODULES = {
    'MatrixExpr': 'matexpr', 'MatrixSymbol': 'matexpr', 'FunctionMatrix': 'funcmatrix',
    'ElementwiseApplyFunction': 'applyfunc', 'MatAdd': 'matadd', 'MatMul': 'matmul',
    'MatPow': 'matpow', 'Inverse': 'inverse', 'Transpose': 'transpose', 'Adjoint': 'adjoint',
    'Identity': 'special', 'ZeroMatrix': 'special', 'OneMatrix': 'special',
    'GenericIdentity': 'special', 'GenericZeroMatrix': 'special', 'DFT': 'fourier', 'IDFT': 'fourier',
    'HadamardProduct': 'hadamard', 'HadamardPower': 'hadamard', 'KroneckerProduct': 'kronecker',
    'DiagMatrix': 'diagonal', 'DiagonalOf': 'diagonal', 'DiagonalMatrix': 'diagonal',
    'BlockMatrix': 'blockmatrix', 'BlockDiagMatrix': 'blockmatrix', 'MatrixSlice': 'slice',
    'PermutationMatrix': 'permutation', 'MatrixPermute': 'permutation',
}


def _register_classes(*classes):
    """Allocate the native nodes of each kind as its public class."""
    for cls in classes:
        kinds = cls.__dict__.get('_kinds')
        if kinds and cls.__name__ in kinds:
            _core._register_matrixexpr_type(cls.__name__, cls)
        module = _SYMPY_MODULES.get(cls.__name__)
        if module is not None:
            cls.__module__ = 'symox.matrices.expressions.' + module


if hasattr(_core, '_register_matrixexpr_type'):
    _register_classes(MatrixExpr, MatrixSymbol, FunctionMatrix, ElementwiseApplyFunction, MatAdd, MatMul,
                      MatPow, Inverse, Transpose, Adjoint, MatrixConjugate, Identity, ZeroMatrix, OneMatrix,
                      GenericIdentity, GenericZeroMatrix, DFT, IDFT,
                      HadamardProduct, HadamardPower, KroneckerProduct, DiagMatrix, DiagonalOf,
                      DiagonalMatrix, BlockMatrix, BlockDiagMatrix, MatrixSlice)


__all__ = ['DotProduct','MatrixExpr','MatrixSymbol','FunctionMatrix','ElementwiseApplyFunction','MatAdd','MatMul','MatPow','Transpose',
           'Adjoint','Inverse','Identity','ZeroMatrix','MatrixElement','MatrixConjugate',
           'Trace','Determinant','OneMatrix','HadamardProduct','HadamardPower',
           'KroneckerProduct','DiagMatrix','DiagonalOf','DiagonalMatrix','BlockDiagMatrix','BlockMatrix','block_collapse','MatrixSlice','hadamard_product','hadamard_power','kronecker_product']

MatMul.identity = GenericIdentity()
MatAdd.identity = GenericZeroMatrix()
