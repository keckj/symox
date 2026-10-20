"""Reconstruction of mutable dense matrix subclasses after native kernels."""
from functools import wraps

from . import _core


def install():
    matrix = _core.Matrix
    native_new = matrix.__new__

    @classmethod
    def new(cls, *args, copy=True, **kwargs):
        if not copy and len(args) != 3:
            raise TypeError("'copy=False' requires a matrix be initialized as rows,cols,[list]")
        # Bypass a user's public constructor, as the reference _new protocol
        # does. The native allocator still allocates the real Python subclass.
        return native_new(cls, *args, **kwargs)

    matrix._new = new
    matrix._class_priority = 4

    def mutable_dense(value):
        return isinstance(value, matrix) and not value._immutable

    def reconstruct(value, cls):
        if cls is matrix:
            return value
        if isinstance(value, matrix):
            if value._immutable or type(value) is cls:
                return value
            return cls._new(value.rows, value.cols,
                            [entry for row in matrix.tolist(value) for entry in row])
        if isinstance(value, tuple):
            return tuple(reconstruct(item, cls) for item in value)
        if isinstance(value, list):
            return [reconstruct(item, cls) for item in value]
        return value

    arithmetic = {'__add__', '__sub__', '__mul__', '__matmul__'}
    right_class = {'LUsolve'}

    def priority(value):
        # Immutable dense matrices are a flag of the native class; the
        # reference gives ImmutableDenseMatrix priority 8.
        return 8 if value._immutable else value._class_priority

    def classof(self, other):
        # The reference classof: the right operand wins ties.
        return self if priority(self) > priority(other) else other
    reflected = {'__radd__', '__rsub__', '__rmul__', '__rmatmul__'}

    def result_class(self, args, name):
        if name in right_class and args and isinstance(args[0], matrix):
            return type(args[0]) if mutable_dense(args[0]) else matrix
        if not mutable_dense(self):
            return matrix
        cls = type(self)
        if name in arithmetic | reflected and args and mutable_dense(args[0]):
            left, right = (args[0], self) if name in reflected else (self, args[0])
            # Addition unifies twice in the reference, retaining the left
            # class on equal priority. Multiplication retains the right one.
            left_wins = (left._class_priority >= right._class_priority
                         if name in {'__add__','__radd__','__sub__','__rsub__'}
                         else left._class_priority > right._class_priority)
            cls = type(left) if left_wins else type(right)
        return cls

    def wrap(name, original):
        @wraps(original)
        def method(self, *args, **kwargs):
            result = original(self, *args, **kwargs)
            if name == 'multiply_elementwise' and args and isinstance(args[0], matrix) and isinstance(result, matrix):
                winner = classof(self, args[0])
                result = result.as_immutable() if winner._immutable else result.as_mutable()
                return reconstruct(result, type(winner) if mutable_dense(winner) else matrix)
            if name in right_class and args and isinstance(args[0], matrix) and isinstance(result, matrix):
                # The operand's class also decides mutability here.
                result = result.as_immutable() if args[0]._immutable else result.as_mutable()
            return reconstruct(result, result_class(self, args, name))
        return method

    # Conversions, derivative arrays, constructors and mutators have separate
    # category rules. Wrap the operations qualified against the reference.
    names = arithmetic | reflected | {'__pow__','__neg__','__pos__',
        '__truediv__','__getitem__','transpose','conjugate','adjoint','extract',
        'reshape','row','col','applyfunc','subs','xreplace','copy','inv','rref',
        'cholesky',
        # The reference builds these results with self._new or with the
        # class unifying self and a matrix operand.
        'expand','evalf','n','adjugate','jacobian','cofactor_matrix','exp','log',
        'analytic_func','col_insert','row_insert','col_join','row_join',
        'diagonal','upper_triangular','lower_triangular','minor_submatrix',
        'cross','normalized','pinv','gauss_jordan_solve','solve','rot90',
        'permute','vec','vech','as_real_imag','multiply_elementwise',
        'nullspace','columnspace','rowspace','eigenvects','left_eigenvects',
        'diagonalize','LUdecomposition','QRdecomposition','jordan_form',
        'inverse_ADJ','inverse_LU','inverse_GE','echelon_form','LUsolve',
        'lower_triangular_solve','upper_triangular_solve',
        'cholesky_solve','LDLsolve','LDLdecomposition',
        'singular_value_decomposition','elementary_row_op','elementary_col_op'}
    for name in names:
        setattr(matrix, name, wrap(name, matrix.__dict__[name]))

    original_rank_decomposition = matrix.__dict__['rank_decomposition']

    @wraps(original_rank_decomposition)
    def rank_decomposition(self, *args, **kwargs):
        # The reference builds C from the input class and F from rref,
        # whose result is a plain mutable dense matrix.
        C, F = original_rank_decomposition(self, *args, **kwargs)
        return reconstruct(C, result_class(self, args, 'rank_decomposition')), F
    matrix.rank_decomposition = rank_decomposition

    def constructor(name):
        original = matrix.__dict__[name].__func__
        @wraps(original)
        def build(klass, /, *args, **kwargs):
            if name == 'diag' and 'cls' in kwargs:
                # The reference builds the result with the requested class.
                target = kwargs.pop('cls')
                result = original(*args, **kwargs)
                return reconstruct(result, target) if isinstance(target, type) and issubclass(target, matrix) else target(result)
            result = original(*args, **kwargs)
            if name in {'hstack', 'vstack'}:
                # The reference reduces with type(args[0]).row_join.
                target = type(args[0]) if args and mutable_dense(args[0]) else (klass if not args else matrix)
            else:
                target = klass
            return reconstruct(result, target) if isinstance(target, type) and issubclass(target, matrix) else result
        return classmethod(build)

    for name in ('zeros', 'eye', 'ones', 'diag', 'hstack', 'vstack'):
        setattr(matrix, name, constructor(name))
    for name in ['T', 'H', 'C']:
        descriptor = matrix.__dict__[name]
        setattr(matrix, name, property(wrap(name, descriptor.__get__)))
