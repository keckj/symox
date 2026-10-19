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
    reflected = {'__radd__', '__rsub__', '__rmul__', '__rmatmul__'}

    def result_class(self, args, name):
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
            return reconstruct(result, result_class(self, args, name))
        return method

    # Conversions, derivative arrays, constructors and mutators have separate
    # category rules. Wrap the operations qualified against the reference.
    names = arithmetic | reflected | {'__pow__','__neg__','__pos__',
        '__truediv__','__getitem__','transpose','conjugate','adjoint','extract',
        'reshape','row','col','applyfunc','subs','xreplace','copy','inv','rref',
        'cholesky'}
    for name in names:
        setattr(matrix, name, wrap(name, matrix.__dict__[name]))
    for name in ['T', 'H', 'C']:
        descriptor = matrix.__dict__[name]
        setattr(matrix, name, property(wrap(name, descriptor.__get__)))
