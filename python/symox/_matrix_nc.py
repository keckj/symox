"""Matrix algorithms on noncommutative entries, as in the reference.

The characteristic polynomial and every algorithm built on it (eigenvalues,
diagonalization, Jordan form, exponential, singular values, pseudoinverse,
nilpotency) build a ``Poly`` of the entries in the reference, which raises
``PolynomialError`` for noncommutative entries. The native kernels would
otherwise treat such entries as commuting and return a wrong value.

The inverse by Gauss elimination or by the adjugate and ``solve`` are
computed by the LU method for noncommutative entries: the reference
methods apply commutative formulas (``A*(A**2 - B*C)**-1`` for the block
``[[A, B], [C, A]]``), and only the LU inverse is the inverse of the block
matrix. Echelon forms with noncommutative entries are not supported.
"""
import functools

from . import _core

_POLY = ('charpoly', 'eigenvals', 'eigenvects', 'is_diagonalizable', 'diagonalize',
         'jordan_form', 'exp', 'singular_values', 'pinv', 'is_nilpotent', 'log',
         'analytic_func', 'left_eigenvects', 'bidiagonal_decomposition')


def _is_nc(matrix):
    try:
        entries = list(matrix)
    except Exception:
        return False
    return any(isinstance(e, _core.Expr) and _core._has_noncommutative_atoms(e) for e in entries)


def _polynomial_error():
    from .polys.polyerrors import PolynomialError
    return PolynomialError('non-commutative expressions are not supported')


def _guard_poly(native):
    @functools.wraps(native)
    def method(self, *args, **kwargs):
        if _is_nc(self):
            raise _polynomial_error()
        return native(self, *args, **kwargs)
    return method


def _guard_inv(native):
    @functools.wraps(native)
    def inv(self, method=None, iszerofunc=None, try_block_diag=False):
        if method in ('GE', 'ADJ') and _is_nc(self):
            method = 'LU'
        kwargs = {}
        if iszerofunc is not None:
            kwargs['iszerofunc'] = iszerofunc
        if try_block_diag:
            kwargs['try_block_diag'] = try_block_diag
        return native(self, method=method, **kwargs)
    return inv


def _guard_solve(native):
    @functools.wraps(native)
    def solve(self, rhs, method='GJ'):
        if _is_nc(self) or _is_nc(rhs):
            if not self.is_square:
                raise NotImplementedError('solve with noncommutative entries needs a square matrix')
            return self.LUsolve(rhs)
        return native(self, rhs, method=method)
    return solve


def _guard_unsupported(native, name):
    @functools.wraps(native)
    def method(self, *args, **kwargs):
        if _is_nc(self):
            raise NotImplementedError('%s of a matrix with noncommutative entries is not implemented' % name)
        return native(self, *args, **kwargs)
    return method


def install(*classes):
    for cls in classes:
        for name in _POLY:
            native = getattr(cls, name, None)
            if native is not None:
                setattr(cls, name, _guard_poly(native))
        if getattr(cls, 'inv', None) is not None:
            cls.inv = _guard_inv(cls.inv)
        if getattr(cls, 'solve', None) is not None:
            cls.solve = _guard_solve(cls.solve)
        for name in ('echelon_form', 'QRdecomposition', 'QRsolve'):
            native = getattr(cls, name, None)
            if native is not None:
                setattr(cls, name, _guard_unsupported(native, name))
