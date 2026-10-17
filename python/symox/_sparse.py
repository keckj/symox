"""Dictionary-backed symbolic sparse matrices.

Construction and ordinary algebra operate on stored nonzero entries. Advanced
factorizations explicitly use the native dense matrix algorithms and convert
matrix results back to sparse storage.
"""
from __future__ import annotations

import operator
from types import MappingProxyType

from ._core import Matrix, sympify
from . import ShapeError, NonSquareMatrixError

_ZERO = sympify(0)
_ONE = sympify(1)


def _zero_product(value):
    # Held sums/products can hide an explicit infinity from scalar fast paths.
    if value.is_Add or value.is_Mul:
        if any(_zero_product(arg) != 0 for arg in value.args):
            return sympify('nan')
    return _ZERO*value


def _integer(value):
    try:
        return operator.index(value)
    except TypeError:
        if getattr(value, 'is_Integer', False):
            return int(value)
        raise TypeError('Matrix dimensions and indices must be integers') from None


def _dimension(value):
    value = _integer(value)
    if value < 0:
        raise ValueError('Matrix dimensions must be nonnegative')
    return value


def _index(value, limit):
    value = _integer(value)
    if value < 0:
        value += limit
    if not 0 <= value < limit:
        raise IndexError('Matrix index out of range')
    return value


def _matrix_like(value):
    return hasattr(value, 'rows') and hasattr(value, 'cols') and hasattr(value, '__getitem__')


def _entries(value):
    if isinstance(value, _SparseBase):
        return dict(value._data)
    if hasattr(value, 'todok'):
        return {(int(i), int(j)): sympify(v) for (i, j), v in value.todok().items() if v != 0}
    return {(i, j): sympify(value[i, j]) for i in range(int(value.rows))
            for j in range(int(value.cols)) if value[i, j] != 0}


def _creation(args):
    if not args:
        return 0, 0, {}
    if len(args) == 1:
        value = args[0]
        if _matrix_like(value):
            return int(value.rows), int(value.cols), _entries(value)
        if hasattr(value, 'tolist'):
            shape = getattr(value, 'shape', ())
            if len(shape) == 2 and not all(shape):
                return int(shape[0]), int(shape[1]), {}
            value = value.tolist()
        try:
            rows = list(value)
        except TypeError:
            raise TypeError('Expected a matrix or sequence') from None
        nested = [list(row) if isinstance(row, (list, tuple)) else [row] for row in rows]
        cols = max((len(row) for row in nested), default=0)
        data = {(i, j): sympify(v) for i, row in enumerate(nested)
                for j, v in enumerate(row) if sympify(v) != 0}
        return len(nested) if cols else 0, cols, data
    if len(args) == 2 and args[0] is None:
        args = (None, None, args[1])
    if len(args) != 3:
        raise TypeError('Expected (), (matrix), or (rows, cols, entries)')
    rows, cols, source = args
    automatic = rows is None and cols is None
    if not automatic:
        if rows is None or cols is None:
            raise ValueError('Both dimensions must be None for autosizing')
        rows, cols = _dimension(rows), _dimension(cols)
    data = {}

    def store(i, j, value):
        i, j = _integer(i), _integer(j)
        if i < 0 or j < 0:
            raise ValueError('Sparse coordinates must be nonnegative')
        if not automatic and (i >= rows or j >= cols):
            raise ValueError('Sparse coordinate is outside matrix shape')
        value = sympify(value)
        if value != 0:
            if (i, j) in data and data[i, j] != value:
                raise ValueError('Overlapping sparse entries disagree')
            data[i, j] = value

    if hasattr(source, 'items'):
        for key, value in source.items():
            if not isinstance(key, tuple) or len(key) != 2:
                raise TypeError('Sparse entries must use (row, col) keys')
            i, j = key
            if _matrix_like(value) or isinstance(value, (list, tuple)):
                block = MutableSparseMatrix(value)
                for (bi, bj), entry in block._data.items():
                    store(i + bi, j + bj, entry)
            else:
                store(i, j, value)
    elif callable(source):
        if automatic:
            raise ValueError('Callable construction requires dimensions')
        for i in range(rows):
            for j in range(cols):
                store(i, j, source(i, j))
    else:
        if automatic:
            raise ValueError('Autosizing requires a coordinate dictionary')
        source = list(source)
        if any(isinstance(row, (list, tuple)) for row in source):
            r, c, source = _creation((source,))
            if (r, c) != (rows, cols):
                raise ValueError('Sequence does not match matrix shape')
            for (i, j), value in source.items():
                store(i, j, value)
        else:
            if len(source) != rows * cols:
                raise ValueError('Sequence does not match matrix shape')
            for pos, value in enumerate(source):
                store(pos // cols, pos % cols, value)
    if automatic:
        rows = max((i for i, _ in data), default=-1) + 1
        cols = max((j for _, j in data), default=-1) + 1
    return rows, cols, data


class _SparseBase:
    __slots__ = ('_rows', '_cols', '_data', '_frozen')
    is_Matrix = True
    is_MatrixExpr = False
    _symox_sparse = True
    _immutable = False
    is_Immutable = False
    __array_priority__ = 10000
    __hash__ = None

    def __init__(self, *args):
        if getattr(self, '_frozen', False):
            raise TypeError('Cannot reinitialize an immutable sparse matrix')
        rows, cols, data = _creation(args)
        object.__setattr__(self, '_rows', rows)
        object.__setattr__(self, '_cols', cols)
        object.__setattr__(self, '_data', MappingProxyType(data) if self._immutable else data)
        object.__setattr__(self, '_frozen', self._immutable)

    def __setattr__(self, name, value):
        if getattr(self, '_frozen', False):
            raise TypeError('Cannot modify an immutable sparse matrix')
        object.__setattr__(self, name, value)

    def __delattr__(self, name):
        if getattr(self, '_frozen', False):
            raise TypeError('Cannot modify an immutable sparse matrix')
        object.__delattr__(self, name)

    @property
    def rows(self): return self._rows
    @property
    def cols(self): return self._cols
    @property
    def shape(self): return self._rows, self._cols
    @property
    def is_square(self): return self._rows == self._cols
    @property
    def is_upper(self): return all(i <= j for i, j in self._data)
    @property
    def is_lower(self): return all(i >= j for i, j in self._data)
    @property
    def is_zero_matrix(self):
        if not self._data: return True
        return False if any(v.is_zero is False for v in self._data.values()) else None
    @property
    def is_hermitian(self): return self.as_dense().is_hermitian
    @property
    def free_symbols(self):
        result = set()
        for value in self._data.values(): result.update(value.free_symbols)
        return result
    @property
    def T(self): return self.transpose()
    @property
    def H(self): return self.conjugate().transpose()
    @property
    def C(self): return self.conjugate()

    def __len__(self): return self.rows * self.cols
    def __bool__(self): return bool(len(self))
    def __iter__(self):
        for i in range(self.rows):
            for j in range(self.cols): yield self._data.get((i, j), _ZERO)

    def _new(self, rows, cols, data, peer=None):
        cls = ImmutableSparseMatrix if self._immutable or getattr(peer, '_immutable', False) else MutableSparseMatrix
        return cls(rows, cols, data)

    def __getitem__(self, key):
        if not isinstance(key, tuple):
            if isinstance(key, slice):
                return [self[pos] for pos in range(*key.indices(len(self))) ]
            pos = _index(key, len(self))
            return self._data.get(divmod(pos, self.cols), _ZERO)
        if len(key) != 2: raise IndexError('Expected two matrix indices')
        i, j = key
        if not isinstance(i, slice) and not isinstance(j, slice):
            return self._data.get((_index(i, self.rows), _index(j, self.cols)), _ZERO)
        ri = list(range(*i.indices(self.rows))) if isinstance(i, slice) else [_index(i, self.rows)]
        cj = list(range(*j.indices(self.cols))) if isinstance(j, slice) else [_index(j, self.cols)]
        rm, cm = {v:k for k,v in enumerate(ri)}, {v:k for k,v in enumerate(cj)}
        data = {(rm[r], cm[c]): v for (r,c),v in self._data.items() if r in rm and c in cm}
        return self._new(len(ri), len(cj), data)

    def _mutable(self):
        if self._immutable: raise TypeError('Cannot modify an immutable sparse matrix')

    def __setitem__(self, key, value):
        self._mutable()
        if not isinstance(key, tuple):
            if isinstance(key, slice):
                positions = list(range(*key.indices(len(self))))
                values = list(value)
                if len(positions) != len(values): raise ShapeError('Slice assignment size mismatch')
                for pos, v in zip(positions, values): self[pos] = v
                return
            key = divmod(_index(key, len(self)), self.cols)
        i, j = key
        if isinstance(i, slice) or isinstance(j, slice):
            ri = list(range(*i.indices(self.rows))) if isinstance(i, slice) else [_index(i, self.rows)]
            cj = list(range(*j.indices(self.cols))) if isinstance(j, slice) else [_index(j, self.cols)]
            block = MutableSparseMatrix(value)
            if block.shape != (len(ri), len(cj)): raise ShapeError('Slice assignment shape mismatch')
            for r in ri:
                for c in cj: self._data.pop((r,c), None)
            for (r,c),v in block._data.items(): self._data[ri[r],cj[c]] = v
            return
        key = _index(i, self.rows), _index(j, self.cols)
        value = sympify(value)
        if value == 0: self._data.pop(key, None)
        else: self._data[key] = value

    def todok(self): return dict(self._data)
    def nnz(self): return len(self._data)
    def row_list(self): return [(i,j,self._data[i,j]) for i,j in sorted(self._data)]
    def col_list(self): return [(i,j,self._data[i,j]) for i,j in sorted(self._data, key=lambda k:(k[1],k[0]))]
    def values(self): return [v for _,_,v in self.row_list()]
    def iter_values(self): return iter(self.values())
    def iter_items(self): return iter(sorted(self._data.items()))
    def tolist(self): return [[self._data.get((i,j), _ZERO) for j in range(self.cols)] for i in range(self.rows)]
    def as_dense(self):
        result = Matrix(self.rows, self.cols, list(self))
        return result.as_immutable() if self._immutable else result
    todense = as_dense
    def as_mutable(self): return MutableSparseMatrix(self)
    def as_immutable(self): return ImmutableSparseMatrix(self)
    def copy(self): return type(self)(self)

    def to_sympy(self):
        import sympy
        cls = sympy.ImmutableSparseMatrix if self._immutable else sympy.SparseMatrix
        return cls(self.rows,self.cols,{k:v.to_sympy() for k,v in self._data.items()})
    _sympy_ = to_sympy
    def __reduce__(self): return type(self), (self.rows,self.cols,self.todok())
    def __repr__(self): return repr(self.as_dense())
    def __str__(self): return str(self.as_dense())
    def __array__(self, dtype=None, copy=None):
        import numpy as np
        if copy is False: raise ValueError('Sparse matrices require a copy to make a dense array')
        return np.array(self.tolist(),dtype=dtype or object).reshape(self.shape)

    def __eq__(self, other):
        if not _matrix_like(other): return False
        return self.shape == (int(other.rows),int(other.cols)) and self.todok() == _entries(other)
    def __ne__(self, other): return not self == other
    def equals(self, other, failing_expression=False):
        if not _matrix_like(other) or self.shape != (int(other.rows),int(other.cols)): return False
        difference = self - other
        result = True
        for entry in difference.values():
            value = entry.equals(0)
            if value is False: return False
            if value is None: result = None
        return result

    def _add(self, other, sign):
        if not _matrix_like(other): return NotImplemented
        if self.shape != (int(other.rows),int(other.cols)): raise ShapeError('Matrix size mismatch')
        if not isinstance(other, _SparseBase):
            return self.as_dense() + other if sign == 1 else self.as_dense() - other
        data = self.todok()
        for key, value in other._data.items():
            result = data.get(key, _ZERO) + sign * value
            if result == 0: data.pop(key, None)
            else: data[key] = result
        return self._new(self.rows,self.cols,data,other)
    def __add__(self, other): return self._add(other,1)
    def __radd__(self, other): return self._add(other,1)
    def __sub__(self, other): return self._add(other,-1)
    def __rsub__(self, other): return (-self)._add(other,1)
    def __neg__(self): return self._new(self.rows,self.cols,{k:-v for k,v in self._data.items()})
    def __pos__(self): return self.copy()

    def __mul__(self, other):
        if _matrix_like(other):
            if self.cols != int(other.rows): raise ShapeError('Matrix size mismatch')
            if not isinstance(other, _SparseBase): return self.as_dense() * other
            right_rows = {}
            for (i,j),value in other._data.items(): right_rows.setdefault(i,[]).append((j,value))
            data = {}
            for (i,k),left in self._data.items():
                for j,right in right_rows.get(k,()):
                    data[i,j] = data.get((i,j),_ZERO) + left*right
            # Zero annihilates finite factors, but 0*oo/0*zoo/0*nan is nan.
            # Fill only affected output rows/columns; neither input is densified.
            for (i,k),left in self._data.items():
                zero_product = _zero_product(left)
                if zero_product != 0:
                    for j in range(other.cols):
                        if (k,j) not in other._data:
                            data[i,j] = data.get((i,j),_ZERO) + zero_product
            for (k,j),right in other._data.items():
                zero_product = _zero_product(right)
                if zero_product != 0:
                    for i in range(self.rows):
                        if (i,k) not in self._data:
                            data[i,j] = data.get((i,j),_ZERO) + zero_product
            return self._new(self.rows,other.cols,{k:v for k,v in data.items() if v != 0},other)
        try: scalar = sympify(other)
        except (TypeError,ValueError): return NotImplemented
        if _ZERO*scalar != 0: return self.applyfunc(lambda value:value*scalar)
        return self._new(self.rows,self.cols,{k:v*scalar for k,v in self._data.items() if v*scalar != 0})
    def __rmul__(self, other):
        if _matrix_like(other): return other * self.as_dense()
        try: scalar = sympify(other)
        except (TypeError,ValueError): return NotImplemented
        if scalar*_ZERO != 0: return self.applyfunc(lambda value:scalar*value)
        return self._new(self.rows,self.cols,{k:scalar*v for k,v in self._data.items() if scalar*v != 0})
    def __matmul__(self, other): return self * other if _matrix_like(other) else NotImplemented
    def __rmatmul__(self, other): return self.__rmul__(other) if _matrix_like(other) else NotImplemented
    def __truediv__(self, other):
        if _matrix_like(other): return NotImplemented
        scalar = sympify(other)
        if _ZERO/scalar != 0: return self.applyfunc(lambda value:value/scalar)
        return self._new(self.rows,self.cols,{k:v/scalar for k,v in self._data.items() if v/scalar != 0})
    def __pow__(self, exponent, modulo=None):
        if modulo is not None: raise TypeError('Modular matrix powers are unsupported')
        exponent = _integer(exponent)
        if not self.is_square: raise NonSquareMatrixError('Matrix must be square')
        base = self
        if exponent < 0: base, exponent = self.inv(), -exponent
        result = type(self).eye(self.rows)
        while exponent:
            if exponent & 1: result = result*base
            exponent >>= 1
            if exponent: base = base*base
        return result

    def transpose(self): return self._new(self.cols,self.rows,{(j,i):v for (i,j),v in self._data.items()})
    def conjugate(self):
        from . import conjugate
        return self.applyfunc(conjugate)
    def adjoint(self): return self.H
    def is_diagonal(self): return all(i == j for i,j in self._data)
    def is_symmetric(self, simplify=True): return self == self.T
    def trace(self):
        if not self.is_square: raise NonSquareMatrixError('Matrix must be square')
        return sum((v for (i,j),v in self._data.items() if i == j),_ZERO)
    def row(self, i): return self[_index(i,self.rows),:]
    def col(self, j): return self[:,_index(j,self.cols)]
    def extract(self, rows, cols):
        rows,cols=list(rows),list(cols)
        ri,cj=[_index(i,self.rows) for i in rows],[_index(j,self.cols) for j in cols]
        return self._new(len(ri),len(cj),{(i,j):self._data[r,c] for i,r in enumerate(ri) for j,c in enumerate(cj) if (r,c) in self._data})
    def reshape(self, rows, cols):
        rows,cols=_dimension(rows),_dimension(cols)
        if rows*cols != len(self): raise ValueError('Invalid reshape size')
        return self._new(rows,cols,{divmod(i*self.cols+j,cols):v for (i,j),v in self._data.items()})
    def multiply_elementwise(self, other):
        if not _matrix_like(other) or self.shape != (int(other.rows),int(other.cols)): raise ShapeError('Matrix size mismatch')
        return self._new(self.rows,self.cols,{k:v*other[k] for k,v in self._data.items() if v*other[k] != 0},other)
    def applyfunc(self, function):
        zero = sympify(function(_ZERO))
        if zero != 0:
            return self._new(self.rows,self.cols,{(i,j):function(self[i,j]) for i in range(self.rows) for j in range(self.cols)})
        return self._new(self.rows,self.cols,{k:result for k,v in self._data.items() if (result := sympify(function(v))) != 0})
    def subs(self, *args, **kwargs): return self.applyfunc(lambda x:x.subs(*args,**kwargs))
    def xreplace(self, rule): return self.applyfunc(lambda x:x.xreplace(rule))
    def diff(self, *args, **kwargs): return self.applyfunc(lambda x:x.diff(*args,**kwargs))
    def expand(self, **kwargs): return self.applyfunc(lambda x:x.expand(**kwargs))
    def simplify(self, **kwargs): return self.applyfunc(lambda x:x.simplify(**kwargs))
    def row_swap(self, i, j):
        self._mutable(); i,j=_index(i,self.rows),_index(j,self.rows)
        self._data = {(j if r==i else i if r==j else r,c):v for (r,c),v in self._data.items()}
    def col_swap(self, i, j):
        self._mutable(); i,j=_index(i,self.cols),_index(j,self.cols)
        self._data = {(r,j if c==i else i if c==j else c):v for (r,c),v in self._data.items()}
    def row_del(self, i):
        self._mutable(); i=_index(i,self.rows)
        self._data = {(r-(r>i),c):v for (r,c),v in self._data.items() if r != i}; self._rows -= 1
    def col_del(self, j):
        self._mutable(); j=_index(j,self.cols)
        self._data = {(r,c-(c>j)):v for (r,c),v in self._data.items() if c != j}; self._cols -= 1
    def row_join(self, other): return self.hstack(self,other)
    def col_join(self, other): return self.vstack(self,other)

    @classmethod
    def zeros(cls, rows, cols=None): return cls(_dimension(rows),_dimension(rows if cols is None else cols),{})
    @classmethod
    def ones(cls, rows, cols=None):
        rows,cols=_dimension(rows),_dimension(rows if cols is None else cols)
        return cls(rows,cols,{(i,j):_ONE for i in range(rows) for j in range(cols)})
    @classmethod
    def eye(cls, rows, cols=None):
        rows,cols=_dimension(rows),_dimension(rows if cols is None else cols)
        return cls(rows,cols,{(i,i):_ONE for i in range(min(rows,cols))})
    @classmethod
    def diag(cls, *entries):
        data={}; rows=cols=0
        for value in entries:
            block = cls(value) if _matrix_like(value) or isinstance(value,(list,tuple)) else cls(1,1,{(0,0):value})
            data.update({(i+rows,j+cols):v for (i,j),v in block._data.items()})
            rows += block.rows; cols += block.cols
        return cls(rows,cols,data)
    @classmethod
    def hstack(cls, *matrices):
        matrices=[m if isinstance(m,_SparseBase) else cls(m) for m in matrices]
        if not matrices: return cls()
        rows=matrices[0].rows; cols=0; data={}
        for matrix in matrices:
            if matrix.rows != rows: raise ShapeError('Matrices must have equal rows')
            data.update({(i,j+cols):v for (i,j),v in matrix._data.items()}); cols += matrix.cols
        target = ImmutableSparseMatrix if cls._immutable or any(m._immutable for m in matrices) else cls
        return target(rows,cols,data)
    @classmethod
    def vstack(cls, *matrices):
        return cls.hstack(*(m.T if isinstance(m,_SparseBase) else cls(m).T for m in matrices)).T

    _dense_methods = frozenset('det inv inverse_GE inverse_LU inverse_ADJ inverse_CH inverse_LDL inverse_QR adjugate minor minor_submatrix cofactor cofactor_matrix rank rref echelon_form nullspace columnspace rowspace rank_decomposition LUdecomposition LUsolve solve gauss_jordan_solve pinv pinv_solve lower_triangular_solve upper_triangular_solve diagonal_solve QRsolve cholesky_solve LDLsolve solve_least_squares charpoly eigenvals eigenvects jacobian norm cholesky LDLdecomposition QRdecomposition diagonalize is_diagonalizable singular_values condition_number jordan_form jordan_cells vec dot cross normalized evalf'.split())
    _bounded_dense_functions = frozenset(('exp', 'log', 'analytic_func', 'singular_value_decomposition', 'inv_mod', 'lll', 'lll_transform'))
    def __getattr__(self, name):
        if name not in self._dense_methods and name not in self._bounded_dense_functions:
            raise AttributeError(name)
        def method(*args,**kwargs):
            if name == 'inv_mod':
                Matrix._inv_mod_preflight(self.rows,self.cols,*args,**kwargs)
            if name in ('lll', 'lll_transform'):
                Matrix._lll_preflight(self.rows,self.cols,*args,**kwargs)
            if name in self._bounded_dense_functions and name not in ('lll', 'lll_transform') and (
                    self.rows*self.cols > 100_000 or self.rows+self.cols > 200_000):
                raise NotImplementedError('Exact matrix function resource limit exceeded before sparse densification')
            result = getattr(self.as_dense(),name)(*args,**kwargs)
            def convert(value):
                if isinstance(value, Matrix): return self._new(value.rows,value.cols,_entries(value),value)
                if isinstance(value, tuple): return tuple(convert(v) for v in value)
                if isinstance(value, list): return [convert(v) for v in value]
                return value
            return convert(result)
        return method


class MutableSparseMatrix(_SparseBase):
    """Mutable sparse matrix storing only nonzero Expr entries."""
    __slots__ = ()


class ImmutableSparseMatrix(_SparseBase):
    """Immutable sparse matrix with protected, hashable sparse storage."""
    __slots__ = ()
    _immutable = True
    is_Immutable = True

    def __hash__(self):
        from ._core import _sparse_matrix_hash
        return _sparse_matrix_hash(self.rows,self.cols,self.todok())


SparseMatrix = MutableSparseMatrix
