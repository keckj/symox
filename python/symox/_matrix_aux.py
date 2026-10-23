"""SymPy ``MatrixBase`` protocol methods and legacy aliases.

The algebraic kernels (Berkowitz vectors, fraction-free LU, Householder
reductions, permanents, Cramer's rule, Gram-Schmidt, ``equals``) are native;
this module only adapts SymPy's signatures, result classes and pure container
protocols (``todok``, ``key2bounds``, ``copyin_*``, ``table``...). Methods are
installed on the native dense class, the immutable facade and the sparse
classes by ``install``.
"""
import random as _random

from ._core import Matrix, sympify

_ALIGN = {'left': 'ljust', 'right': 'rjust', 'center': 'center',
          '<': 'ljust', '>': 'rjust', '^': 'center'}


def _errors():
    from . import NonSquareMatrixError, ShapeError
    return NonSquareMatrixError, ShapeError


def _is_sparse(value):
    return getattr(value, '_symox_sparse', False)


def _a2idx(j, n=None):
    if type(j) is not int:
        jindex = getattr(j, '__index__', None)
        if jindex is not None:
            j = jindex()
        else:
            raise IndexError("Invalid index a[%r]" % (j,))
    if n is not None:
        if j < 0:
            j += n
        if not (0 <= j < n):
            raise IndexError("Index out of range: a[%s]" % (j,))
    return int(j)


def _is_sequence(value):
    from .utilities.iterables import is_sequence
    return is_sequence(value)


def _convert(cls, value):
    """Return a native dense result in the class ``cls`` asked for."""
    if cls is Matrix or cls is None:
        return value
    return cls(value)


# --- container protocol -------------------------------------------------

def iter_items(self):
    for i in range(self.rows):
        for j in range(self.cols):
            value = self[i, j]
            if value != 0:
                yield (i, j), value


def iter_values(self):
    return (value for value in self if value != 0)


def todok(self):
    return {(i, j): v for (i, j), v in iter_items(self)}


def todod(self):
    rows = {}
    for (i, j), v in iter_items(self):
        rows.setdefault(i, {})[j] = v
    return rows


def from_dok(cls, rows, cols, dok):
    if _is_sparse(cls):
        return cls(rows, cols, dok)
    data = [sympify(0)] * (rows * cols)
    for (i, j), value in dok.items():
        data[i * cols + j] = sympify(value)
    return _convert(cls, Matrix(rows, cols, data))


def key2bounds(self, keys):
    islice, jslice = [isinstance(k, slice) for k in keys]
    if islice:
        rlo, rhi = (0, 0) if not self.rows else keys[0].indices(self.rows)[:2]
    else:
        rlo = _a2idx(keys[0], self.rows)
        rhi = rlo + 1
    if jslice:
        clo, chi = (0, 0) if not self.cols else keys[1].indices(self.cols)[:2]
    else:
        clo = _a2idx(keys[1], self.cols)
        chi = clo + 1
    return rlo, rhi, clo, chi


def key2ij(self, key):
    if _is_sequence(key):
        if not len(key) == 2:
            raise TypeError('key must be a sequence of length 2')
        return [_a2idx(i, n) if not isinstance(i, slice) else i
                for i, n in zip(key, self.shape)]
    if isinstance(key, slice):
        return key.indices(len(self))[:2]
    return divmod(_a2idx(key, len(self)), self.cols)


def copyin_matrix(self, key, value):
    _, ShapeError = _errors()
    rlo, rhi, clo, chi = self.key2bounds(key)
    if value.shape != (rhi - rlo, chi - clo):
        raise ShapeError("The Matrix `value` doesn't have the same dimensions "
                         "as the in sub-Matrix given by `key`.")
    for i in range(value.rows):
        for j in range(value.cols):
            self[i + rlo, j + clo] = value[i, j]


def _block_assignment(self, key, value):
    """SymPy's ``_setitem`` rule for a scalar key and a matrix value.

    ``M[i, j] = B`` (or a flat ``M[k] = B``) copies ``B`` with its top-left
    corner at ``(i, j)``; a sequence that is not an expression becomes a
    Matrix first. Returns the ``(slices, matrix)`` to copy in, or None.
    """
    i, j = self.key2ij(key)
    if isinstance(i, slice) or isinstance(j, slice):
        return None
    if not getattr(value, 'is_Matrix', False):
        from ._core import Expr
        if isinstance(value, Expr) or not _is_sequence(value):
            return None
        value = Matrix(value)
    return (slice(i, i + value.rows), slice(j, j + value.cols)), value


_native_setitem = None
_native_has = None


def _matrix_has(self, *patterns):
    # SymPy's ``_eval_has``: any entry has any pattern. The native test takes
    # expressions only, so a pattern class (``A.has(Float)``) asks each entry.
    if any(isinstance(p, type) for p in patterns):
        return any(entry.has(*patterns) for entry in self)
    return _native_has(self, *patterns)


def _matrix_setitem(self, key, value):
    try:
        return _native_setitem(self, key, value)
    except ValueError as exc:
        if str(exc) != 'size mismatch in assignment':
            raise
        block = _block_assignment(self, key, value)
        if block is None:
            raise
    self.copyin_matrix(*block)


def copyin_list(self, key, value):
    if not _is_sequence(value):
        raise TypeError("`value` must be an ordered iterable, not %s." % type(value))
    return self.copyin_matrix(key, type(self)(value))


def fill(self, value):
    value = sympify(value)
    for i in range(self.rows):
        for j in range(self.cols):
            self[i, j] = value


def row_add(self, s, t, k):
    for j in range(self.cols):
        self[t, j] += k * self[s, j]


def row_mult(self, i, factor):
    for j in range(self.cols):
        self[i, j] *= factor


def table(self, printer, rowstart='[', rowend=']', rowsep='\n', colsep=', ', align='right'):
    if 0 in self.shape:
        return '[]'
    res = []
    maxlen = [0] * self.cols
    for i in range(self.rows):
        res.append([])
        for j in range(self.cols):
            s = printer._print(self[i, j])
            res[-1].append(s)
            maxlen[j] = max(len(s), maxlen[j])
    align = _ALIGN[align]
    for i, row in enumerate(res):
        for j, elem in enumerate(row):
            row[j] = getattr(elem, align)(maxlen[j])
        res[i] = rowstart + colsep.join(row) + rowend
    return rowsep.join(res)


def print_nonzero(self, symb="X"):
    lines = []
    for i in range(self.rows):
        lines.append("[%s]" % ''.join(" " if self[i, j] == 0 else str(symb) for j in range(self.cols)))
    print('\n'.join(lines))


# --- entrywise protocol ---------------------------------------------------

def atoms(self, *types):
    values = self.values()
    found = set(sympify(0).atoms(*types)) if len(values) < self.rows * self.cols else set()
    for value in values:
        found |= set(value.atoms(*types))
    return found


def doit(self, **hints):
    return self.applyfunc(lambda x: x.doit(**hints))


def integrate(self, *args, **kwargs):
    return self.applyfunc(lambda x: x.doit().integrate(*args, **kwargs))


def limit(self, *args):
    return self.applyfunc(lambda x: x.doit().limit(*args))


def replace(self, F, G, map=False, simultaneous=True, exact=None):
    kwargs = {}
    if not simultaneous:
        kwargs['simultaneous'] = simultaneous
    if exact is not None:
        kwargs['exact'] = exact
    if map:
        found = {}

        def entry(value):
            value, mapping = value.replace(F, G, map=True, **kwargs)
            found.update(mapping)
            return value
        return self.applyfunc(entry), found
    return self.applyfunc(lambda x: x.replace(F, G, **kwargs))


def refine(self, assumptions=True):
    import importlib
    function = getattr(importlib.import_module(__package__), 'refine', None)
    if function is None:
        raise NotImplementedError('refine is not available in symox')
    return self.applyfunc(lambda x: function(x, assumptions))


def as_real_imag(self, deep=True, **hints):
    from . import re, im
    return self.applyfunc(re), self.applyfunc(im)


# --- algebra and legacy aliases -------------------------------------------

def add(self, b):
    return self + b


def rmultiply(self, other, dotprodsimp=None):
    _, ShapeError = _errors()
    if hasattr(other, 'shape') and hasattr(other, 'rows'):
        if self.shape[0] != other.shape[1]:
            raise ShapeError("Matrix size mismatch.")
        return other * self
    return other * self


def pow(self, exp, method=None):
    NonSquareMatrixError, _ = _errors()
    if method is not None and method not in ['multiply', 'mulsimp', 'jordan', 'cayley']:
        raise TypeError('No such method')
    if self.rows != self.cols:
        raise NonSquareMatrixError()
    exp = sympify(exp)
    integral = exp.is_Number and exp == int(exp)
    if method in ('multiply', 'mulsimp', 'cayley') and not integral:
        raise ValueError("%s method is only valid for integer powers" % method)
    if not integral and not exp.is_zero and exp != 1 and self.is_diagonal():
        return _diag_power(self, exp)
    if not integral:
        identity = self.eye(self.rows)
        if (self - identity).is_nilpotent():
            from . import Dummy
            variable = Dummy('matrix_power')
            return self.analytic_func(variable**exp, variable)
    return self ** exp


def _diag_power(self, exp):
    out = Matrix.zeros(self.rows, self.cols)
    for i in range(self.rows):
        out[i, i] = self[i, i] ** exp
    return _like(self, _convert(type(self), out))


def _like(self, value):
    if getattr(self, '_immutable', False) and hasattr(value, 'as_immutable'):
        return value.as_immutable()
    return value


def vech(self, diagonal=True, check_symmetry=True):
    NonSquareMatrixError, _ = _errors()
    if not self.is_square:
        raise NonSquareMatrixError
    if check_symmetry and not self.is_symmetric():
        raise ValueError("The matrix is not symmetric.")
    return self._vech(diagonal)


def project(self, v):
    return v * (self.dot(v) / v.dot(v))


def LUdecompositionFF(self):
    from ._sparse import MutableSparseMatrix
    dense = self.as_dense() if _is_sparse(self) else self
    P, L, DD, U = dense._lu_ff()
    U = MutableSparseMatrix(U) if _is_sparse(self) else U
    return MutableSparseMatrix(P), MutableSparseMatrix(L), MutableSparseMatrix(DD), U


def berkowitz_det(self):
    return self.det(method='berkowitz')


def det_bareis(self):
    return self.det(method='bareiss')


def det_LU_decomposition(self):
    return self.det(method='lu')


def berkowitz_charpoly(self, x=None, simplify=None):
    if x is None:
        from . import Dummy
        x = Dummy('lambda')
    return self.charpoly(x)


def berkowitz_eigenvals(self, **flags):
    return self.eigenvals(**flags)


def cofactorMatrix(self, method="berkowitz"):
    return self.cofactor_matrix(method=method)


def minorEntry(self, i, j, method="berkowitz"):
    return self.minor(i, j, method=method)


def minorMatrix(self, i, j):
    return self.minor_submatrix(i, j)


def left_eigenvects(self, **flags):
    return [(val, mult, [b.transpose() for b in basis])
            for val, mult, basis in self.transpose().eigenvects(**flags)]


def jordan_cell(self, eigenval, n):
    return _like(self, type(self).jordan_block(size=n, eigenvalue=eigenval))


# --- class constructors ----------------------------------------------------

def _as_int(value):
    import operator
    try:
        return operator.index(value)
    except TypeError:
        if getattr(value, 'is_Integer', False):
            return int(value)
        raise ValueError('%s is not an integer' % (value,)) from None


def jordan_block(cls, size=None, eigenvalue=None, *, band='upper', **kwargs):
    klass = kwargs.pop('cls', cls)
    eigenval = kwargs.get('eigenval', None)
    if eigenvalue is None and eigenval is None:
        raise ValueError("Must supply an eigenvalue")
    elif eigenvalue != eigenval and None not in (eigenval, eigenvalue):
        raise ValueError("Inconsistent values are given: 'eigenval'={}, "
                         "'eigenvalue'={}".format(eigenval, eigenvalue))
    elif eigenval is not None:
        eigenvalue = eigenval
    if size is None:
        raise ValueError("Must supply a matrix size")
    size = _as_int(size)
    return _convert(klass, Matrix._jordan_block(size, eigenvalue, band == 'lower'))


def companion(cls, poly):
    from . import Poly
    poly = sympify(poly) if not isinstance(poly, Poly) else poly
    if not isinstance(poly, Poly):
        raise ValueError("{} must be a Poly instance.".format(poly))
    if not poly.is_monic:
        raise ValueError("{} must be a monic polynomial.".format(poly))
    if not poly.is_univariate:
        raise ValueError("{} must be a univariate polynomial.".format(poly))
    if not poly.degree() >= 1:
        raise ValueError("{} must have degree not less than 1.".format(poly))
    return _convert(cls, Matrix._companion(poly.all_coeffs()))


def wilkinson(cls, n, **kwargs):
    klass = kwargs.get('cls', cls)
    minus, plus = Matrix._wilkinson(_as_int(n))
    return _convert(klass, minus), _convert(klass, plus)


def irregular(cls, ntop, *matrices, **kwargs):
    ntop = _as_int(ntop)
    explicit = [m.as_explicit() if hasattr(m, 'as_explicit') else m for m in matrices]
    explicit = [m.as_dense() if _is_sparse(m) else m for m in explicit]
    return _convert(cls, Matrix._irregular(ntop, explicit))


def orthogonalize(cls, *vecs, normalize=False, rankcheck=False):
    dense = [v.as_dense() if _is_sparse(v) else v for v in vecs]
    return [_convert(cls, q) for q in Matrix._orthogonalize(dense, normalize, rankcheck)]


def __abs__(self):
    """Returns a new matrix with entry-wise absolute values."""
    from . import Abs
    return self.applyfunc(Abs)


def __mod__(self, other):
    return self.applyfunc(lambda x: x % other)


_INSTANCE = ('iter_items iter_values todok todod key2bounds key2ij table print_nonzero atoms doit '
             'integrate limit replace refine as_real_imag add rmultiply pow vech project LUdecompositionFF '
             'berkowitz_det det_bareis det_LU_decomposition berkowitz_charpoly berkowitz_eigenvals '
             'cofactorMatrix minorEntry minorMatrix left_eigenvects jordan_cell __abs__ __mod__').split()
_MUTATORS = 'copyin_matrix copyin_list fill row_add row_mult'.split()
_CLASS = 'from_dok jordan_block companion wilkinson irregular orthogonalize'.split()
# Native dense kernels reached by sparse matrices through densification.
SPARSE_DENSE = frozenset((
    'lower_triangular upper_triangular rot90 _vech hat vee dual per berkowitz berkowitz_minors '
    'LUdecomposition_Simple is_nilpotent is_symbolic get_diag_blocks connected_components inverse_BLOCK '
    'upper_hessenberg_decomposition bidiagonal_decomposition bidiagonalize rref_rhs is_anti_symmetric '
    'cramer_solve strongly_connected_components col_insert row_insert diagonal flat has n').split())


# --- sparse-only helpers -----------------------------------------------------

def liupc(self):
    # Algorithm 2.4 of Van Grondelle (1999), ported from SymPy.
    R = [[] for _ in range(self.rows)]
    for r, c, _ in self.row_list():
        if c <= r:
            R[r].append(c)
    inf = len(R)
    parent = [inf] * self.rows
    virtual = [inf] * self.rows
    for r in range(self.rows):
        for c in R[r][:-1]:
            while virtual[c] < r:
                t = virtual[c]
                virtual[c] = r
                c = t
            if virtual[c] == inf:
                parent[c] = virtual[c] = r
    return R, parent


def row_structure_symbolic_cholesky(self):
    R, parent = self.liupc()
    inf = len(R)
    Lrow = [list(r) for r in R]
    for k in range(self.rows):
        for j in R[k]:
            while j != inf and j != k:
                Lrow[k].append(j)
                j = parent[j]
        Lrow[k] = sorted(set(Lrow[k]))
    return Lrow


def _sparse_dual(self):
    return self.as_dense().dual()


def _sparse_lu_simple(self, iszerofunc=None, simpfunc=None, rankcheck=False):
    from ._sparse import MutableSparseMatrix
    lu, swaps = self.as_dense().LUdecomposition_Simple(iszerofunc, simpfunc, rankcheck)
    cls = type(self) if not len(self) else MutableSparseMatrix
    return cls(lu), swaps


def _sparse_rref_rhs(self, rhs):
    from ._sparse import _entries
    rational = all(value.is_Rational for value in self._data.values()) and all(sympify(value).is_Rational for value in rhs)
    dense_rhs = rhs.as_dense() if _is_sparse(rhs) else rhs
    result = self.as_dense().rref_rhs(dense_rhs)
    if rational:
        return result
    return tuple(self._new(value.rows, value.cols, _entries(value), value) for value in result)


def scalar_multiply(self, scalar):
    return scalar * self


def multiply(self, other, dotprodsimp=None):
    return self * other


_SPARSE_ONLY = 'liupc row_structure_symbolic_cholesky scalar_multiply multiply'.split()


def _immutable_constructor(name):
    def constructor(*args, **kwargs):
        result = getattr(Matrix, name)(*args, **kwargs)
        if isinstance(result, tuple):
            return tuple(r.as_immutable() for r in result)
        if isinstance(result, list):
            return [r.as_immutable() for r in result]
        return result.as_immutable()
    constructor.__name__ = name
    return staticmethod(constructor)


def _sparse_equals(self, other, failing_expression=False):
    if getattr(other, 'shape', None) != self.shape:
        return False
    other = other.as_dense() if _is_sparse(other) else other
    return self.as_dense().equals(other, failing_expression)


class _MutableMethod:
    def __init__(self, function):
        self.function = function

    def __get__(self, instance, owner):
        if instance is None:
            return self.function
        if instance._immutable:
            raise AttributeError(self.function.__name__)
        return self.function.__get__(instance, owner)


def _getitem_iter(self):
    # SymPy matrices define no __iter__: iteration uses __getitem__ with flat
    # indices, so an overridden __getitem__ also controls iteration.
    for index in range(len(self)):
        yield self[index]


def _matrix_init_subclass(cls, **kwargs):
    super(Matrix, cls).__init_subclass__(**kwargs)
    if '__getitem__' in cls.__dict__ and '__iter__' not in cls.__dict__:
        cls.__iter__ = _getitem_iter


_native_eigenvects = None


def _eigenvects(self, error_when_incomplete=True, iszerofunc=None, *, chop=False, **flags):
    """SymPy's ``_eigenvects`` flag handling around the native kernel.

    An explicit ``simplify`` flag also requests primitive eigenvectors,
    ``v/gcd(v)``, exactly as in the reference implementation.
    """
    from types import FunctionType
    from . import gcd, simplify as _simplify
    primitive = flags.pop('simplify', False)
    flags.pop('multiple', None)
    if iszerofunc is not None:
        flags['iszerofunc'] = iszerofunc
    if chop:
        flags['chop'] = chop
    ret = _native_eigenvects(self, error_when_incomplete, **flags)
    if primitive:
        simpfunc = primitive if isinstance(primitive, FunctionType) else _simplify
        ret = [(val, mult, [(v/gcd(list(v))).applyfunc(simpfunc) for v in es])
               for val, mult, es in ret]
    return ret


def _diagonalize(self, reals_only=False, sort=False, normalize=False):
    """Port of SymPy's ``_diagonalize`` over primitive eigenvectors."""
    from . import NonSquareMatrixError, MatrixError
    from ._compat import default_sort_key
    if not self.is_square:
        raise NonSquareMatrixError()
    eigenvecs = self.eigenvects(simplify=True)
    for val, mult, basis in eigenvecs:
        if (reals_only and not val.is_real) or mult != len(basis):
            raise MatrixError('Matrix is not diagonalizable')
    if sort:
        eigenvecs = sorted(eigenvecs, key=default_sort_key)
    p_cols, diag = [], []
    for val, mult, basis in eigenvecs:
        diag += [val]*mult
        p_cols += basis
    if normalize:
        p_cols = [v/v.norm() for v in p_cols]
    P, D = self.hstack(*p_cols), self.diag(*diag)
    if getattr(self, '_immutable', False):
        P, D = P.as_immutable(), D.as_immutable()
    return P, D


def install(immutable_facade, sparse_base, mutable_sparse):
    namespace = globals()
    global _native_eigenvects
    if _native_eigenvects is None:
        _native_eigenvects = Matrix.eigenvects
        Matrix.eigenvects = _eigenvects
        Matrix.diagonalize = _diagonalize
    Matrix.__init_subclass__ = classmethod(_matrix_init_subclass)
    Matrix.__array__ = _matrix_array
    Matrix.__array_priority__ = 11
    global _native_setitem
    if _native_setitem is None:
        _native_setitem = Matrix.__setitem__
        Matrix.__setitem__ = _matrix_setitem
    global _native_has
    if _native_has is None:
        _native_has = Matrix.has
        Matrix.has = _matrix_has
    sparse_base.__array__ = _matrix_array
    for name in _INSTANCE:
        if not hasattr(Matrix, name):
            setattr(Matrix, name, namespace[name])
    for name in _MUTATORS:
        setattr(Matrix, name, _MutableMethod(namespace[name]))
    for name in _CLASS:
        setattr(Matrix, name, classmethod(namespace[name]))
        setattr(immutable_facade, name, _immutable_constructor(name))
        setattr(sparse_base, name, classmethod(namespace[name]))
    for name in _INSTANCE:
        if name not in sparse_base.__dict__ and name not in ('iter_items', 'iter_values', 'todok'):
            setattr(sparse_base, name, namespace[name])
    for name in _MUTATORS:
        setattr(mutable_sparse, name, namespace[name])
    for name in _SPARSE_ONLY:
        setattr(sparse_base, name, namespace[name])
    sparse_base.dual = _sparse_dual
    sparse_base.LUdecomposition_Simple = _sparse_lu_simple
    sparse_base.equals = _sparse_equals
    sparse_base.rref_rhs = _sparse_rref_rhs
    sparse_base._dense_methods = sparse_base._dense_methods | SPARSE_DENSE


# --- top-level functions ---------------------------------------------------

def GramSchmidt(vlist, orthonormal=False):
    return Matrix.orthogonalize(*vlist, normalize=orthonormal, rankcheck=True)


def casoratian(seqs, n, zero=True):
    return Matrix._casoratian([sympify(s) for s in seqs], sympify(n), zero)


def jordan_cell_function(eigenval, n):
    return Matrix.jordan_block(size=n, eigenvalue=eigenval)


def matrix_multiply_elementwise(A, B):
    return A.multiply_elementwise(B)


def list2numpy(l, dtype=object):
    from numpy import empty
    a = empty(len(l), dtype)
    for i, s in enumerate(l):
        a[i] = s
    return a


def _matrix_array(self, dtype=object, copy=None):
    """Expose explicit matrix shape and exact entries to NumPy consumers."""
    if copy is not None and not copy:
        raise TypeError("Cannot implement copy=False when converting Matrix to ndarray")
    return matrix2numpy(self, dtype=dtype)


def matrix2numpy(m, dtype=object):
    from numpy import empty
    a = empty(m.shape, dtype)
    for i in range(m.rows):
        for j in range(m.cols):
            a[i, j] = m[i, j]
    return a


def symarray(prefix, shape, **kwargs):
    from numpy import empty, ndindex
    from . import Symbol
    arr = empty(shape, dtype=object)
    for index in ndindex(shape):
        arr[index] = Symbol('%s_%s' % (prefix, '_'.join(map(str, index))), **kwargs)
    return arr


def randMatrix(r, c=None, min=0, max=99, seed=None, symmetric=False, percent=100, prng=None):
    prng = prng or _random.Random(seed)
    if c is None:
        c = r
    if symmetric and r != c:
        raise ValueError('For symmetric matrices, r must equal c, but %i != %i' % (r, c))
    ij = range(r * c)
    if percent != 100:
        ij = prng.sample(ij, int(len(ij) * percent // 100))
    m = Matrix.zeros(r, c)
    if not symmetric:
        for ijk in ij:
            i, j = divmod(ijk, c)
            m[i, j] = prng.randint(min, max)
    else:
        for ijk in ij:
            i, j = divmod(ijk, c)
            if i <= j:
                m[i, j] = m[j, i] = prng.randint(min, max)
    return m


def rot_givens(i, j, theta, dim=3):
    from . import cos, sin
    if not isinstance(dim, int) or dim < 2:
        raise ValueError('dim must be an integer biggen than one, got {}.'.format(dim))
    if i == j:
        raise ValueError('i and j must be different, got ({}, {})'.format(i, j))
    for ij in [i, j]:
        if not isinstance(ij, int) or ij < 0 or ij > dim - 1:
            raise ValueError('i and j must be integers between 0 and '
                             '{}, got i={} and j={}.'.format(dim - 1, i, j))
    theta = sympify(theta)
    c, s = cos(theta), sin(theta)
    M = Matrix.eye(dim)
    M[i, i] = c
    M[j, j] = c
    M[i, j] = s
    M[j, i] = -s
    return M


def rot_axis1(theta):
    return rot_givens(1, 2, theta, dim=3)


def rot_axis2(theta):
    return rot_givens(2, 0, theta, dim=3)


def rot_axis3(theta):
    return rot_givens(0, 1, theta, dim=3)


def rot_ccw_axis1(theta):
    return rot_givens(2, 1, theta, dim=3)


def rot_ccw_axis2(theta):
    return rot_givens(0, 2, theta, dim=3)


def rot_ccw_axis3(theta):
    return rot_givens(1, 0, theta, dim=3)


def banded(*args, **kwargs):
    """SymPy's ``banded``: a sparse matrix from ``{diagonal: values}``."""
    from ._sparse import SparseMatrix
    from . import Dummy
    try:
        if len(args) not in (1, 2, 3):
            raise TypeError
        if not isinstance(args[-1], dict):
            raise TypeError
        if len(args) == 1:
            rows = kwargs.get('rows', None)
            cols = kwargs.get('cols', None)
            rows = _as_int(rows) if rows is not None else None
            cols = _as_int(cols) if cols is not None else None
        elif len(args) == 2:
            rows = cols = _as_int(args[0])
        else:
            rows, cols = map(_as_int, args[:2])
        _ = all(_as_int(k) for k in args[-1])
    except (ValueError, TypeError):
        raise TypeError('unrecognized input to banded: expecting [[row,] col,] {int: value}')

    def rc(d):
        r = -d if d < 0 else 0
        return r, (0 if r else d)

    smat = {}
    undone = []
    tba = Dummy()
    for d, v in args[-1].items():
        r, c = rc(d)
        if isinstance(v, (list, tuple)):
            extra = 0
            for i, vi in enumerate(v):
                i += extra
                if _is_sequence(vi) or hasattr(vi, 'shape'):
                    vi = SparseMatrix(vi)
                    smat[r + i, c + i] = vi
                    extra += min(vi.shape) - 1
                else:
                    smat[r + i, c + i] = vi
        elif _is_sequence(v) or hasattr(v, 'shape'):
            v = SparseMatrix(v)
            rv, cv = v.shape
            if rows and cols:
                nr, xr = divmod(rows - r, rv)
                nc, xc = divmod(cols - c, cv)
                x = xr or xc
                do = min(nr, nc)
            elif rows:
                do, x = divmod(rows - r, rv)
            elif cols:
                do, x = divmod(cols - c, cv)
            else:
                do, x = 1, 0
            if x:
                raise ValueError('\nsequence does not fit an integral number of times in the matrix')
            j = min(v.shape)
            for i in range(do):
                smat[r, c] = v
                r += j
                c += j
        elif v:
            smat[r, c] = tba
            undone.append((d, v))
    expanded, srows, scols = {}, 0, 0

    def expand(i, j, v):
        # SymPy's SparseMatrix(None, smat) refuses overlapping blocks.
        if v is tba or v:
            if (i, j) in expanded and v != expanded[i, j]:
                raise ValueError('There is a collision at {} for {} and {}.'.format((i, j), v, expanded[i, j]))
            expanded[i, j] = v
    for (i, j), v in smat.items():
        if hasattr(v, 'shape') and hasattr(v, 'todok'):
            for (a, b), w in v.todok().items():
                expand(i + a, j + b, w)
            srows, scols = builtin_max(srows, i + v.rows), builtin_max(scols, j + v.cols)
        else:
            expand(i, j, sympify(v) if v is not tba else v)
            srows, scols = builtin_max(srows, i + 1), builtin_max(scols, j + 1)
    smat = {k: v for k, v in expanded.items() if v is tba or v != 0}
    if rows is not None and rows < srows:
        raise ValueError('Designated rows %s < needed %s' % (rows, srows))
    if cols is not None and cols < scols:
        raise ValueError('Designated cols %s < needed %s' % (cols, scols))
    if rows is None and cols is None:
        rows, cols = srows, scols
    elif rows is not None and cols is None:
        cols = builtin_max(rows, scols)
    elif cols is not None and rows is None:
        rows = builtin_max(cols, srows)

    def update(i, j, v):
        if v:
            if (i, j) in smat and smat[i, j] is not tba and smat[i, j] != v:
                raise ValueError('collision at %s' % ((i, j),))
            smat[i, j] = sympify(v)
    for d, vi in undone:
        r, c = rc(d)
        value = vi if callable(vi) else (lambda _: vi)
        i = 0
        while r + i < rows and c + i < cols:
            update(r + i, c + i, value(i))
            i += 1
    return SparseMatrix(rows, cols, {k: v for k, v in smat.items() if v is not tba})


builtin_max = max


class _MatrixBaseMeta(type):
    def __getattr__(cls,name):
        from ._constructor_methods import instance_method
        return instance_method(Matrix,name)

    def __instancecheck__(cls, obj):
        return isinstance(obj, Matrix) or _is_sparse(obj)

    def __subclasscheck__(cls, sub):
        return (isinstance(sub, type) and issubclass(sub, Matrix)) or getattr(sub, '_symox_sparse', False) is True


class MatrixBase(metaclass=_MatrixBaseMeta):
    """Abstract base of explicit matrices, for ``isinstance`` checks."""

    @classmethod
    def _new(cls, *args, **kwargs):
        raise NotImplementedError('Subclasses must implement this.')

    @classmethod
    def hstack(cls, *args):
        from functools import reduce
        if not args:
            return cls._new()
        return reduce(type(args[0]).row_join, args)

    @classmethod
    def vstack(cls, *args):
        from functools import reduce
        if not args:
            return cls._new()
        return reduce(type(args[0]).col_join, args)
