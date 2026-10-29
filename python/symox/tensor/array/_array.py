"""Owned dense/sparse arrays following SymPy's N-dim array protocol.

Arrays are built in ``__new__`` through the class method ``_new`` (dense) or
directly in ``__new__`` (sparse), exactly like SymPy, so user subclasses that
override ``__new__`` or ``__init__`` keep working. Operations that SymPy
implements with ``type(self)(...)`` do the same here, so they preserve the
subclass. Immutable arrays are structural ``Basic`` objects; mutable arrays
are not.
"""
import itertools
import math
import operator
import sys
from collections.abc import Iterable, Mapping
from types import MappingProxyType

from ... import _core
from ..._basic import Basic

MAX_ENTRIES = 100_000
MAX_RANK = 32


def _shape_tuple(shape):
    if isinstance(shape, _core.Expr) and getattr(shape.func, '__name__', '') == 'Tuple':
        shape = shape.args
    try:
        shape = (operator.index(shape),)
    except TypeError:
        pass
    try:
        shape = tuple(operator.index(n) for n in shape)
    except TypeError:
        raise TypeError('Shape should contain integers only.') from None
    if len(shape) > MAX_RANK:
        raise ValueError('Array rank exceeds 32')
    if any(n < 0 for n in shape):
        raise ValueError('Array dimensions must be nonnegative')
    if math.prod(shape) > sys.maxsize:
        raise ValueError('Array size exceeds native indexing capacity')
    return shape


def _native(value):
    return isinstance(value, (_core.Expr, _core.AlgebraicNumber))


def _check_entry(value):
    if isinstance(value, _core.Expr) and getattr(value.func, '__name__', '') == 'Tuple':
        raise TypeError('A tuple is array structure, not a scalar entry')
    return value


def _scalar(value):
    """A native scalar expression; anything else is rejected."""
    from ... import sympify
    if not _native(value):
        value = sympify(value)
    if not _native(value):
        raise TypeError('Array entries must be native scalar expressions')
    return _check_entry(value)


def _dense_entry(value):
    # SymPy stores dense entries in a Tuple, which sympifies non-strictly:
    # undefined function classes, None and other sympified objects that are
    # not scalar expressions remain valid entries.
    from ... import sympify
    if _native(value):
        return _check_entry(value)
    return _check_entry(sympify(value))


def _sparse_entry(value):
    # SymPy's sparse arrays strictly sympify every stored entry.
    from ... import sympify, SympifyError
    if not _native(value):
        value = sympify(value)
        if not _native(value):
            raise SympifyError(value)
    return _check_entry(value)


def _is_zero(value):
    from ... import S
    try:
        return bool(value == S.Zero)
    except TypeError:
        return False


def _unravel(index, shape):
    coordinates = []
    for dimension in reversed(shape):
        index, position = divmod(index, dimension)
        coordinates.append(position)
    return tuple(reversed(coordinates))


def _ravel(coordinates, shape):
    index = 0
    for position, dimension in zip(coordinates, shape):
        index = index*dimension+position
    return index


def _structured(value):
    """Whether ``value`` is array structure rather than a scalar entry."""
    if isinstance(value, NDimArray) or isinstance(value, _core._MatrixExpr):
        return True
    if hasattr(value, 'shape') and hasattr(value, 'tolist'):
        return True
    if isinstance(value, _core.Expr):
        return getattr(value.func, '__name__', '') == 'Tuple'
    if isinstance(value, (str, bytes, Mapping)):
        return False
    try:
        iter(value)
    except TypeError:
        return False
    return True


def _walk_input(value, put, depth=0, active=None):
    """Put the raw leaves of a nested input in row-major order; return its shape."""
    if depth > MAX_RANK:
        raise ValueError('Array input nesting exceeds 32')
    if active is None:
        active = set()
    if isinstance(value, NDimArray):
        for entry in value._flat_values():
            put(entry)
        return value.shape
    if isinstance(value, _core._MatrixExpr):
        value = value.as_explicit()
    if hasattr(value, 'shape') and hasattr(value, 'tolist'):
        shape = _shape_tuple(value.shape)
        if not shape:
            # A NumPy scalar array: its single entry.
            put(value.tolist() if not hasattr(value, 'item') else value[()])
            return ()
        _walk_input(value.tolist(), put, depth, active)
        return shape
    if isinstance(value, _core.Expr):
        sequence = value.args if getattr(value.func, '__name__', '') == 'Tuple' else None
    elif isinstance(value, (str, bytes, Mapping)):
        sequence = None
    else:
        try:
            sequence = iter(value)
        except TypeError:
            sequence = None
    if sequence is None:
        put(value)
        return ()
    identity = id(value)
    if identity in active:
        raise ValueError('Cyclic array input')
    active.add(identity)
    child_shape, count = None, 0
    try:
        for child in sequence:
            shape = _walk_input(child, put, depth+1, active)
            if child_shape is not None and shape != child_shape:
                raise ValueError('could not determine shape unambiguously')
            child_shape, count = shape, count+1
    finally:
        active.remove(identity)
    return (count,) + (child_shape or ())


def _array_class(sparse=False, mutable=False):
    if sparse:
        return MutableSparseNDimArray if mutable else ImmutableSparseNDimArray
    return MutableDenseNDimArray if mutable else ImmutableDenseNDimArray


def _flat_key(key, shape, size):
    if isinstance(key, (tuple, list)) or (isinstance(key, _core.Expr) and
                                          getattr(key.func, '__name__', '') == 'Tuple'):
        key = tuple(key)
        if len(key) != len(shape):
            raise ValueError('Coordinate rank does not match the array')
        positions = tuple(NDimArray._index(i, n) for i, n in zip(key, shape))
        return _ravel(positions, shape)
    return NDimArray._index(key, size)


class NDimArray:
    """N-dimensional array.

    ``NDimArray(...)`` builds an ``ImmutableDenseNDimArray``, as in SymPy.
    """
    __array_priority__ = 1000
    _sparse = False
    _mutable = False
    _diff_wrt = True
    is_scalar = False
    is_NDimArray = True
    is_Matrix = is_Number = is_Symbol = is_Add = is_Mul = is_Pow = is_Function = False
    is_Atom = is_number = False
    is_commutative = None

    def __new__(cls, iterable, shape=None, **kwargs):
        return ImmutableDenseNDimArray(iterable, shape, **kwargs)

    @property
    def _symox_(self):
        # An array is never a native scalar: hide Basic's conversion hook so
        # scalar operators return NotImplemented and reach the array methods.
        raise AttributeError('_symox_')

    # Construction protocol -------------------------------------------------

    @classmethod
    def _scan_iterable_shape(cls, iterable):
        values = []
        shape = _walk_input(iterable, values.append)
        return values, shape

    @classmethod
    def _handle_ndarray_creation_inputs(cls, iterable=None, shape=None, **kwargs):
        if shape is None:
            if iterable is None:
                return (), ()
            if isinstance(iterable, SparseNDimArray):
                return iterable._shape, dict(iterable._sparse_array)
            if isinstance(iterable, NDimArray):
                shape = iterable.shape
            elif isinstance(iterable, Mapping):
                raise ValueError('A sparse mapping needs an explicit shape')
            elif _structured(iterable):
                if hasattr(iterable, 'shape') and not isinstance(iterable, _core._MatrixExpr):
                    shape = iterable.shape
                else:
                    iterable, shape = cls._scan_iterable_shape(iterable)
            else:
                return (), (iterable,)
        shape = _shape_tuple(shape)
        if isinstance(iterable, Mapping):
            size = math.prod(shape)
            iterable = {_flat_key(key, shape, size): value for key, value in iterable.items()}
        return shape, iterable

    @classmethod
    def _check_special_bounds(cls, flat_list, shape):
        if shape == () and len(flat_list) != 1:
            raise ValueError('arrays without shape need one scalar value')
        if shape == (0,) and len(flat_list) > 0:
            raise ValueError('if array shape is (0,) there cannot be elements')

    @classmethod
    def _collect(cls, iterable, shape):
        """Validated ``(shape, storage)`` for ``cls``: a list (dense) or a dict (sparse)."""
        sparse = cls._sparse
        entry = _sparse_entry if sparse else _dense_entry
        if iterable is None and shape is None:
            iterable, shape = (), ()
        shape = _shape_tuple(shape) if shape is not None else None
        if hasattr(iterable, 'shape'):
            source_shape = _shape_tuple(iterable.shape)
            if shape is None:
                shape = source_shape
            elif math.prod(shape) != math.prod(source_shape):
                raise ValueError('Shape does not match the number of entries')
        size = math.prod(shape) if shape is not None else None
        if not sparse and size is not None and size > MAX_ENTRIES:
            raise ValueError('Dense arrays are limited to 100000 entries')
        if isinstance(iterable, NDimArray):
            shape = iterable.shape if shape is None else shape
            size = math.prod(shape)
            if size != len(iterable):
                raise ValueError('Shape does not match the number of entries')
            if not sparse and size > MAX_ENTRIES:
                raise ValueError('Dense arrays are limited to 100000 entries')
            if sparse:
                return shape, {index: entry(value) for index, value in iterable._items()}
            from ... import S
            values = [S.Zero]*size
            for index, value in iterable._items():
                values[index] = value
            return shape, values
        if isinstance(iterable, Mapping) or (sparse and getattr(iterable, '_symox_sparse', False)):
            if not isinstance(iterable, Mapping):
                original_shape = _shape_tuple(iterable.shape)
                iterable = {_ravel(position, original_shape): value
                            for position, value in iterable.todok().items()}
            if shape is None:
                raise ValueError('A sparse mapping needs an explicit shape')
            size = math.prod(shape)
            if sparse:
                values = {}
            else:
                from ... import S
                values = [S.Zero]*size
            for key, value in iterable.items():
                index = _flat_key(key, shape, size)
                value = entry(value)
                if sparse:
                    if not _is_zero(value):
                        values[index] = value
                    else:
                        values.pop(index, None)
                else:
                    values[index] = value
            if sparse and len(values) > MAX_ENTRIES:
                raise ValueError('Sparse arrays are limited to 100000 stored entries')
            return shape, values
        values = {} if sparse else []
        count = 0

        def put(value):
            nonlocal count
            if size is not None and count >= size:
                raise ValueError('Too many entries for array shape')
            value = entry(value)
            if sparse:
                if not _is_zero(value):
                    values[count] = value
                    if len(values) > MAX_ENTRIES:
                        raise ValueError('Sparse arrays are limited to 100000 stored entries')
            else:
                if count >= MAX_ENTRIES:
                    raise ValueError('Dense arrays are limited to 100000 entries')
                values.append(value)
            count += 1
        inferred = _walk_input(() if iterable is None else iterable, put)
        shape = _shape_tuple(inferred) if shape is None else shape
        if count != math.prod(shape):
            raise ValueError('Shape does not match the number of entries')
        return shape, values

    def _set_storage(self, shape, values):
        from ... import Integer
        size = math.prod(shape)
        self._shape = shape
        self._rank = len(shape)
        self._size = size
        self._loop_size = size if self._mutable else Integer(size)
        if self._sparse:
            self._sparse_array = values if self._mutable else MappingProxyType(values)
        else:
            self._array = values
        self._hash_cache = None
        return self

    @staticmethod
    def _index(value, dimension):
        value = operator.index(value)
        if not -dimension <= value < dimension:
            raise ValueError('Array index out of range')
        return value % dimension

    # SymPy's index helpers --------------------------------------------------

    def _parse_index(self, index):
        try:
            index = operator.index(index)
        except TypeError:
            pass
        else:
            if index >= len(self):
                raise ValueError('Only a tuple index is accepted')
            return index
        if len(self) == 0:
            raise ValueError('Index not valid with an empty array')
        if len(index) != self._rank:
            raise ValueError('Wrong number of array axes')
        real_index = 0
        for i, dimension in zip(index, self.shape):
            i = operator.index(i)
            if i >= dimension or i < -dimension:
                raise ValueError('Index ' + str(index) + ' out of border')
            real_index = real_index*dimension + i % dimension
        return real_index

    def _get_tuple_index(self, integer_index):
        index = []
        for dimension in reversed(self.shape):
            index.append(integer_index % dimension)
            integer_index //= dimension
        index.reverse()
        return tuple(index)

    def _check_index_for_getitem(self, index):
        if isinstance(index, slice):
            index = (index,)
        else:
            try:
                index = (operator.index(index),)
            except TypeError:
                pass
        if len(index) < self.rank():
            index = tuple(index) + tuple(slice(None) for i in range(len(index), self.rank()))
        if len(index) > self.rank():
            raise ValueError('Dimension of index greater than rank of array')
        return index

    def _slice_expand(self, s, dim):
        if not isinstance(s, slice):
            return (s,)
        start, stop, step = s.indices(dim)
        return list(range(start, stop, step))

    def _get_slice_data_for_array_access(self, index):
        sl_factors = [self._slice_expand(i, dim) for (i, dim) in zip(index, self.shape)]
        eindices = itertools.product(*sl_factors)
        return sl_factors, eindices

    def _get_slice_data_for_array_assignment(self, index, value):
        if not isinstance(value, NDimArray):
            value = type(self)(value)
        sl_factors, eindices = self._get_slice_data_for_array_access(index)
        slice_offsets = [min(i) if isinstance(i, list) else None for i in sl_factors]
        return value, eindices, slice_offsets

    @classmethod
    def _setter_iterable_check(cls, value):
        if isinstance(value, (NDimArray, _core.Matrix, _core._MatrixExpr)) or (
                not isinstance(value, _core.Expr) and isinstance(value, Iterable)):
            raise NotImplementedError

    # Structure ----------------------------------------------------------------

    @property
    def _data(self):
        return self._sparse_array if self._sparse else self._array

    @property
    def shape(self):
        return self._shape

    @property
    def kind(self):
        from ..._kind import scalar_kind, NumberKind, UndefinedKind
        from .kind import ArrayKind
        kinds = set()
        stored = 0
        for _, value in self._items():
            kinds.add(scalar_kind(value) if _native(value) else
                      getattr(value, 'kind', UndefinedKind))
            stored += 1
        if stored < len(self):
            kinds.add(NumberKind)
        return ArrayKind(next(iter(kinds)) if len(kinds) == 1 else UndefinedKind)

    @property
    def func(self):
        return type(self)

    @property
    def args(self):
        from ... import Dict, Tuple
        self._materialize_limit()
        shape = Tuple(*self.shape)
        try:
            if self._sparse:
                return Dict(dict(self._sparse_array)), shape
            return Tuple(*self._array), shape
        except (TypeError, ValueError):
            # Entries that native containers cannot hold (an undefined
            # function class, None) keep a plain tuple of arguments.
            return tuple(self._flat_values()), shape

    @property
    def free_symbols(self):
        found = set()
        for _, value in self._items():
            symbols = getattr(value, 'free_symbols', None)
            if isinstance(symbols, (set, frozenset)):
                found |= symbols
        return found

    def rank(self):
        return self._rank

    def __len__(self):
        return self._size

    def _items(self):
        if self._sparse:
            return iter(self._sparse_array.items())
        return ((i, value) for i, value in enumerate(self._array) if not _is_zero(value))

    def _get(self, index):
        if self._sparse:
            from ... import S
            return self._sparse_array.get(index, S.Zero)
        return self._array[index]

    def _flat_values(self):
        if not self._sparse:
            return iter(self._array)
        return (self._get(i) for i in range(len(self)))

    def _materialize_limit(self, nested=False):
        positions = len(self)
        if nested:
            product, positions = 1, 0
            for dimension in self.shape:
                product *= dimension
                positions += product
        if positions > MAX_ENTRIES:
            raise ValueError('Array materialization exceeds 100000 entries or containers')

    def _selection(self, key):
        from ... import sympify
        key = key if isinstance(key, tuple) else (key,)
        if len(key) > self.rank():
            raise ValueError('Dimension of index greater than rank of array')
        key += (slice(None),)*(self.rank()-len(key))
        coordinates, shape = [], []
        for index, dimension in zip(key, self.shape):
            if isinstance(index, slice):
                selected = range(*index.indices(dimension))
                coordinates.append(selected)
                shape.append(len(selected))
            else:
                try:
                    selected = self._index(index, dimension)
                except TypeError:
                    selected = sympify(index)
                    if not isinstance(selected, _core.Expr) or selected.is_integer is False:
                        raise ValueError('Array indices must be integers') from None
                    if selected.is_number:
                        raise ValueError('Array index is not a concrete integer')
                coordinates.append((selected,))
        return coordinates, tuple(shape)

    def __getitem__(self, key):
        coordinates, shape = self._selection(key)
        if not shape:
            return self._entry(tuple(axis[0] for axis in coordinates))
        if 0 in shape:
            return type(self)({}, shape)
        if math.prod(shape) > MAX_ENTRIES:
            if self._sparse and all(isinstance(axis, range) or isinstance(axis[0], int) for axis in coordinates):
                data = {}
                for index, value in self._items():
                    old = _unravel(index, self.shape)
                    if all(i in axis for i, axis in zip(old, coordinates)):
                        new = tuple((i-axis.start)//axis.step for i, axis in zip(old, coordinates)
                                    if isinstance(axis, range))
                        data[_ravel(new, shape)] = value
                return type(self)(data, shape)
            raise ValueError('Array slice exceeds 100000 entries')
        return type(self)([self._entry(indices) for indices in itertools.product(*coordinates)], shape)

    def _entry(self, coordinates):
        if all(isinstance(i, int) for i in coordinates):
            return self._get(_ravel(coordinates, self.shape))
        from ... import S, Piecewise, Eq, And, Contains
        normalized, bounds = [], []
        for index, dimension in zip(coordinates, self.shape):
            if isinstance(index, int):
                normalized.append(index)
                continue
            value = (index if index.is_nonnegative is True else index+dimension
                     if index.is_negative is True else Piecewise((index+dimension, index < 0), (index, True)))
            normalized.append(value)
            bounds.extend((value >= 0, value < dimension))
            if index.is_integer is not True:
                bounds.append(Contains(value, S.Integers))
        cases = []
        for index, value in self._items():
            position = _unravel(index, self.shape)
            if any(isinstance(i, int) and i != j for i, j in zip(normalized, position)):
                continue
            condition = And(*(Eq(i, j) for i, j in zip(normalized, position) if not isinstance(i, int)))
            cases.append((value, condition))
            if len(cases) > 4096:
                raise NotImplementedError('Symbolic array indexing exceeds 4096 cases')
        cases.extend(((S.Zero, And(*bounds)), (S.NaN, True)))
        return Piecewise(*cases)

    def __iter__(self):
        if not self.shape:
            yield self._get(0)
        else:
            for i in range(self.shape[0]):
                yield self[i]

    def tolist(self):
        if not self.shape:
            # SymPy's nested-list conversion indexes the first axis.
            raise IndexError('tuple index out of range')
        self._materialize_limit(nested=True)

        def build(prefix, axis):
            if axis == self.rank():
                return self._get(_ravel(prefix, self.shape))
            return [build(prefix+(i,), axis+1) for i in range(self.shape[axis])]
        return build((), 0)

    def __repr__(self):
        if not self.shape:
            return str(self._get(0))
        if 0 in self.shape:
            return '%s([], %r)' % (type(self).__name__, self.shape)
        if self._sparse and len(self) > MAX_ENTRIES:
            return '%s(%r, %r)' % (type(self).__name__, dict(self._sparse_array), self.shape)
        self._materialize_limit(nested=True)

        def build(prefix, axis):
            if axis == self.rank():
                return str(self._get(_ravel(prefix, self.shape)))
            return '[' + ', '.join(build(prefix+(i,), axis+1) for i in range(self.shape[axis])) + ']'
        return build((), 0)

    __str__ = __repr__

    def __eq__(self, other):
        return (isinstance(other, NDimArray) and self.shape == other.shape
                and dict(self._items()) == dict(other._items()))

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        if self._hash_cache is None:
            self._hash_cache = hash((self.shape, tuple(sorted(self._items(), key=operator.itemgetter(0)))))
        return self._hash_cache

    def sort_key(self, order=None):
        from ... import S, Tuple, default_sort_key
        if order not in (None, 'lex'):
            raise NotImplementedError('Array sort keys support lexicographic order')
        one = S.One.sort_key()

        def basic_key(name, children):
            return (5, 0, name), (len(children), tuple(children)), one, S.One
        if self._mutable:
            self._materialize_limit(nested=True)
            children = tuple(default_sort_key(child) for child in self)
            return (10, 0, type(self).__name__), (len(children), children), one, S.One
        if self._sparse:
            pairs = tuple(Tuple(index, value).sort_key() for index, value in sorted(self._items()))
            storage = basic_key('Dict', pairs)
        else:
            storage = Tuple(*self._array).sort_key()
        return basic_key(type(self).__name__, (storage, Tuple(*self.shape).sort_key()))

    def __reduce__(self):
        values = dict(self._sparse_array) if self._sparse else list(self._array)
        return type(self), (values, self.shape)

    @classmethod
    def zeros(cls, *shape):
        from ... import S
        if len(shape) == 1 and isinstance(shape[0], (tuple, list)):
            shape = shape[0]
        shape = _shape_tuple(shape)
        if cls._sparse:
            return cls({}, shape)
        if math.prod(shape) > MAX_ENTRIES:
            raise ValueError('Dense arrays are limited to 100000 entries')
        return cls._new([S.Zero]*math.prod(shape), shape)

    @classmethod
    def ones(cls, *shape):
        from ... import S
        if len(shape) == 1 and isinstance(shape[0], (tuple, list)):
            shape = shape[0]
        shape = _shape_tuple(shape)
        if math.prod(shape) > MAX_ENTRIES:
            raise ValueError('An all-ones array exceeds 100000 stored entries')
        return cls([S.One]*math.prod(shape), shape)

    def reshape(self, *shape):
        if len(shape) == 1 and isinstance(shape[0], (tuple, list)):
            shape = shape[0]
        shape = _shape_tuple(shape)
        if math.prod(shape) != len(self):
            raise ValueError('Expecting reshape size to %d but got prod(%s) = %d'
                             % (len(self), str(shape), math.prod(shape)))
        return type(self)(dict(self._sparse_array) if self._sparse else self._array, shape)

    def copy(self):
        return type(self)(self, self.shape)

    __copy__ = copy

    def __deepcopy__(self, memo):
        return self.copy()

    def as_immutable(self):
        return _array_class(self._sparse)(self) if self._mutable else self

    def as_mutable(self):
        return _array_class(self._sparse, True)(self)

    # Arithmetic ---------------------------------------------------------------

    def _binary(self, other, subtract=False):
        from ... import S
        if not isinstance(other, NDimArray):
            return NotImplemented
        if self.shape != other.shape:
            raise ValueError('array shape mismatch')
        if not self._sparse:
            if subtract:
                values = [a-b for a, b in zip(self._array, other._flat_values())]
            else:
                values = [a+b for a, b in zip(self._array, other._flat_values())]
            return type(self)(values, self.shape)
        left, right = dict(self._items()), dict(other._items())
        data = {i: left.get(i, S.Zero)-right.get(i, S.Zero) if subtract else
                left.get(i, S.Zero)+right.get(i, S.Zero) for i in left.keys() | right.keys()}
        return type(self)(data, self.shape)

    def __add__(self, other):
        return self._binary(other)

    def __sub__(self, other):
        return self._binary(other, True)

    def __neg__(self):
        return self.applyfunc(operator.neg)

    def __pos__(self):
        return self.copy()

    @staticmethod
    def _scalar_operand(value):
        if isinstance(value, (NDimArray, _core.Matrix, _core._MatrixExpr)) or (
                not isinstance(value, _core.Expr) and getattr(value, '_symox_sparse', False)) or (
                not _native(value) and isinstance(value, Iterable)):
            raise ValueError('scalar expected, use tensorproduct(...) for tensorial product')
        return _scalar(value)

    def __mul__(self, scalar):
        try:
            scalar = self._scalar_operand(scalar)
        except TypeError:
            return NotImplemented
        return self.applyfunc(lambda entry: entry*scalar)

    def __rmul__(self, scalar):
        try:
            scalar = self._scalar_operand(scalar)
        except TypeError:
            return NotImplemented
        return self.applyfunc(lambda entry: scalar*entry)

    def __truediv__(self, scalar):
        if isinstance(scalar, (NDimArray, _core.Matrix, _core._MatrixExpr)) or (
                not _native(scalar) and isinstance(scalar, Iterable)):
            raise ValueError('scalar expected')
        scalar = _scalar(scalar)
        return self.applyfunc(lambda entry: entry/scalar)

    def __rtruediv__(self, other):
        raise NotImplementedError('unsupported operation on NDimArray')

    def applyfunc(self, f):
        from ... import S
        if not callable(f):
            raise TypeError('applyfunc requires a callable')
        if not len(self):
            return type(self)({}, self.shape)
        if self._sparse:
            zero = _scalar(f(S.Zero))
            if zero == S.Zero:
                return type(self)({i: f(v) for i, v in self._items()}, self.shape)
            self._materialize_limit()
        return type(self)([f(v) for v in self._flat_values()], self.shape)

    def diff(self, *variables, **options):
        if options.get('evaluate', True) is False:
            raise NotImplementedError('Held array derivatives require lazy array expressions')
        if variables and options.get('evaluate', True) is True:
            # Symbolic orders by scalar Symbols: SymPy's held ArrayDerivative.
            from .array_derivatives import symbolic_order_derivative
            held = symbolic_order_derivative(self, variables)
            if held is not None:
                return held
        from .arrayop import derive_by_array, _validate_variable
        from ._adapters import _array_variable, _check_derivative_options
        if not variables:
            symbols = self.free_symbols
            if len(symbols) != 1:
                raise ValueError('Supply differentiation variables explicitly')
            variables = (next(iter(symbols)),)
        if len(variables) == 1 and _array_variable(variables[0]):
            _check_derivative_options(options)
            return derive_by_array(self, variables[0])
        seen_variable = False
        for variable in variables:
            if isinstance(variable, tuple):
                if len(variable) != 2:
                    raise ValueError('Derivative specification must contain a variable and order')
                _validate_variable(variable[0])
                seen_variable = True
            else:
                try:
                    operator.index(variable)
                except TypeError:
                    _validate_variable(variable)
                    seen_variable = True
                else:
                    if not seen_variable:
                        raise ValueError('Derivative order without a variable')
        scalar_variables = [v[0] if isinstance(v, tuple) else v for v in variables]
        # SymPy differentiates ``self.as_immutable()`` entrywise.
        source = self.as_immutable()
        if any(getattr(getattr(v, 'func', None), '__name__', '') == 'MatrixElement'
               and not getattr(v, 'is_UndefinedFunction', False) for v in scalar_variables):
            from .arrayop import _expand_matrix_scalars
            return source.applyfunc(lambda entry: _expand_matrix_scalars(entry).diff(*variables, **options))
        return source.applyfunc(lambda entry: entry.diff(*variables, **options))

    def subs(self, *args, **options):
        return self.applyfunc(lambda entry: entry.subs(*args, **options))

    def xreplace(self, rule):
        for old, new in rule.items():
            if isinstance(old, NDimArray) and self == old:
                return new
        rule = {k: v for k, v in rule.items() if not isinstance(k, NDimArray)}
        return self.applyfunc(lambda entry: entry.xreplace(rule))

    def has(self, *patterns):
        if any(isinstance(p, NDimArray) and p == self for p in patterns):
            return True
        return any(value.has(*patterns) for _, value in self._items() if hasattr(value, 'has'))

    def evalf(self, n=15, **options):
        return self.applyfunc(lambda entry: entry.evalf(n, **options))

    n = evalf

    def simplify(self, **options):
        return self.applyfunc(lambda entry: entry.simplify(**options))

    def _eval_simplify(self, **options):
        from ... import simplify
        return self.applyfunc(lambda entry: simplify(entry, **options))

    def doit(self, **options):
        return self.applyfunc(lambda entry: entry.doit(**options))

    def conjugate(self):
        from ... import conjugate
        return self.applyfunc(conjugate)

    def transpose(self):
        from .arrayop import permutedims
        if self.rank() != 2:
            raise ValueError('array rank not 2')
        return permutedims(self, (1, 0))

    T = property(transpose)

    def adjoint(self):
        return self.transpose().conjugate()

    H = property(adjoint)

    def tomatrix(self):
        from ... import Matrix, SparseMatrix
        if self.rank() != 2:
            raise ValueError('Dimensions must be of size of 2')
        if self._sparse:
            return SparseMatrix(*self.shape, {_unravel(i, self.shape): v for i, v in self._items()})
        return Matrix(*self.shape, list(self._array))

    def to_sympy(self):
        import sympy
        name = _array_class(self._sparse, self._mutable).__name__
        values = ({i: v.to_sympy() for i, v in self._items()} if self._sparse else
                  [v.to_sympy() for v in self._array])
        if self._sparse and not self.shape:
            # SymPy 1.14 rejects an empty sparse dict for a scalar zero tensor.
            values = [self._get(0).to_sympy()]
        return getattr(sympy, name)(values, self.shape)

    _sympy_ = to_sympy

    def __array__(self, dtype=None, copy=None):
        import numpy
        if copy is False:
            raise ValueError('Native expressions cannot share a NumPy storage buffer')
        self._materialize_limit()
        return numpy.array(list(self._flat_values()), dtype=object if dtype is None else dtype).reshape(self.shape)


class DenseNDimArray(NDimArray):

    def __new__(cls, *args, **kwargs):
        return ImmutableDenseNDimArray(*args, **kwargs)


class SparseNDimArray(NDimArray):
    _sparse = True

    def __new__(cls, *args, **kwargs):
        return ImmutableSparseNDimArray(*args, **kwargs)


class MutableNDimArray(NDimArray):
    _mutable = True
    __hash__ = None

    def as_immutable(self):
        return _array_class(self._sparse)(self)

    def as_mutable(self):
        return self

    def __setitem__(self, key, value):
        coordinates, shape = self._selection(key)
        if any(not isinstance(axis, range) and not isinstance(axis[0], int) for axis in coordinates):
            raise ValueError('Mutable array assignment requires concrete indices')
        entry = _sparse_entry if self._sparse else _dense_entry
        if not shape:
            self._setter_iterable_check(value)
            entries = [(_ravel(tuple(axis[0] for axis in coordinates), self.shape), entry(value))]
        else:
            source = value if isinstance(value, NDimArray) else Array(value)
            if source.shape != shape:
                raise ValueError('Slice assignment shapes must match')
            if 0 in shape:
                return
            if math.prod(shape) > MAX_ENTRIES:
                raise ValueError('Array slice assignment exceeds 100000 entries')
            entries = [(_ravel(index, self.shape), entry(value)) for index, value in
                       zip(itertools.product(*coordinates), source._flat_values())]
        if self._sparse:
            data = self._sparse_array
            count = len(data)
            for index, value in entries:
                count += int(not _is_zero(value))-int(index in data)
            if count > MAX_ENTRIES:
                raise ValueError('Sparse arrays are limited to 100000 stored entries')
            for index, value in entries:
                if _is_zero(value):
                    data.pop(index, None)
                else:
                    data[index] = value
        else:
            for index, value in entries:
                self._array[index] = value


class ImmutableNDimArray(NDimArray, Basic):
    # Basic precedes object but follows NDimArray: the array protocol wins.
    __hash__ = NDimArray.__hash__

    def __setitem__(self, index, value):
        raise TypeError('immutable N-dim array')

    def as_immutable(self):
        return self

    def as_mutable(self):
        return _array_class(self._sparse, True)(self)

    def atoms(self, *types):
        # SymPy's Basic.atoms over the stored arguments: the flat entries
        # (a Dict of the nonzero entries when sparse) and the shape.
        from ... import Dict, Integer, Tuple
        types = tuple(t if isinstance(t, type) else type(t) for t in types)
        found = Tuple(*self.shape).atoms(*types)
        if self._sparse:
            try:
                found |= Dict(dict(self._sparse_array)).atoms(*types)
            except (TypeError, ValueError):
                for key, value in self._sparse_array.items():
                    found |= Integer(key).atoms(*types) | self._leaf_atoms(value, types)
        else:
            try:
                found |= Tuple(*self._array).atoms(*types)
            except (TypeError, ValueError):
                for value in self._array:
                    found |= self._leaf_atoms(value, types)
        if types and isinstance(self, types):
            found.add(self)
        return found

    @staticmethod
    def _leaf_atoms(value, types):
        if _native(value) or isinstance(value, Basic):
            return value.atoms(*types)
        return {value} if types and isinstance(value, types) else set()


class ImmutableDenseNDimArray(DenseNDimArray, ImmutableNDimArray):

    def __new__(cls, iterable, shape=None, **kwargs):
        return cls._new(iterable, shape, **kwargs)

    @classmethod
    def _new(cls, iterable, shape, **kwargs):
        if kwargs:
            raise TypeError('Unknown array options: '+', '.join(sorted(kwargs)))
        shape, values = cls._collect(iterable, shape)
        return object.__new__(cls)._set_storage(shape, values)


class MutableDenseNDimArray(DenseNDimArray, MutableNDimArray):

    def __new__(cls, iterable=None, shape=None, **kwargs):
        return cls._new(iterable, shape, **kwargs)

    @classmethod
    def _new(cls, iterable, shape, **kwargs):
        # SymPy's mutable arrays ignore extra keyword arguments.
        shape, values = cls._collect(iterable, shape)
        return object.__new__(cls)._set_storage(shape, values)


class ImmutableSparseNDimArray(SparseNDimArray, ImmutableNDimArray):

    def __new__(cls, iterable=None, shape=None, **kwargs):
        if kwargs:
            raise TypeError('Unknown array options: '+', '.join(sorted(kwargs)))
        shape, values = cls._collect(iterable, shape)
        return object.__new__(cls)._set_storage(shape, values)


class MutableSparseNDimArray(MutableNDimArray, SparseNDimArray):

    def __new__(cls, iterable=None, shape=None, **kwargs):
        shape, values = cls._collect(iterable, shape)
        return object.__new__(cls)._set_storage(shape, values)


Array = ImmutableDenseNDimArray
