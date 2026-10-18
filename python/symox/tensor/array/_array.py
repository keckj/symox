"""Owned dense/sparse arrays of native scalar expressions."""
import itertools
import math
import operator
import sys
from collections.abc import Mapping
from types import MappingProxyType

from ... import _core

MAX_ENTRIES = 100_000
MAX_RANK = 32


def _shape_tuple(shape):
    if isinstance(shape, _core.Expr) and getattr(shape.func, '__name__', '') == 'Tuple':
        shape = shape.args
    try:
        shape = tuple(operator.index(n) for n in shape)
    except TypeError:
        raise TypeError('Array dimensions must be concrete integers') from None
    if len(shape) > MAX_RANK:
        raise ValueError('Array rank exceeds 32')
    if any(n < 0 for n in shape):
        raise ValueError('Array dimensions must be nonnegative')
    if math.prod(shape) > sys.maxsize:
        raise ValueError('Array size exceeds native indexing capacity')
    return shape


def _scalar(value):
    from ... import sympify
    if not isinstance(value,(_core.Expr,_core.AlgebraicNumber)):
        value = sympify(value)
    if not isinstance(value, (_core.Expr, _core.AlgebraicNumber)):
        raise TypeError('Array entries must be native scalar expressions')
    if isinstance(value, _core.Expr) and getattr(value.func, '__name__', '') == 'Tuple':
        raise TypeError('A tuple is array structure, not a scalar entry')
    return value


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


def _walk_input(value, put, depth=0, active=None):
    if depth > MAX_RANK:
        raise ValueError('Array input nesting exceeds 32')
    if active is None:
        active = set()
    if isinstance(value, NDimArray):
        for entry in value._flat_values():
            put(entry)
        return value.shape
    if isinstance(value,_core._MatrixExpr):
        value=value.as_explicit()
    if hasattr(value, 'shape') and hasattr(value, 'tolist'):
        shape = _shape_tuple(value.shape)
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
        put(_scalar(value))
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
                raise ValueError('Array input must be rectangular')
            child_shape, count = shape, count+1
    finally:
        active.remove(identity)
    return (count,) + (child_shape or ())


def _array_class(sparse=False, mutable=False):
    if sparse:
        return MutableSparseNDimArray if mutable else ImmutableSparseNDimArray
    return MutableDenseNDimArray if mutable else ImmutableDenseNDimArray


class NDimArray:
    __slots__ = ('_shape', '_data', '_hash_cache')
    __array_priority__ = 1000
    _sparse = False
    _mutable = False
    is_NDimArray = True
    is_Matrix = is_Number = is_Symbol = is_Add = is_Mul = is_Pow = is_Function = False
    is_Atom = is_number = False
    is_commutative = None

    def __new__(cls, iterable=None, shape=None, **options):
        if cls in (NDimArray, DenseNDimArray, ImmutableNDimArray):
            cls = ImmutableDenseNDimArray
        elif cls is SparseNDimArray:
            cls = ImmutableSparseNDimArray
        elif cls is MutableNDimArray:
            cls = MutableDenseNDimArray
        return object.__new__(cls)

    def __init__(self, iterable=None, shape=None, **options):
        from ... import S
        if options:
            raise TypeError('Unknown array options: '+', '.join(sorted(options)))
        shape = _shape_tuple(shape) if shape is not None else None
        if hasattr(iterable,'shape'):
            source_shape=_shape_tuple(iterable.shape)
            if shape is None:
                shape=source_shape
            elif math.prod(shape)!=math.prod(source_shape):
                raise ValueError('Shape does not match the number of entries')
        size = math.prod(shape) if shape is not None else None
        if not self._sparse and size is not None and size > MAX_ENTRIES:
            raise ValueError('Dense arrays are limited to 100000 entries')
        values = {} if self._sparse else []
        if isinstance(iterable, NDimArray):
            shape = iterable.shape if shape is None else shape
            size = math.prod(shape)
            if size != len(iterable):
                raise ValueError('Shape does not match the number of entries')
            if not self._sparse and size > MAX_ENTRIES:
                raise ValueError('Dense arrays are limited to 100000 entries')
            values = dict(iterable._items()) if self._sparse else [S.Zero]*size
            if not self._sparse:
                for index, value in iterable._items():
                    values[index] = value
        elif isinstance(iterable, Mapping) or (
                self._sparse and getattr(iterable, '_symox_sparse', False)):
            if not isinstance(iterable, Mapping):
                original_shape=_shape_tuple(iterable.shape)
                iterable={_ravel(position,original_shape):value for position,value in iterable.todok().items()}
            if shape is None:
                raise ValueError('A sparse mapping needs an explicit shape')
            size = math.prod(shape)
            if not self._sparse:
                values = [S.Zero]*size
            for key, value in iterable.items():
                if isinstance(key, tuple):
                    if len(key) != len(shape):
                        raise ValueError('Coordinate rank does not match the array')
                    positions = tuple(self._index(i,n) for i,n in zip(key,shape))
                    index = _ravel(positions,shape)
                else:
                    index = self._index(key,size)
                value = _scalar(value)
                if self._sparse:
                    if value != S.Zero:
                        values[index] = value
                    else:
                        values.pop(index,None)
                else:
                    values[index] = value
            if self._sparse and len(values) > MAX_ENTRIES:
                raise ValueError('Sparse arrays are limited to 100000 stored entries')
        else:
            count = 0
            def put(value):
                nonlocal count
                if size is not None and count >= size:
                    raise ValueError('Too many entries for array shape')
                if self._sparse:
                    if value != S.Zero:
                        values[count] = value
                        if len(values) > MAX_ENTRIES:
                            raise ValueError('Sparse arrays are limited to 100000 stored entries')
                else:
                    if count >= MAX_ENTRIES:
                        raise ValueError('Dense arrays are limited to 100000 entries')
                    values.append(value)
                count += 1
            inferred = _walk_input([] if iterable is None else iterable, put)
            shape = _shape_tuple(inferred) if shape is None else shape
            if count != math.prod(shape):
                raise ValueError('Shape does not match the number of entries')
        self._shape = shape
        self._data = (values if self._mutable else MappingProxyType(values)) if self._sparse else (
            values if self._mutable else tuple(values))
        self._hash_cache = None

    @staticmethod
    def _index(value, dimension):
        value = operator.index(value)
        if not -dimension <= value < dimension:
            raise ValueError('Array index out of range')
        return value % dimension

    @property
    def shape(self):
        return self._shape

    @property
    def func(self):
        return type(self)

    @property
    def args(self):
        self._materialize_limit()
        return tuple(self._flat_values()), self.shape

    @property
    def free_symbols(self):
        return set().union(*(value.free_symbols for _,value in self._items()))

    def rank(self):
        return len(self.shape)

    def __len__(self):
        return math.prod(self.shape)

    def _items(self):
        from ... import S
        if self._sparse:
            return iter(self._data.items())
        return ((i,value) for i,value in enumerate(self._data) if value != S.Zero)

    def _get(self, index):
        from ... import S
        return self._data.get(index,S.Zero) if self._sparse else self._data[index]

    def _flat_values(self):
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
        key = key if isinstance(key,tuple) else (key,)
        if len(key) > self.rank():
            raise ValueError('Too many array indices')
        key += (slice(None),)*(self.rank()-len(key))
        coordinates, shape = [], []
        for index, dimension in zip(key,self.shape):
            if isinstance(index,slice):
                selected = range(*index.indices(dimension))
                coordinates.append(selected)
                shape.append(len(selected))
            else:
                try:
                    selected = self._index(index,dimension)
                except TypeError:
                    selected = sympify(index)
                    if not isinstance(selected,_core.Expr) or selected.is_integer is False:
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
            return type(self)({},shape)
        if math.prod(shape) > MAX_ENTRIES:
            if self._sparse and all(isinstance(axis,range) or isinstance(axis[0],int) for axis in coordinates):
                data = {}
                for index,value in self._items():
                    old = _unravel(index,self.shape)
                    if all(i in axis for i,axis in zip(old,coordinates)):
                        new = tuple((i-axis.start)//axis.step for i,axis in zip(old,coordinates)
                                    if isinstance(axis,range))
                        data[_ravel(new,shape)] = value
                return type(self)(data,shape)
            raise ValueError('Array slice exceeds 100000 entries')
        return type(self)((self._entry(indices) for indices in itertools.product(*coordinates)),shape)

    def _entry(self, coordinates):
        if all(isinstance(i,int) for i in coordinates):
            return self._get(_ravel(coordinates,self.shape))
        from ... import S, Piecewise, Eq, And, Contains
        normalized, bounds = [], []
        for index, dimension in zip(coordinates,self.shape):
            if isinstance(index,int):
                normalized.append(index)
                continue
            value = (index if index.is_nonnegative is True else index+dimension
                     if index.is_negative is True else Piecewise((index+dimension,index<0),(index,True)))
            normalized.append(value)
            bounds.extend((value>=0,value<dimension))
            if index.is_integer is not True:
                bounds.append(Contains(value,S.Integers))
        cases = []
        for index,value in self._items():
            position = _unravel(index,self.shape)
            if any(isinstance(i,int) and i!=j for i,j in zip(normalized,position)):
                continue
            condition = And(*(Eq(i,j) for i,j in zip(normalized,position) if not isinstance(i,int)))
            cases.append((value,condition))
            if len(cases) > 4096:
                raise NotImplementedError('Symbolic array indexing exceeds 4096 cases')
        cases.extend(((S.Zero,And(*bounds)),(S.NaN,True)))
        return Piecewise(*cases)

    def __iter__(self):
        if not self.shape:
            yield self._get(0)
        else:
            for i in range(self.shape[0]):
                yield self[i]

    def tolist(self):
        self._materialize_limit(nested=True)
        def build(prefix,axis):
            if axis == self.rank():
                return self._get(_ravel(prefix,self.shape))
            return [build(prefix+(i,),axis+1) for i in range(self.shape[axis])]
        return build((),0)

    def __repr__(self):
        if 0 in self.shape:
            return '%s([], %r)' % (type(self).__name__,self.shape)
        if self._sparse and len(self) > MAX_ENTRIES:
            return '%s(%r, %r)' % (type(self).__name__,dict(self._data),self.shape)
        return str(self.tolist())

    __str__ = __repr__

    def __eq__(self, other):
        return isinstance(other,NDimArray) and self.shape == other.shape and dict(self._items()) == dict(other._items())

    def __hash__(self):
        if self._hash_cache is None:
            self._hash_cache = hash((self.shape,tuple(sorted(self._items()))))
        return self._hash_cache

    def sort_key(self, order=None):
        from ... import S, Tuple, default_sort_key
        if order not in (None,'lex'):
            raise NotImplementedError('Array sort keys support lexicographic order')
        one=S.One.sort_key()
        def basic_key(name,children):
            return (5,0,name),(len(children),tuple(children)),one,S.One
        if self._mutable:
            self._materialize_limit(nested=True)
            children=tuple(default_sort_key(child) for child in self)
            return (10,0,type(self).__name__),(len(children),children),one,S.One
        if self._sparse:
            pairs=tuple(Tuple(index,value).sort_key() for index,value in sorted(self._items()))
            storage=basic_key('Dict',pairs)
        else:
            storage=Tuple(*self._data).sort_key()
        return basic_key(type(self).__name__,(storage,Tuple(*self.shape).sort_key()))

    def __reduce__(self):
        values = dict(self._data) if self._sparse else list(self._data)
        return type(self),(values,self.shape)

    @classmethod
    def zeros(cls, *shape):
        from ... import S
        if len(shape)==1 and isinstance(shape[0],(tuple,list)):
            shape=shape[0]
        shape=_shape_tuple(shape)
        return cls({},shape) if cls._sparse else cls(itertools.repeat(S.Zero,math.prod(shape)),shape)

    @classmethod
    def ones(cls, *shape):
        from ... import S
        if len(shape)==1 and isinstance(shape[0],(tuple,list)):
            shape=shape[0]
        shape=_shape_tuple(shape)
        if math.prod(shape)>MAX_ENTRIES:
            raise ValueError('An all-ones array exceeds 100000 stored entries')
        return cls(itertools.repeat(S.One,math.prod(shape)),shape)

    def reshape(self, *shape):
        if len(shape)==1 and isinstance(shape[0],(tuple,list)):
            shape=shape[0]
        shape=_shape_tuple(shape)
        if math.prod(shape)!=len(self):
            raise ValueError('Reshape must preserve the number of entries')
        return type(self)(self,shape)

    def copy(self):
        return type(self)(self,self.shape)

    __copy__ = copy

    def __deepcopy__(self, memo):
        return self.copy()

    def as_immutable(self):
        return _array_class(self._sparse)(self) if self._mutable else self

    def as_mutable(self):
        return _array_class(self._sparse,True)(self)

    def _binary(self, other, subtract=False):
        from ... import S
        if not isinstance(other,NDimArray):
            return NotImplemented
        if self.shape!=other.shape:
            raise ValueError('Array shapes must match')
        left,right=dict(self._items()),dict(other._items())
        data={i:left.get(i,S.Zero)-right.get(i,S.Zero) if subtract else
              left.get(i,S.Zero)+right.get(i,S.Zero) for i in left.keys()|right.keys()}
        return type(self)(data,self.shape)

    def __add__(self, other):
        return self._binary(other)

    def __sub__(self, other):
        return self._binary(other,True)

    def __neg__(self):
        return self.applyfunc(operator.neg)

    def __pos__(self):
        return self.copy()

    def __mul__(self, scalar):
        if isinstance(scalar,NDimArray):
            raise ValueError('Use tensorproduct for multiplication of arrays')
        try:
            scalar=_scalar(scalar)
        except TypeError:
            return NotImplemented
        return self.applyfunc(lambda entry: entry*scalar)

    def __rmul__(self, scalar):
        try:
            scalar=_scalar(scalar)
        except TypeError:
            return NotImplemented
        return self.applyfunc(lambda entry: scalar*entry)

    def __truediv__(self, scalar):
        scalar=_scalar(scalar)
        return self.applyfunc(lambda entry: entry/scalar)

    def applyfunc(self, function):
        from ... import S
        if not callable(function):
            raise TypeError('applyfunc requires a callable')
        if not len(self):
            return type(self)({},self.shape)
        if self._sparse:
            zero=_scalar(function(S.Zero))
            if zero==S.Zero:
                return type(self)({i:function(v) for i,v in self._items()},self.shape)
            self._materialize_limit()
        return type(self)((function(v) for v in self._flat_values()),self.shape)

    def diff(self, *variables, **options):
        if options.get('evaluate',True) is False:
            raise NotImplementedError('Held array derivatives require lazy array expressions')
        from .arrayop import derive_by_array, _validate_variable
        from ._adapters import _array_variable, _check_derivative_options
        if not variables:
            symbols=self.free_symbols
            if len(symbols)!=1:
                raise ValueError('Supply differentiation variables explicitly')
            variables=(next(iter(symbols)),)
        if len(variables)==1 and _array_variable(variables[0]):
            _check_derivative_options(options)
            return derive_by_array(self,variables[0])
        seen_variable=False
        for variable in variables:
            if isinstance(variable,tuple):
                if len(variable)!=2:
                    raise ValueError('Derivative specification must contain a variable and order')
                _validate_variable(variable[0])
                seen_variable=True
            else:
                try:
                    operator.index(variable)
                except TypeError:
                    _validate_variable(variable)
                    seen_variable=True
                else:
                    if not seen_variable:
                        raise ValueError('Derivative order without a variable')
        scalar_variables=[v[0] if isinstance(v,tuple) else v for v in variables]
        if any(getattr(getattr(v,'func',None),'__name__','')=='MatrixElement'
               and not getattr(v,'is_UndefinedFunction',False) for v in scalar_variables):
            from .arrayop import _expand_matrix_scalars
            result=self.applyfunc(lambda entry: _expand_matrix_scalars(entry).diff(*variables,**options))
        else:
            result=self.applyfunc(lambda entry: entry.diff(*variables,**options))
        return result.as_immutable()

    def subs(self, *args, **options):
        return self.applyfunc(lambda entry: entry.subs(*args,**options))

    def xreplace(self, mapping):
        for old,new in mapping.items():
            if isinstance(old,NDimArray) and self==old:
                return new
        mapping={k:v for k,v in mapping.items() if not isinstance(k,NDimArray)}
        return self.applyfunc(lambda entry: entry.xreplace(mapping))

    def evalf(self, n=15, **options):
        return self.applyfunc(lambda entry: entry.evalf(n,**options))

    n = evalf

    def simplify(self, **options):
        return self.applyfunc(lambda entry: entry.simplify(**options))

    def doit(self, **options):
        return self.applyfunc(lambda entry: entry.doit(**options))

    def conjugate(self):
        from ... import conjugate
        return self.applyfunc(conjugate)

    def transpose(self):
        from .arrayop import permutedims
        if self.rank()!=2:
            raise ValueError('Transpose requires rank two; use permutedims for other ranks')
        return permutedims(self,(1,0))

    T = property(transpose)

    def adjoint(self):
        return self.transpose().conjugate()

    H = property(adjoint)

    def tomatrix(self):
        from ... import Matrix, SparseMatrix
        if self.rank()!=2:
            raise ValueError('Only rank-two arrays convert to matrices')
        if self._sparse:
            return SparseMatrix(*self.shape,{_unravel(i,self.shape):v for i,v in self._items()})
        return Matrix(*self.shape,list(self._data))

    def to_sympy(self):
        import sympy
        name=_array_class(self._sparse,self._mutable).__name__
        values=({i:v.to_sympy() for i,v in self._items()} if self._sparse else
                [v.to_sympy() for v in self._data])
        if self._sparse and not self.shape:
            # SymPy 1.14 rejects an empty sparse dict for a scalar zero tensor.
            values=[self._get(0).to_sympy()]
        return getattr(sympy,name)(values,self.shape)

    _sympy_ = to_sympy

    def __array__(self, dtype=None, copy=None):
        import numpy
        if copy is False:
            raise ValueError('Native expressions cannot share a NumPy storage buffer')
        self._materialize_limit()
        return numpy.array(list(self._flat_values()),dtype=object if dtype is None else dtype).reshape(self.shape)


class DenseNDimArray(NDimArray):
    __slots__ = ()


class SparseNDimArray(NDimArray):
    __slots__ = ()
    _sparse = True


class MutableNDimArray(NDimArray):
    __slots__ = ()
    _mutable = True
    __hash__ = None

    def __setitem__(self, key, value):
        from ... import S
        coordinates, shape=self._selection(key)
        if any(not isinstance(axis,range) and not isinstance(axis[0],int) for axis in coordinates):
            raise ValueError('Mutable array assignment requires concrete indices')
        if not shape:
            entries=[(_ravel(tuple(axis[0] for axis in coordinates),self.shape),_scalar(value))]
        else:
            source=value if isinstance(value,NDimArray) else Array(value)
            if source.shape!=shape:
                raise ValueError('Slice assignment shapes must match')
            if 0 in shape:
                return
            if math.prod(shape)>MAX_ENTRIES:
                raise ValueError('Array slice assignment exceeds 100000 entries')
            entries=[(_ravel(index,self.shape),value) for index,value in
                     zip(itertools.product(*coordinates),source._flat_values())]
        if self._sparse:
            count=len(self._data)
            for index,value in entries:
                count += int(value!=S.Zero)-int(index in self._data)
            if count>MAX_ENTRIES:
                raise ValueError('Sparse arrays are limited to 100000 stored entries')
        for index,value in entries:
            if self._sparse and value==S.Zero:
                self._data.pop(index,None)
            else:
                self._data[index]=value


class ImmutableNDimArray(NDimArray):
    __slots__ = ()


class ImmutableDenseNDimArray(DenseNDimArray,ImmutableNDimArray):
    __slots__ = ()


class MutableDenseNDimArray(DenseNDimArray,MutableNDimArray):
    __slots__ = ()
    _mutable = True
    __hash__ = None


class ImmutableSparseNDimArray(SparseNDimArray,ImmutableNDimArray):
    __slots__ = ()


class MutableSparseNDimArray(SparseNDimArray,MutableNDimArray):
    __slots__ = ()
    _mutable = True
    __hash__ = None


Array = ImmutableDenseNDimArray
